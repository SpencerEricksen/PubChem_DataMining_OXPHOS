#!/home/ssericksen/anaconda2/envs/py36_chem/bin/python3.6

import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_score
# use custom scoring function in CV
from sklearn.metrics import make_scorer
from ast import literal_eval

def convert_str_to_list( input_str ):
    '''pandas stores descriptor list as a literal string, use
       literal_eval to get list object back for parsing, etc.'''
    try:
        l = literal_eval( input_str )
    except:
        l = np.nan
        pass
    return l

def enrichment_factor( labels_arr, scores_arr, percentile ):
    '''calculate the enrichment factor based on some upper fraction
       of library ordered by docking scores. upper fraction is determined
       by percentile (actually a fraction of value 0.0-1.0). Also calculate
       the theoretical maximum EF given the percentile threshold'''
    sample_size = int(labels_arr.shape[0] * percentile)  # determine number mols in subset
    indices = np.argsort(scores_arr)[::-1][:sample_size] # get the index positions for these in library
    n_actives = np.nansum(labels_arr)                    # count number of positive labels in library
    n_actives_sample = np.nansum( labels_arr[indices] )    # count number of positive labels in subset
    if n_actives > 0.0:
       ef = float(n_actives_sample) / n_actives / percentile # calc EF at percentile
       ef_max = min(n_actives, sample_size) / (n_actives * percentile)
    else:
       ef = np.nan
       ef_max = np.nan
    return n_actives, ef, ef_max

def normalized_enrichment_factor_single(labels_arr, scores_arr, percentile):
    n_actives, ef, ef_max = enrichment_factor(labels_arr, scores_arr, percentile)
    return ef/ef_max

df = pd.read_csv('training_data.csv')
# shuffle it up order of rows in dataframe
df = df.sample(frac=1, random_state=42).reset_index( drop=True )

feat_type = 'morgan'
#feat_type = 'descriptors'
#feat_type = 'both'
#max_feat = 'log2'
max_feat = 'sqrt'
n_est = 8000

# features 
if feat_type == 'descriptors':
    df['descriptors'] = df['descriptors'].apply( convert_str_to_list )
    #df['descriptors'] = df['descriptors'].apply( lambda x: np.array(x) )
    df = df.dropna()
    X = np.stack( df['descriptors'].values )
    X = np.delete( X, 32, 1)
    #df_feats = pd.DataFrame( X )
elif feat_type == 'morgan':
    df['morgan'] = df['morgan_bitstring'].apply( lambda x: [int(i) for i in list(x)] )
    df = df.dropna()
    X = np.stack( df['morgan'].values )
elif feat_type == 'both':
    df['descriptors'] = df['descriptors'].apply( convert_str_to_list )
    df['morgan'] = df['morgan_bitstring'].apply( lambda x: [int(i) for i in list(x)] )
    df = df.dropna()
    X1 = np.stack( df['descriptors'].values )
    X1 = np.delete( X1, 32, 1 )
    X2 = np.stack( df['morgan'].values )
    X = np.concatenate( (X1, X2), axis=1 )

# labels
y = df['label'].values

# split data set into test and train set
skf = StratifiedKFold( n_splits=5, random_state=42, shuffle=True )
for train_idx, test_idx in skf.split(X,y):
    print("TRAIN:", train_idx, "TEST:", test_idx)
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

# build and train RF
rfc = RandomForestClassifier( max_features=max_feat,
                              min_samples_leaf=1, 
                              n_estimators=n_est, 
                              oob_score=False, 
                              class_weight='balanced', 
                              n_jobs=18, 
                              verbose=1 )

# this doesn't seem to be returning correct NEF value on 5-fold CV
#nef_scorer = make_scorer( normalized_enrichment_factor_single, percentile=0.01 )

# run cross validation (5-fold split training set)
#scores = cross_validate(rfc,X_train,y_train,scoring=(nef_scorer),cv=5, return_estimator=True)
scores = cross_validate(rfc,X_train,y_train,scoring=('accuracy','roc_auc','f1'),cv=5, return_estimator=True)

# check performance
print('+'*50)
print( 'fit_time', scores['fit_time'], np.mean(scores['fit_time']) )
print( 'score_time', scores['score_time'], np.mean(scores['score_time']) )
#print( scores['train_score'], np.mean(scores['train_score'])
#print( 'NEF test_scores', scores['test_score'], np.mean(scores['test_score']), np.std(scores['test_score']) )  
print( 'test_accuracy:', scores['test_accuracy'], np.mean(scores['test_accuracy']), np.std(scores['test_accuracy']) )
print( 'test_roc_auc:', scores['test_roc_auc'], np.mean(scores['test_roc_auc']), np.std(scores['test_roc_auc']) )
print( 'test_f1:', scores['test_f1'], np.mean(scores['test_f1']), np.std(scores['test_f1']) )
print('+'*50)

# train on full train set (80% of data)
rfc.fit(X_train, y_train)

# save the model to disk
fname = 'rf_'+str(max_feat)+'_'+str(n_est)+'_'+feat_type+'_scramble.sav'
pickle.dump( rfc, open(fname, 'wb'))

# get probability scores from RFC model on test set (prioritize test set of molecules)
probs = rfc.predict_proba( X_test )
prob_scores = probs[:,1]

# build dataframe for test set molecules with true labels and model scores
prob_series = pd.DataFrame(prob_scores)
prob_series.reset_index( inplace=True, drop=True )

# get the molids and labels from datafile
q = df.iloc[ test_idx ][['PUBCHEM_CID','label']]
q.reset_index( inplace=True, drop=True )

# merge molecule scores with cids and labels
df_test_scores = pd.concat([q, prob_series], axis=1 )

# dmp to a CSV
df_test_scores.rename( columns={0:'score'}, inplace=True )
df_test_scores.sort_values( by='score', ascending=False, inplace=True )

# add individual estimators...
probs_0 = scores['estimator'][0].predict_proba(X_test)
probs_1 = scores['estimator'][1].predict_proba(X_test)
probs_2 = scores['estimator'][2].predict_proba(X_test)
probs_3 = scores['estimator'][3].predict_proba(X_test)
probs_4 = scores['estimator'][4].predict_proba(X_test)
probs_0_scores = probs_0[:,1]
probs_1_scores = probs_1[:,1]
probs_2_scores = probs_2[:,1]
probs_3_scores = probs_3[:,1]
probs_4_scores = probs_4[:,1]
prob_df = pd.DataFrame( [probs_0_scores, probs_1_scores, probs_2_scores, probs_3_scores, probs_4_scores ] )
prob_df = prob_df.T
prob_df.reset_index( inplace=True, drop=True )
df_test_scores_2 = pd.concat([ df_test_scores, prob_df ], axis=1 )
df_test_scores_2.sort_values( by='score', ascending=False, inplace=True )
df_test_scores_2.to_csv('rf_'+str(max_feat)+'_'+str(n_est)+'_'+feat_type+'_scores_w_cv_estimators_scramble.csv', index=False )

'''
# to use pickeled model later...
# load model from disk
loaded_model = pickle.load( open(fname, 'rb') )
y_pred = loaded_model.predict( X_test )
#result = loaded_model.score( X_test, y_test )
#print(result)
'''
