
import sys
import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_validate
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler

try:
    # feature types -- 'fps', 'desc' or 'both'
    feat_type = sys.argv[1]
    # type of regularization (norm of the penalty) l1 or l2
    reg = sys.argv[2]
except:
    print("provide feature type (desc or fps) and regularization (l1 or l2) as input arguments")
    sys.exit(1)

# which solver ('lbfgs' for L2 and 'liblinear' for l1)
if reg == 'l1':
    #sr = 'liblinear'
    sr = 'saga'
elif reg == 'l2':
    #sr = 'lbfgs'
    sr = 'saga'

# max number of iterations of optimization
n_iter = 200

# load data, get features 
df = pd.read_pickle('../data/training_data.pkl')
# in descriptors, had to replace infinite values with NaNs
df = df.replace( [np.inf, -np.inf], np.nan )
# get molids
molid_list = df.PUBCHEM_CID.tolist()

# features
col_list = df.columns.tolist()
if feat_type == 'desc':
    feat_cols = [ c for c in col_list[2:] if 'desc' in str(c) ]
    feat_cols.remove('desc_Ipc')
    df_X = df[feat_cols]
    # replace NaNs with column means (impute)
    df_X = df_X.fillna( df_X.mean() )
    # load the scaler previously fit on the training data
    scaler = pickle.load( open('../data/stand_scaler_desc_training_data.pkl', 'rb') )
    X = scaler.transform( df_X )
elif feat_type == 'fps':
    feat_cols = [ c for c in col_list[2:] if 'desc' not in str(c) ]
    X = df[feat_cols].values
elif feat_type == 'both':
    # get dataframe of desc features
    desc_cols = [ c for c in col_list[2:] if 'desc' in str(c) ]
    desc_cols.remove('desc_Ipc')
    df_X1 = df[desc_cols]
    # replace NaNs with column means (impute)
    df_X1 = df_X1.fillna( df_X1.mean() )
    scaler = pickle.load( open('../data/stand_scaler_desc_training_data.pkl', 'rb') )
    X1 = scaler.transform( df_X1 )
    df_X1 = pd.DataFrame( X1, columns=desc_cols, index=molid_list )
    # get data frame of fps features
    fps_cols = [ c for c in col_list[2:] if 'desc' not in str(c) ]
    X2 = df[fps_cols].values
    df_X2 = pd.DataFrame( X2, columns=fps_cols, index=molid_list )
    # merge desc and fps into single dataframe, get feature array
    df_X12 = pd.concat( [df_X1, df_X2], axis=1 )
    X = df_X12.values


# labels
y = df.label.values
#weird issue with int type in y, had to set type to match X
y = np.array( y, np.int64 )

# search hyperparams for LR model (tests 10 different C values from 1E-4 to 1E4.
# C is inverse of regularization strength (small value => strong regularization)
# note that default cross-validation generator is Stratified K-folds (5-folds)
skf = StratifiedKFold( n_splits=5, shuffle=True, random_state=42 )
model = LogisticRegressionCV( Cs=10,
                              fit_intercept=True,
                              cv=skf,
                              penalty=reg,
                              scoring='balanced_accuracy',
                              solver=sr,
                              max_iter=n_iter,
                              tol=1e-4,
                              class_weight='balanced',
                              n_jobs=-1, 
                              random_state=42 ).fit( X, y )
# scoring='average_precision'

# test best hyperpararm configuration in 5-fold CV (default cross-validation generator is 5-fold Stratified K-Folds)
#cv = StratifiedKFold( n_splits=5, random_state=42)
scores = cross_validate( model, X, y, scoring=('balanced_accuracy','average_precision','roc_auc','f1'), cv=5, return_estimator=True )

# print performance in 5-CV
metrics = list(scores.keys())[-4:]
print( "model,"+','.join( metrics ))
model_name = reg + "_" + sr + "_C" + '{:.5f}'.format( model.C_[0] ) + "_" + feat_type
print( model_name+","+",".join( [ '{:.3f}'.format( np.mean(scores[i]) ) for i in metrics ] ))

# train on full train set (80% of total data)
model.fit(X, y)
# save the model to disk
fname = "model_"+model_name+".pkl"
pickle.dump( model, open(fname, 'wb'))

# to use pickled model later...
# load model from disk
#loaded_model = pickle.load( open(fname, 'rb') )
#y_pred = loaded_model.predict( X_test )
#result = loaded_model.score( X_test, y_test )
#print(result)

