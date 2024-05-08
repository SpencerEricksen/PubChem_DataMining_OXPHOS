
import pandas as pd
import numpy as np
import pickle
import math
import sys
from sklearn.preprocessing import StandardScaler
# explore random forest bootstrap sample size on performance
from sklearn.model_selection import cross_val_score
#from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import RepeatedKFold
from sklearn.ensemble import RandomForestRegressor
#from matplotlib import pyplot

# load potency data set for analogs
df = pd.read_pickle('./data/potency_chem_structs_20230201_feats.pkl')

df = df.drop( columns='desc_Ipc' )
# remove outlier
df = df.loc[ ~(df['ID'] == 'RZLT-215') ]

# keep only data set with measured IC50s
df = df.dropna( subset=['pIC50'] )

# get rdkit mol descriptors
desc_no_frag_cols = [ i for i in df.columns if ( 'desc_' in i ) and ( '_fr_' not in i ) ]
Z = np.stack( df[desc_no_frag_cols].values )
Z_scaled = StandardScaler().fit_transform(Z)
df_Z = pd.DataFrame(Z_scaled, columns=desc_no_frag_cols, index=df['ID'])

# get ECFP6 features
df['morgan'] = df['ECFP6_bitstr'].apply( lambda x: [int(i) for i in list(x)] )
Y = np.stack( df['morgan'].values )
# compress fingerprint to meaningful bit indices
df_Y = pd.DataFrame(Y, index=df['ID'])
for c in df_Y.columns:
    # remove constant bits
    if (df_Y[c].sum() == 63) or (df_Y[c].sum() == 0):
        df_Y.drop( columns=c, inplace=True )
#bit_indices = df_Y.columns.to_list()

# merge, get feature array
df_X = pd.concat( [df_Z, df_Y], axis=1 )
X = df_X.values

# labels
y = df['pIC50'].values

# generate a list of models to evaluate
def get_models():
    # n_estimators = (int, deafult=100)
    # max_feat = 'log2', 'sqrt', 'log2', None, float (default=1.0)
    # max_depth = None (int, default=None)
    # criterion = 'absolute error', 'squared error', 'poisson' (default='squared error')
    models = dict()

    # n_trees
    for i in [ 800, 1200, 1600, 2000 ]:

        # max_samples (bootstrap true, fraction of samples to use for each tree)
        for j in [ None ]:
        #for j in [0.2, 0.4, 0.6, 0.8, None]:
            if j == None:
                key_j = 1.0
            else:
                key_j = j
            key_j = '%.1f' % key_j

            # max depth
            #for d in [1, 2, 3, 6, None]:
            for d in [ None ]:
                if d == None:
                    key_d = "none"
                else:
                    key_d = str(d) 

                # explore max_feats -- fraction from 10% to 100% in 10% increments
                for k in [ 0.20, 0.25, 0.30, 0.33 ]:
                    key_k = '%.2f' % k 
                    key = str(i)+"_"+str(key_j)+"_"+key_d+"_"+str(key_k)

                    # left out "criterion" setting
                    models[key] = RandomForestRegressor( n_estimators=i, 
                                             max_depth=d, 
                                             max_samples=j,
                                             min_samples_split=2, 
                                             min_samples_leaf=1, 
                                             max_features=k, 
                                             max_leaf_nodes=None, 
                                             min_impurity_decrease=0.0, 
                                             bootstrap=True, 
                                             oob_score=False, 
                                             n_jobs=18, 
                                             verbose=1 )
    return models
 
# evaluate a given model using cross-validation
def evaluate_model(model, X, y):
    # define the evaluation procedure
    cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
    # evaluate the model and collect the results
    scores = cross_val_score(model, X, y, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1, error_score='raise')
    return scores
 
# get the models to evaluate
models = get_models()
# evaluate the models and store results
results, names = list(), list()
print('>ntrees_nsamples_ndepth_nfeats,MAE_avg,MAE_std')
for name, model in models.items():
    # evaluate the model
    scores = evaluate_model(model, X, y)
    # store the results
    results.append(scores)
    names.append(name)
    # summarize the performance along the way
    print('>%s,%.3f,(%.3f)' % (name, np.mean(scores), np.std(scores)))

# plot model performance for comparison
#pyplot.boxplot(results, labels=names, showmeans=True)
#pyplot.savefig('rfr_explore_maxdepth.png', dpi=600 )
#pyplot.show()
