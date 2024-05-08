
import sys
import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
#from sklearn.metrics import balanced_accuracy_score
#from sklearn.metrics import roc_auc_score
#from sklearn.metrics import f1_score
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_score

try:
    # feature types -- 'fps', 'desc' or 'both'
    feat_type = sys.argv[1]
    # type of regularization (norm of the penalty) 'l1' or 'l2'
    #reg = sys.argv[2]
except:
    print("provide feature type (desc or fps) and XXX as input arguments")
    sys.exit(1)

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

# labels
y = df.label.values
#weird issue with int type in y, had to set type to match X
y = np.array( y, np.int64 )

# generate a list of models to evaluate
def get_models():
    # n_estimators = (int, deafult=100)
    # max_feat = 'log2', 'sqrt', 'log2', None, float (default=1.0)
    # max_depth = None (int, default=None)
    # criterion = 'absolute error', 'squared error', 'poisson' (default='squared error')
    models = dict()

    # n_trees
    for i in [ 200, 600 ]:

        # max_samples (bootstrap true, fraction of samples to use for each tree)
        #for j in [ None ]
        for j in [0.2, 0.5, None]:
            if j == None:
                key_j = 1.0
            else:
                key_j = j
            key_j = '%.1f' % key_j

            # max depth
            for d in [ 2, 6, None]:
                if d == None:
                    key_d = "none"
                else:
                    key_d = str(d)

                # explore max_feats -- fraction from 10% to 100% in 10% increments
                for k in [ 0.25, 0.45 ]:
                    key_k = '%.2f' % k
                    key = str(i)+"_"+str(key_j)+"_"+key_d+"_"+str(key_k)

                    # left out "criterion" setting
                    models[key] = RandomForestClassifier( 
                                             n_estimators=i,
                                             criterion='gini',
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
                                             class_weight='balanced',
                                             verbose=1 )
    return models

# evaluate a given model using cross-validation
def evaluate_model(model, X, y):
    # define the evaluation procedure
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
    # evaluate the model and collect the results
    scores = cross_val_score(model, X, y, scoring='balanced_accuracy', cv=cv, n_jobs=-1, error_score='raise')
    return scores

# print performanced in 5-CV

# get the models to evaluate
models = get_models()
# evaluate the models and store results
results, names = list(), list()
print('>ntrees_nsamples_ndepth_nfeats,bal_acc_mean,bal_acc_std')
for name, model in models.items():
    # evaluate the model
    scores = evaluate_model(model, X, y)
    # store the results
    results.append(scores)
    names.append(name)
    # summarize the performance along the way
    print('>%s,%.3f,(%.3f)' % (name, np.mean(scores), np.std(scores)))
