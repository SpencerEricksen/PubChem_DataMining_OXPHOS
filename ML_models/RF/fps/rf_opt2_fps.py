
import sys
import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate

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

# generate a list of models to evaluate

# opt1 best params
# model, 		bal_acc (std),	avg_prc	()	roc_auc,()	F1,() 
# 4000_1.0_none_0.0009,	0.913	0.006,	0.801, 	0.017	0.970	0.004	0.171,0.004
# 4000_1.0_none_0.0014,	0.921	0.006,	0.811, 	0.016	0.972	0.003	0.211,0.003
# 4000_1.0_none_0.0018,	0.921	0.002,	0.813, 	0.015	0.973	0.003	0.199,0.001
# 4000_1.0_none_0.0023,	0.915	0.004,	0.814, 	0.016	0.972	0.005	0.185,0.001
# 4000_1.0_none_0.0027,	0.915	0.004,	0.813, 	0.015	0.972	0.004	0.178,0.002

def get_models():
    # n_estimators = (int, deafult=100)
    # max_feat = 'log2', 'sqrt', 'log2', None, float (default=1.0)
    # max_depth = None (int, default=None)
    models = dict()

    # n_trees
    for i in [ 4000, 6000, 8000 ]: 

        # max_samples (bootstrap true, fraction of samples to use for each tree)
        for j in [ 0.0002, 0.002, 0.006, 0.012, 0.05, 0.1, 0.2, 0.5, None ]:
            if j == None:
                key_j = 1.0
            else:
                key_j = j
            key_j = '%.1f' % key_j

            # max depth
            for d in [ None ]:
                if d == None:
                    key_d = "none"
                else:
                    key_d = str(d)

                # explore max_feats -- fraction from 10% to 100% in 10% increments
                # 0.020 = sqrt(N_feats) and 0.005 = log2(N_feats)
                for k in [ 0.0014, 0.0018, 0.0023 ]:
                    key_k = '%.4f' % k
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
                                             n_jobs=-1,
                                             class_weight='balanced',
                                             verbose=1 )
    return models


# get the models to evaluate
models = get_models()

f = open('rf_opt2_fps_loggers.txt','w')
print( 'ntrees_nsamples_ndepth_nfeats,bal_acc,(bal_acc),avg_prec,(avg_prec),roc_auc,(roc_auc),f1,(f1)', file=f)

for name, model in models.items():
    # evaluate the model
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
    scores = cross_validate( model, X, y, scoring=('balanced_accuracy','average_precision','roc_auc','f1'), cv=5 )
    metrics = list(scores.keys())[-4:]
    metrics_avgs = [ '{:.3f}'.format( np.mean(scores[i]) ) for i in metrics ]
    metrics_stds = [ '{:.3f}'.format( np.std(scores[i]) ) for i in metrics ]
    metrics_list = []
    for (avg, std) in zip( metrics_avgs, metrics_stds ):
        metrics_list.append( avg )
        metrics_list.append( std )
    # summarize the performance along the way
    print( "{},{}".format( name, ",".join( metrics_list ) ), file=f )

f.close()

