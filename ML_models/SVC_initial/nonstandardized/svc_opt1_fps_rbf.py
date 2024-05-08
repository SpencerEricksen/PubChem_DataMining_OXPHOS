
import sys
import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.model_selection import cross_validate

try:
    # feature types -- 'fps', 'desc' or 'both'
    feat_type = sys.argv[1]
except:
    print("provide feature type (desc or fps) and XXX as input arguments")
    sys.exit(1)

# load data, get features
df = pd.read_pickle('../data/training_data.pkl')
# in descriptors, had to replace infinite values with NaNs
df = df.replace( [np.inf, -np.inf], np.nan )
df = df.sample( frac=0.33 )
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
def get_models():
    # SVC params
    # C = regularization param--strength of reg inversely proportional to C (default C=1.0)
    # kernel = 'rbf' (default), 'poly', 'sigmoid', 'linear' --  something like Jaccard or Tanimoto might make sense
    # degree = 3 (default) -- degree of polynomial kernel--only used if kernel = 'poly'
    # gamma = kernel coef for 'rbf', 'poly', or 'sigmoid', 
    #       = 'scale' (default) -> 1/(n_feat * X.var()),
    #       = 'auto' -> 1/n_feat, or
    #       =  float (must be non-negative)
    models = dict()
    # reg param
    for c in [0.001, 0.01, 0.1, 1.0, 10., 100., 1000.]:
        key_c = '%.3f' % c
        #for k in ['rbf', 'poly', 'sigmoid', 'linear']:
        for k in ['rbf']:
            key_k = k
            for g in [ 'scale', 'auto', 0.1, 0.5 ]:
                key_g = str(g)
                key = str(key_c) + "_" + str(key_k) + "_" + str(key_g) + "_" + feat_type
                models[key] = SVC( C=c, kernel=k, gamma=g, class_weight='balanced', probability=True )
    return models

# get the models to evaluate
models = get_models()

f = open('svc_opt1_fps_loggers.txt','w')
print( 'c_kernel_gamma_feats,bal_acc,(bal_acc),avg_prec,(avg_prec),roc_auc,(roc_auc),f1,(f1)', file=f)

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


