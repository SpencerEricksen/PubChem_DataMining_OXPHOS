
import sys
import pandas as pd
import numpy as np
import pickle
from sklearn.svm import SVC


try:
    feat_type = sys.argv[1]    # feature types -- 'fps', 'desc' or 'both'
    kernel_type = sys.argv[2]  # 'linear','rbf','sigmoid','poly'
    c = float(sys.argv[3])     # 'regularization factor C'
    g = sys.argv[4]     # kernel coef for 'rbf', 'poly', 'sigmoid'
    cf = float(sys.argv[5])    # term added to 'sigmoid' or 'poly' kernels
    d = int(sys.argv[6])       # degree of kernel function in 'poly'
except:
    print("provide feature type: 'desc', 'fps', or 'both'")
    sys.exit(1)

# load data, get features
df = pd.read_pickle('../data/training_data.pkl')
# in descriptors, had to replace infinite values with NaNs
df = df.replace( [np.inf, -np.inf], np.nan )
# to get around issues of memory overflow
#df = df.sample( frac=0.33 )
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
    scaler = pickle.load( open('../data/stand_scaler_fps_training_data.pkl', 'rb') )
    X = scaler.transform(X)
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
    scaler = pickle.load( open('../data/stand_scaler_fps_training_data.pkl', 'rb') )
    X2 = scaler.transform(X2)
    df_X2 = pd.DataFrame( X2, columns=fps_cols, index=molid_list )
    # merge desc and fps into single dataframe, get feature array
    df_X12 = pd.concat( [df_X1, df_X2], axis=1 )
    X = df_X12.values

# labels
y = df.label.values
#weird issue with int type in y, had to set type to match X
y = np.array( y, np.int64 )

# SVC params
# C = regularization param--strength of reg inversely proportional to C (default C=1.0)
# kernel = 'rbf' (default), 'poly', 'sigmoid', 'linear' --  something like Jaccard or Tanimoto might make sense
# degree = 3 (default) -- degree of polynomial kernel--only used if kernel = 'poly'
# gamma = kernel coef for 'rbf', 'poly', or 'sigmoid',
#       = 'scale' (default) -> 1/(n_feat * X.var()),
#       = 'auto' -> 1/n_feat, or
#       =  float (must be non-negative)

key_c = '{:.3f}'.format(c)

if type(g) is str:
    key_g = g
    if g == 'none':
        g = None
else:
    key_g = '{:.3f}'.format(g)

if cf == 'none':
    key_cf = cf
    cf = None
else:
    key_cf = '{:.3f}'.format(cf)

if d == 'none':
    key_d = d
    d = None
else:
    key_d = str(d)


model_name = "svc_" + kernel_type + "_" + feat_type + "_C_" + key_c + "_gamma_" + key_g + "_coef0_" + key_cf + "_degree_" + key_d
model = SVC( C=c, kernel=kernel_type, gamma=g, coef0=cf, degree=d, class_weight='balanced', probability=True )

# train on full train set (80% of the total data)
model.fit( X, y )

# save the model to disk
model_filename = "model_" + model_name + ".pkl"
pickle.dump( model, open( model_filename, 'wb'))

