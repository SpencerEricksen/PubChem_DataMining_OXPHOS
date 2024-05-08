
import sys
import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.model_selection import cross_validate
#from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import f1_score, average_precision_score, roc_auc_score, balanced_accuracy_score

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
#df = df.sample( frac=0.50 )
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

'''
#dists = pdist( X, metric='jaccard' )
dists = pairwise_distances( X, metric='jaccard', n_jobs=-1 )
sims = 1 - pdist
# get similarity matrix as SVC kernel
#similarity_mat = squareform( sims )
clf = SVC( kernel=sims)
clf.fit(X,y)
from scipy.spatial.distance import jaccard
clf = SVC( kernel=jaccard, class_weight='balanced' )
clf.fit(X,y)
'''

# SVC params
# C = regularization param--strength of regularization is inversely proportional to C
# kernel = linear?  something like Jaccard or Tanimoto would make more sense
# degree = 3 (default) -- degree of polynomial kernel--only used if kernel = 'poly'
# gamma = kernel coef for kernels 'rbf', 'poly', or 'sigmoid'
# class_weight = 'balanced'

#clf = SVC( kernel='linear', C=1.0, class_weight='balanced', probability=True )
#svc_linear = clf.fit(X,y)

clf = SVC( kernel='rbf', gamma=0.001, C=100.0, class_weight='balanced', probability=True )
svc_rbf = clf.fit(X,y)

# load test data, evaluate model
df = pd.read_pickle('../data/test_data.pkl')

if feat_type == 'fps':
    feat_cols = [ c for c in col_list[2:] if 'desc' not in str(c) ]
    X = df[feat_cols].values

y = df.label.values
#weird issue with int type in y, had to set type to match X
y = np.array( y, np.int64 )

y_pred = svc_rbf.predict(X)
y_score = svc_rbf.predict_proba(X)
y_score = y_score[:,1]

# binary classification metrics
f1 = f1_score( y, y_pred )
bac = balanced_accuracy_score( y, y_pred )

# continuous classification metrics
avgprec = average_precision_score( y, y_score )
rocauc = roc_auc_score( y, y_score )

# print results
print('model,balanced_acc,avg_prec,roc_auc,F1')
print( 'svc_rbf_gamma0.001_C100.0', bac, avgprec, rocauc, f1)

