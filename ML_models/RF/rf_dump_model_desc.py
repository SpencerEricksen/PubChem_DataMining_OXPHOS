
import sys
import pandas as pd
import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier

try:
    # feature types -- 'fps', 'desc' or 'both'
    feat_type = sys.argv[1]
except:
    print("provide feature type: 'desc', 'fps', or 'both'")
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

# best hyperparams for desc 
# ntrees	8000
# nsamples	1.0
# ndepth	none
# nfeats	0.050

# generate model with optimal hyperparams 
model = RandomForestClassifier( 
          n_estimators=8000,
          criterion='gini',
          max_depth=None,
          max_samples=None,
          min_samples_split=2,
          min_samples_leaf=1,
          max_features=0.050,
          max_leaf_nodes=None,
          min_impurity_decrease=0.0,
          bootstrap=True,
          oob_score=False,
          n_jobs=-1,
          class_weight='balanced',
          verbose=1 )

# train on full train set (80% of the total data)
model.fit( X, y )

# save the model to disk
model_filename = "model_RF_desc_8000_1.0_none_0.050_3way.pkl"
pickle.dump( model, open( model_filename, 'wb'))

