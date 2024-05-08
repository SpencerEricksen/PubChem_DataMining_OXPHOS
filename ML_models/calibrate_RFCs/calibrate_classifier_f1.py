import sys
import pickle
import pandas as pd
import numpy as np
from sklearn.metrics import f1_score, matthews_corrcoef

try:
    feat_type = sys.argv[1]    # 'desc', 'fps', or 'both'
    model_name = sys.argv[2]   # "model_l2_lbfgs_C0.00599.pkl"
except:
    print('usage:    python calibrate.py   feat_type   model_name')
    sys.exit(1)

# apply threshold to positive probabilities to create labels
def to_labels(pos_probs, threshold):
    return (pos_probs >= threshold).astype('int')

# load calibration data for model calibration based on F1 or MCC
df = pd.read_pickle('../data/calibrate_data.pkl')
molid_list = df.PUBCHEM_CID.tolist()

# features
col_list = df.columns.tolist()
if feat_type == 'desc':
    feat_cols = [ c for c in col_list[2:] if 'desc' in str(c) ]
    feat_cols.remove('desc_Ipc')
    df_X = df[feat_cols]
    # replace infinities with np.nan
    df_X.replace( [np.inf, -np.inf], np.nan, inplace=True )
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
    # replace infinities with np.nan
    df_X1.replace( [np.inf, -np.inf], np.nan, inplace=True )
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

# load model
mymodel = pickle.load( open(model_name, 'rb') )

# predict probabilities
yhat = mymodel.predict_proba(X)

# keep probabilitites for the positive outcomes
probs = yhat[:,1]

# define thresholds
thresholds = np.arange(0, 1, 0.001)

# evaluate each threshold
scores = [ f1_score( y, to_labels(probs, t) ) for t in thresholds ]
#scores = [ matthews_corrcoef(  y, to_labels(probs, t) ) for t in thresholds ]

# get best threshold
ix = np.argmax(scores)

print('model:{}, feat_type:{}, threshold:{:.3f}, F1:{:.3f}'.format( model_name, feat_type, thresholds[ix], scores[ix] ) )
#print('model:{}, feat_type:{}, threshold:{:.3f}, MCC:{:.3f}'.format( model_name, feat_type, thresholds[ix], scores[ix] ) )

#print('Threshold=%.3f, F-Score=%.5f' % (thresholds[ix], scores[ix]))

y_pred = to_labels( probs, thresholds[ix] ) 
df2 = pd.DataFrame( {'PUBCHEM_CID':molid_list, 'y':y, 'yhat_'+str(thresholds[ix]):y_pred, 'probs':probs } )
df2.to_csv( "_".join( model_name.split('/')[-1].split('.')[0:-1] ) + "cal_labels_probs_optF1.csv", index=False ) 
