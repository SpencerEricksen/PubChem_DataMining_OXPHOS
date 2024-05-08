
import sys
import pickle
import pandas as pd
import numpy as np
from sklearn.metrics import f1_score, average_precision_score, roc_auc_score, balanced_accuracy_score
import sys

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
    return ef, ef/ef_max

try:
    feat_type = sys.argv[1]  # 'desc', 'fps', or 'both'
    model_name = sys.argv[2] # "model_l2_lbfgs_C0.00599.pkl"
except:
    print('provide feat_type and model_file (pkl) as arguments')
    sys.exit(1)

# load model
mymodel = pickle.load( open(model_name, 'rb') )

# load held-out test data for model evaluation
df = pd.read_pickle('../data/test_data.pkl')
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

# get predicted labels (0/1)
y_pred = mymodel.predict(X)

# get predicted probabilities (needed for rocauc and prauc)
y_score = mymodel.predict_proba(X)[:,1]
#y_score = loaded_model.decision_function(X)

# binary classification metrics
f1 = f1_score( y, y_pred )
bac = balanced_accuracy_score( y, y_pred )

# continuous classification metrics
avgprec = average_precision_score( y, y_score )
rocauc = roc_auc_score( y, y_score )
ef1, nef1 = normalized_enrichment_factor_single( y, y_score, 0.01 )
print('model,balanced_acc,avg_prec,roc_auc,F1,EF1,NEF1')
print(model_name, bac, avgprec, rocauc, f1, ef1, nef1)
