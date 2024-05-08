import numpy as np
#import matplotlib.pyplot as plt
import pickle
import pandas as pd
from sklearn.inspection import permutation_importance
import sys
import time

# see https://scikit-learn.org/stable/auto_examples/ensemble/plot_forest_importances.html
'''
try:
    feat_type = sys.argv[1] # "desc" or "fps"
    model_file = sys.argv[2] # "../RF/model_RF_desc_8000_1.0_none_0.050.pkl"
except:
    print('')
    print('input   feat_type      and   model_file                                      arguments')
    print( '       (desc or fps)  and   "../RF/model_RF_desc_8000_1.0_none_0.050.pkl"           ')
    print('')
    exit()
'''

# specify feature type to analyze for importances
feat_type = 'desc'

# load model
model_file = '../calibrate_RFCs/model_RF_desc_8000_1_0_none_0_050_3way_calibrated.pkl'
#model_file = '../RF/model_RF_desc_8000_1.0_none_0.050.pkl'
mymodel = pickle.load( open(model_file, 'rb') )

# load test data (just need this for feat col names for now)
df = pd.read_pickle('../data/test_data.pkl')

molid_list = df.PUBCHEM_CID.tolist()

y = df.label.values
y = np.array( y, dtype='int' )

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

# calculate permutation-based feature importances
start_time = time.time()
result = permutation_importance( mymodel, X, y, n_repeats=10, random_state=42, scoring='average_precision' )
elapsed_time = time.time() - start_time
print(f"Elapsed time to compute the importances: {elapsed_time:.3f} seconds")
importances = result.importances_mean
std = result.importances_std
df3 = pd.DataFrame( [importances, std, feat_cols] ).T
df3.columns = ['perm_feat_importance','stdev','feature']
df4 = df3.sort_values( by='perm_feat_importance', ascending=False )

df4.to_csv('RF_calibr_desc_perm_feat_importances_test_avgprec.csv', index=False )


'''
# get top 20 most important features
df_top20 = df4.head(20)
df_top20.set_index( 'feature', inplace=True )

# plot top 20 permutation-based importances
fig, ax = plt.subplots()
df_top20['perm_feat_importance'].plot.bar(yerr=df_top20.stdev, ax=ax)
ax.set_title( "Feature importances by Permutation "+feat_type )
ax.set_ylabel("Mean permutation importance")
fig.tight_layout()
plt.savefig('RFC_'+feat_type+'_feat_importances_permutation_top20.png', dpi=1600 )
plt.close()

# plot feature importances using MDI for fragment-based features
df_fr = df4.loc[ df4['feature'].str.contains('_fr_') ]

fig, ax = plt.subplots()
df_fr.set_index('feature', inplace=True )
df_fr['perm_feat_importance'].plot.bar( ax=ax )
ax.set_title("Feature importances by Permutation")
ax.set_ylabel("Mean permutation importance")
ax.set_xticklabels( df_fr.index, fontsize=1)
fig.tight_layout()
plt.savefig('RFC_desc_feat_importances_permutation_frags.png', dpi=1600 )
plt.close()

# plot property feature importances
df_non_fr = df4.loc[ ~df4['feature'].str.contains('_fr_') ]

fig, ax = plt.subplots()
df_non_fr.set_index('feature', inplace=True )
df_non_fr['perm_feat_importance'].plot.bar( ax=ax )
ax.set_title("Feature importances by Permutation")
ax.set_ylabel("Mean permutation importance")
ax.set_xticklabels( df_non_fr.index, fontsize=1)
fig.tight_layout()
plt.savefig('RFC_desc_feat_importances_permutation_props.png', dpi=1600 )
plt.close()

'''
