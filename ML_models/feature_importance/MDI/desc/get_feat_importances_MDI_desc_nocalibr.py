import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd

# see https://scikit-learn.org/stable/auto_examples/ensemble/plot_forest_importances.html


# load model for feature importance calcs
model_name = "../../../RF/model_RF_desc_8000_1.0_none_0.050.pkl"
mymodel = pickle.load( open(model_name, 'rb') )

feat_type = "desc"

# load test data (just need this for feat col names for now)
df = pd.read_pickle('../../../data/old_85_15_split/test_data.pkl')

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
    scaler = pickle.load( open('../../../data/old_85_15_split/stand_scaler_desc_training_data.pkl', 'rb') )
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
    scaler = pickle.load( open('../../../data/old_85_15_split/stand_scaler_desc_training_data.pkl', 'rb') )
    X1 = scaler.transform( df_X1 )
    df_X1 = pd.DataFrame( X1, columns=desc_cols, index=molid_list )
    # get data frame of fps features
    fps_cols = [ c for c in col_list[2:] if 'desc' not in str(c) ]
    X2 = df[fps_cols].values
    df_X2 = pd.DataFrame( X2, columns=fps_cols, index=molid_list )
    # merge desc and fps into single dataframe, get feature array
    df_X12 = pd.concat( [df_X1, df_X2], axis=1 )
    X = df_X12.values

# get impurity-based importance values for each feature (and stdev based on individual trees)
importances = mymodel.feature_importances_
std = np.std( [ tree.feature_importances_ for tree in mymodel.estimators_ ], axis=0 )

# plot feature importances using MDI
df3 = pd.DataFrame( [importances, std, feat_cols] ).T
df3.columns = ['mean_decrease_impurity','stdev','feature']
df4 = df3.sort_values( by='mean_decrease_impurity', ascending=False )

df4.to_csv('RF_nocalibr_desc_MDI.csv', index=False )

# plot feature importances using MDI for top20 features (any)
df_top20 = df4.head(20)
df_top20.set_index( 'feature', inplace=True )

fig, ax = plt.subplots()
df_top20['mean_decrease_impurity'].plot.bar(yerr=df_top20.stdev, ax=ax)
ax.set_title("MDI feature importances: top 20 molecular descriptors")
ax.set_ylabel("mean decrease in impurity")
fig.tight_layout()
plt.savefig('RFC_nocalibr_all_desc_feat_MDI_top20.png', dpi=1600 )
plt.close()

#forest_importances = pd.Series( importances, index=feat_cols )
#forest_importances.plot.bar(yerr=std, ax=ax)

# plot feature importances using MDI for fragment-based features
df_fr = df4.loc[ df4['feature'].str.contains('_fr_') ]

fig, ax = plt.subplots()
df_fr.set_index('feature', inplace=True )
df_fr['mean_decrease_impurity'].plot.bar( ax=ax )
ax.set_title("MDI feature importances: ")
ax.set_ylabel("mean decrease in impurity")
ax.set_xticklabels( df_fr.index, fontsize=1)
fig.tight_layout()
plt.savefig('RFC_nocalibr_desc_feat_MDI_frags.png', dpi=1600 )
plt.close()

# plot property feature importances
df_non_fr = df4.loc[ ~df4['feature'].str.contains('_fr_') ]

fig, ax = plt.subplots()
df_non_fr.set_index('feature', inplace=True )
df_non_fr['mean_decrease_impurity'].plot.bar( ax=ax )
ax.set_title("Feature importances using MDI")
ax.set_ylabel("Mean decrease in impurity")
ax.set_xticklabels( df_non_fr.index, fontsize=1)
fig.tight_layout()
plt.savefig('RFC_nocalibr_desc_feat_MDI_props.png', dpi=1600 )
plt.close()
