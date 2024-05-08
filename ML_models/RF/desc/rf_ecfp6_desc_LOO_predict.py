
import pandas as pd
import numpy as np
import pickle
import math
import sys
# explore random forest bootstrap sample size on performance
#from sklearn.model_selection import cross_val_score
#from sklearn.model_selection import RepeatedStratifiedKFold
#from sklearn.model_selection import RepeatedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.ensemble import RandomForestRegressor
#from matplotlib import pyplot

# load potency data set for analogs
df = pd.read_pickle('./data/potency_chem_structs_20230201_feats.pkl')

# remove outlier
df = df.loc[ ~(df['ID'] == 'RZLT-215') ]
df = df.loc[ ~(df['ID'] == 'RZ-32') ]
df = df.loc[ ~(df['ID'] == 'RZLT-306') ]
df = df.loc[ ~(df['ID'] == 'RZLT-307') ]

df = df.drop( columns='desc_Ipc' )

# keep only data set with measured IC50s
df = df.dropna( subset=['pIC50'] )

# get rdkit mol descriptors
desc_no_frag_cols = [ i for i in df.columns if ( 'desc_' in i ) and ( '_fr_' not in i ) ]
Z = np.stack( df[desc_no_frag_cols].values )
Z_scaled = StandardScaler().fit_transform(Z)
df_Z = pd.DataFrame(Z_scaled, columns=desc_no_frag_cols, index=df['ID'])

# get ECFP6 features
df['morgan'] = df['ECFP6_bitstr'].apply( lambda x: [int(i) for i in list(x)] )
Y = np.stack( df['morgan'].values )
# compress fingerprint to meaningful bit indices
df_Y = pd.DataFrame(Y, index=df['ID'])
for c in df_Y.columns:
    # remove constant bits
    if (df_Y[c].sum() == 63) or (df_Y[c].sum() == 0):
        df_Y.drop( columns=c, inplace=True )
#bit_indices = df_Y.columns.to_list()

# merge rdkit mol desc and ecfp, get feature array
df_X = pd.concat( [df_Z, df_Y], axis=1 )
X = df_X.values

# labels
y = df['pIC50'].values

# molids list
molids = df['ID'].tolist()

# from hyperparameter search--best config
# >ntrees_nsamples_ndepth_nfeats,MAE_avg,MAE_std
# 2000_1.0_none_0.20,-0.483,(0.116)
model = RandomForestRegressor( n_estimators=2000,
                               max_depth=None,
                               max_samples=None,
                               min_samples_split=2,
                               min_samples_leaf=1,
                               max_features=0.20,
                               max_leaf_nodes=None,
                               min_impurity_decrease=0.0,
                               bootstrap=True,
                               oob_score=False,
                               n_jobs=18,
                               verbose=1,
                               random_state=None )

# create loocv procedure
cv = LeaveOneOut()

# enumerate splits
y_true, y_pred = list(), list()

for train_ix, test_ix in cv.split(X):
    #split data
    X_train, X_test = X[train_ix, :], X[test_ix, :]
    y_train, y_test = y[train_ix], y[test_ix]
    # fit model
    model.fit(X_train, y_train)
    # evaluate model
    yhat = model.predict(X_test)
    y_true.append( y_test[0] )
    y_pred.append( yhat[0] )
    print( y_test[0], yhat[0] )


# dump LOO prediction data
df2 = pd.DataFrame( data=zip( molids, y_true, y_pred ), columns=['molid','pIC50','pred_pIC50'] )
df2.to_csv('loo_evaluations.csv', index=False )

print( '*' * 50 )
print('pIC50,pred_pIC50')
for i in zip( y_true, y_pred ):
    print( '{:.3f},{:.3f}'.format(i[0],i[1]) )


