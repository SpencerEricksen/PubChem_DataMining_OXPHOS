
import sys
import pickle
import pandas as pd
import numpy as np

# load model
feat_type = 'fps'
model_file = 'model_RF_fps_8000_1.0_none_0.0023.pkl'
mymodel = pickle.load( open(model_file, 'rb') )

# load held-out test data for model evaluation
df = pd.read_pickle('../data/old_85_15_split/test_data.pkl')
molid_list = df.PUBCHEM_CID.tolist()

# features
col_list = df.columns.tolist()
feat_cols = [ c for c in col_list[2:] if 'desc' not in str(c) ]
X = df[feat_cols].values
df_fps = df[ ['PUBCHEM_CID','label'] + feat_cols ]
df_fps['PUBCHEM_CID'] = df_fps['PUBCHEM_CID'].astype('int')

# labels
#y = df.label.values
#weird issue with int type in y, had to set type to match X
#y = np.array( y, np.int64 )

# get predicted labels (0/1)
#y_pred = mymodel.predict(X)
# get predicted probabilities (needed for rocauc and prauc)
y_score = mymodel.predict_proba(X)[:,1]

df_scores = pd.DataFrame( y_score, molid_list, columns=['score'] )
df_scores.reset_index( inplace=True )
df_scores.rename( columns={'index':'PUBCHEM_CID'}, inplace=True )
df_scores['PUBCHEM_CID'] = df_scores['PUBCHEM_CID'].astype('int')

df_all = df_scores.merge( df_fps, on='PUBCHEM_CID', how='left' )
df_all.sort_values( by=['label','score'], ascending=[False,False], inplace=True )


# read in smiles for all oxphos cpds
df_smi = pd.read_pickle('../../unique_pubchem_cids.pkl')

df_all.merge( df_smi, on='PUBCHEM_CID', how='left', inplace=True)

df_all.to_pickle( 'scores_RF_fps_test_set.pkl' )


