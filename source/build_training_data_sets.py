
import pandas as pd
import numpy as np

# read cpd data
df1 = pd.read_pickle('unique_pubchem_cids_all_wclusts.pkl')

# get cpds will clear +/- labels, pass Lipinksi-like filter, and have been in ETC-linked AIDs
df2 = df1.loc[ ((df1['label'] == 0) | (df1['label'] == 1)) & (df1['lipinski'] == 1) & (df1['num_etc_linked_aids'] > 0) & (df1['pains'] == False) ]

# on second thought, skip the lipinski filter on the training/test data
#df2 = df1.loc[ ((df1['label'] == 0) | (df1['label'] == 1)) & (df1['num_etc_linked_aids'] > 0) ]

# get fps dataset
df_fps = df2[['PUBCHEM_CID','label','ECFP6_bitstr']]
df_fps['morgan'] = df_fps['ECFP6_bitstr'].apply( lambda x: [int(i) for i in list(x) ] )
X = np.stack( df_fps['morgan'].values )
df_X = pd.DataFrame( X, index=df_fps['PUBCHEM_CID'] )
labels = df_fps['label'].values
df_X['label'] = labels
df_X = df_X.reset_index()
df_X = df_X[['PUBCHEM_CID','label'] + list(range(2048)) ]
df_X.to_pickle('data_set_fps.pkl')

# get desc dataset
desc_cols = [ i for i in df2.columns if 'desc_' in i ]
df_desc = df2[ ['PUBCHEM_CID','label'] + desc_cols ]
df_desc.to_pickle('data_set_desc.pkl')

