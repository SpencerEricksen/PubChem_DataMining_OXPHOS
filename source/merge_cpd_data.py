import pandas as pd

df1 = pd.read_pickle('unique_pubchem_cids_complete_assay_results_w_labels.pkl')
df2 = pd.read_pickle('unique_pubchem_cids_desc.pkl')
df3 = pd.read_pickle('unique_pubchem_cids_fps.pkl')
df4 = pd.read_pickle('unique_pubchem_cids_lipinski.pkl')
df5 = pd.read_pickle('unique_pubchem_cids_npscores.pkl')
df6 = pd.read_pickle('unique_pubchem_cids_pains.pkl')
df7 = pd.read_pickle('unique_pubchem_cids_scaffolds.pkl')

desc_cols = [ i for i in df2.columns if "desc_" in i ]

dfx = df1.merge( df2[['PUBCHEM_CID'] + desc_cols ], on='PUBCHEM_CID', how='left' )
dfx = dfx.merge( df3[['PUBCHEM_CID','ECFP6_bitstr']], on='PUBCHEM_CID', how='left' )
dfx = dfx.merge( df4[['PUBCHEM_CID','lipinski']], on='PUBCHEM_CID', how='left' )
dfx = dfx.merge( df5[['PUBCHEM_CID','npscore']], on='PUBCHEM_CID', how='left' )
dfx = dfx.merge( df6[['PUBCHEM_CID','pains']], on='PUBCHEM_CID', how='left' )
dfx = dfx.merge( df7[['PUBCHEM_CID','murcko','gen_murcko']], on='PUBCHEM_CID', how='left' )

dfx.to_pickle('unique_pubchem_cids_all.pkl')
