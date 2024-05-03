
import pandas as pd
import numpy as np

df = pd.read_pickle('oxphos_merged_aids_cids_clnsmi_assaydesc_ETCassay_pmid.pkl')

df2 = df[['PUBCHEM_CID','rdkit_smiles_cln']]

df3 = df2.drop_duplicates( subset='PUBCHEM_CID', keep='first' )
df3.to_pickle('unique_pubchem_cids.pkl')

