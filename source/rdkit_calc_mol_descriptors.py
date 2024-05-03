import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem

def get_desc_list( m ):
    '''return list of RDKit mol descriptors from rdkit mol object'''
    desc_list = []
    for tup in Descriptors.descList:
        try:
            desc_list.append(tup[1](m))
        except:
            desc_list.append(None)
            pass
    return desc_list

df = pd.read_pickle('unique_pubchem_cids.pkl')

# get rdkit mol objects based on rdkit smiles in dataframe
df['rdkit_mol'] = df['rdkit_smiles_cln'].apply( lambda x: Chem.MolFromSmiles(x) if x is not None else None )

# generate feature headers
desc_headers = []
for tup in Descriptors.descList:
    desc_headers.append( tup[0] )
desc_headers = [ "desc_" + d for d in desc_headers ]

# generate features for each rdkit mol object
df['rdkit_desc'] = df['rdkit_mol'].apply( lambda x: get_desc_list(x) if x is not None else None )

# store array of features as a list of lists (300,000+ cpds by 208 features)
X_desc = list( df['rdkit_desc'] )

# build a list of None values to use when a compound can't be featurized (208 features)
k = []
for i in range(208):
    k.append( None )

# if raw is None, replace with a list of 208 None values (missing features)
for i,j in enumerate(X_desc):
    if j == None:
        X_desc[i] = k

# create dataframe for features for cpd set
df_desc = pd.DataFrame( X_desc, columns=desc_headers )
df.drop( columns=['rdkit_mol', 'rdkit_desc'], inplace=True )

# add PUBCHEM_CIDs to descriptor dataframe (df_desc)
molid = df['PUBCHEM_CID'].tolist()
df_desc.index = molid
df_desc = df_desc.reset_index()
df_desc.rename( columns={'index':'PUBCHEM_CID'}, inplace=True )

# merge descriptor dataframe (df_desc) on to cpd smiles dataframe (df)
df = df.merge( df_desc, on='PUBCHEM_CID', how='left')

# save as a pickle
df.to_pickle('unique_pubchem_cids_desc.pkl')

