import pandas as pd
import numpy as np
from rdkit import Chem
#from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import AllChem


def get_morgan_r3(m):
    '''return fingerprint object from rdkit mol'''
    try:
        return AllChem.GetMorganFingerprintAsBitVect( m, 3, nBits=2048 )
    except:
        return None


def make_bitstring(fp):
    try:
        return "".join( [ str(i) for i in fp ] )
    except:
        return None

df = pd.read_pickle('unique_pubchem_cids.pkl')

# get rdkit mol objects based on rdkit smiles in dataframe
df['rdkit_mol'] = df['rdkit_smiles_cln'].apply( lambda x: Chem.MolFromSmiles(x) if x is not None else None )

# generate bitstring fingerprints from mol objects
df['ECFP6_bitvec'] = df['rdkit_mol'].apply( lambda m: get_morgan_r3(m) if m is not None else None )
df['ECFP6_bitstr'] = df['ECFP6_bitvec'].apply( lambda fp: make_bitstring(fp) if fp is not None else None )

df.drop( columns=['rdkit_mol','ECFP6_bitvec'], inplace=True )

# save as a pickle
df.to_pickle('unique_pubchem_cids_fps.pkl')

