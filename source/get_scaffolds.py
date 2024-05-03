
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from rdkit.Chem.Scaffolds.MurckoScaffold import MakeScaffoldGeneric

# determine Bemis-Murcko scaffolds and reduced Murcko scaffolds (carbon frames)
#
# usage: python ./source/get_scaffolds.py
#    in: 'unique_pubchem_cids.pkl'
#   out: 'unique_pubchem_cids_scaffolds.pkl'


def get_murcko(smi):
    try:
        #return Chem.Scaffolds.MurckoScaffold.MurckoScaffoldSmilesFromSmiles(x)
        return MurckoScaffoldSmiles(smi)
    except:
        return None


def get_gen_murcko(smi):
    try:
        return Chem.MolToSmiles( MakeScaffoldGeneric( Chem.MolFromSmiles( smi ) ) )
        #return Chem.MolToSmiles( Chem.Scaffolds.MurckoScaffold.MakeScaffoldGeneric( Chem.MolFromSmiles(murcko) ) )
        #return Chem.MolToSmiles( Chem.Scaffolds.MurckoScaffold.MakeScaffoldGeneric(x) )
    except:
        return None

df = pd.read_pickle('unique_pubchem_cids.pkl')

# get murcko smiles from rdkit_smiles
print('getting murcko scaffolds')
df['murcko'] = df['rdkit_smiles_cln'].apply( lambda x: get_murcko(x) if x is not None else None )

print('getting generic murcko scaffolds')
df['gen_murcko'] = df['murcko'].apply( lambda x: get_gen_murcko(x) if x is not None else None )


df.to_pickle('unique_pubchem_cids_scaffolds.pkl')
