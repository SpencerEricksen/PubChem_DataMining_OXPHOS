
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize

# need an alternative FRAGMENT REMOVER since some
# non-standard counterions are not RDKit's default
# salt fragment list: $RDBASE/Data/Salts.txt

def alt_frag_remover( m ):
    # split m into mol fragments, keep fragment with highest num atoms
    mols = list(Chem.GetMolFrags( m, asMols=True ))
    if (mols):
        mols.sort(reverse=True, key=lambda x: x.GetNumAtoms() )
        mol = mols[0]
    else:
        mol = None
    return mol


# read in smiles file as dataframe
df = pd.read_csv( 'Smiles_FDA.csv', low_memory=False )

# remove lines with missing smiles strings (2 out of 100482)
df.dropna( subset=['smiles'], inplace=True )

# make column to store mol objects
PandasTools.AddMoleculeColumnToFrame( df, 'smiles', 'rdkit_mol', includeFingerprints=False )

# remove lines where smiles string could not generate mol object
df.dropna( subset=['rdkit_mol'], inplace=True )

_saltRemover = SaltRemover.SaltRemover()
df['rdkit_mol'] = df['rdkit_mol'].apply( lambda x: _saltRemover.StripMol(x) if x is not None else None)

# apply secondary salt remover
df['rdkit_mol'] = df['rdkit_mol'].apply( lambda x: alt_frag_remover(x) if x is not None else None )

# standardize charges
un = rdMolStandardize.Uncharger()
df['rdkit_mol'] = df['rdkit_mol'].apply( lambda x: un.uncharge(x) if x is not None else None)

# replace smiles with clean smiles
df['rdkit_smiles'] = df['rdkit_mol'].apply( lambda x: Chem.MolToSmiles(x) if x is not None else None )

# dump to CSV
# drop the mol objects columns if dumping to CSV
df.drop( columns=['rdkit_mol'], inplace=True )
df.to_csv( 'Smiles_FDA_cln.csv', index=False )

