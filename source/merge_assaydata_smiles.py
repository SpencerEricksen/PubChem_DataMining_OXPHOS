
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

# substructure to search
#frag = Chem.MolFromSmarts('C(=O)C=CC')
#frag = Chem.MolFromSmarts('C(=O)')
#frag = Chem.MolFromSmarts('C(=O)C')
#frag = Chem.MolFromSmarts('C(=O)C=[C,N]')
#frag = Chem.MolFromSmarts('C(=O)C=C')


# get AID data
df_act = pd.read_pickle( './oxphos_merged_aids_records.pkl' )

# get smiles
df_smi = pd.read_csv('./fetch_CGIs_bulk/pcba_oxphos_all_cids.smi', header=None, sep='\t', names=['PUBCHEM_CID','smiles'] )

# add smiles to the AID data
df = df_act.merge( df_smi, how='left', left_on='PUBCHEM_CID', right_on='PUBCHEM_CID' )

# make column to store mol objects
print('generating rdkit mol objects from pubchem smiles')
#PandasTools.AddMoleculeColumnToFrame( df, 'smiles', 'rdkit_mol', includeFingerprints=False )
df['rdkit_mol'] = df['smiles'].apply( lambda x: Chem.MolFromSmiles(x) if x is not None else None )

# remove cpds that didn't process in rdkit (these were confirmed to be inorganics)
df = df.loc[ ~df['rdkit_mol'].isnull() ]

# generate rdkit smiles based on PubChem smiles (as salts)
print('getting rdkit canonicalized smiles')
df['rdkit_smiles'] = df['rdkit_mol'].apply( lambda x: Chem.MolToSmiles(x) if x is not None else None )

# generate de-salted (clean) rdkit smiles (rdkit_smiles_cln)
print('removing salt counterions')
_saltRemover = SaltRemover.SaltRemover()
df['rdkit_mol'] = df['rdkit_mol'].apply( lambda x: _saltRemover.StripMol(x, dontRemoveEverything=True) if x is not None else None)

# apply secondary salt remover
print('removing salt counterions missed by 1st saltremover')
df['rdkit_mol'] = df['rdkit_mol'].apply( lambda x: alt_frag_remover(x) if x is not None else None )

# standardize charges
#un = rdMolStandardize.Uncharger()
#df['rdkit_mol'] = df['rdkit_mol'].apply( lambda x: un.uncharge(x) if x is not None else None)

# replace smiles with clean smiles
print('getting rdkit canonicalized smiles')
df['rdkit_smiles_cln'] = df['rdkit_mol'].apply( lambda x: Chem.MolToSmiles(x) if x is not None else None )

# boolean for detection of substructure
#df['match'] = df[ df['rdkit_mol'] >= frag ]
#df['match'] = df['rdkit_mol'].apply( lambda x: x.HasSubstructMatch(frag) if x is not None else None )

# drop the mol objects columns if dumping to CSV
df.drop( columns=['rdkit_mol','PUBCHEM_RESULT_TAG'], inplace=True )

# get murcko smiles from rdkit_smiles
#print('getting murcko scaffolds')
#df['murcko'] = df['rdkit_smiles_cln'].apply( lambda x: get_murcko(x) if x is not None else None )
#print('getting generic murcko scaffolds')
#df['gen_murcko'] = df['murcko'].apply( lambda x: get_gen_murcko(x) if x is not None else None )
    
# dump to CSV
#df = df.replace(to_replace='None', value=np.nan).dropna( subset=['rdkit_mol'] )
print('dumping pickle')
df.to_pickle('oxphos_merged_aids_cids_clnsmi.pkl')


# dump to spreadsheet (with images):
#df = df.reset_index()
#df.fillna(value='', inplace=True)
#PandasTools.SaveXlsxFromFrame( df, 'aid_'+aid+'.xlsx', molCol='rdkit_mol', size=(200,200) )

