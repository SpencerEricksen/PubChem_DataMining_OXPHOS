
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import rdMolDraw2D

try:
    import Image
except ImportError:
    from PIL import Image
from io import BytesIO

def DrawMolsZoomed( mols, molsPerRow=3, subImgSize=(200, 200), legends=None ):
    nRows = len(mols) // molsPerRow
    if len(mols) % molsPerRow: nRows += 1
    fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
    full_image = Image.new('RGBA', fullSize )
    for ii, mol in enumerate(mols):
        if mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(mol)
        column = ii % molsPerRow
        row = ii // molsPerRow
        offset = ( column * subImgSize[0], row * subImgSize[1] )
        d2d = rdMolDraw2D.MolDraw2DCairo(subImgSize[0], subImgSize[1] )
        d2d.drawOptions().legendFontSize=20
        d2d.DrawMolecule(mol, legend=legends[ii] )
        d2d.FinishDrawing()
        sub = Image.open(BytesIO(d2d.GetDrawingText()))
        full_image.paste(sub,box=offset)
    return full_image

# read in 'unique_pubchem_cids_all_wclusts.pkl'
incsv = sys.argv[1]
#incsv = './ML_models/UMAP/diff_map/hitlist_6_cpds_umap_smi.csv'
df1 = pd.read_csv( incsv )
df1['clustid_0.750'] = df1['clustid_0.750'].astype('int')
df1['clustpop_0.750'] = df1['clustpop_0.750'].astype('int')

ms = []
ms_titles = []

for pccid in df1['PUBCHEM_CID'].tolist():
    clst_pop = df1.loc[ df1['PUBCHEM_CID'] == pccid, 'clustpop_0.750'  ].iloc[0]
    clst_id =  df1.loc[ df1['PUBCHEM_CID'] == pccid, 'clustid_0.750'   ].iloc[0]
    cid     =  df1.loc[ df1['PUBCHEM_CID'] == pccid, 'PUBCHEM_CID'     ].iloc[0]
    smiles =   df1.loc[ df1['PUBCHEM_CID'] == pccid, 'rdkit_smiles_cln'].iloc[0]
    ms.append( Chem.MolFromSmiles( smiles ) )
    wrapped_string = "clustID: "+str(clst_id)+"  size: "+str(clst_pop)+"  CID: "+str(cid)
    ms_titles.append( wrapped_string )
    print( "clust_id:{}, clust_pop:{}, pubchem_CID:{}, smiles:{}".format( clst_id, clst_pop, cid, smiles) )

#img = Draw.MolsToGridImage( ms, molsPerRow=7, subImgSize=(600,600), legends=ms_titles )
img = DrawMolsZoomed( ms, molsPerRow=6, subImgSize=(500,500), legends=ms_titles )
img.save( "gridmols_6cpds_landscape.png")


