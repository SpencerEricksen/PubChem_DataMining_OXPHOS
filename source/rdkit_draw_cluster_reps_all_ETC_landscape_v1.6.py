
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
inpkl = 'unique_pubchem_cids_all_wclusts.pkl'
df1 = pd.read_pickle( inpkl )
#df2 = df1.loc[ (df1['clustmedoid_0.750'] == 1) & (df1['clustpop_0.750'] >= 3.0) ]
df2 = df1.loc[ df1['clustmedoid_0.750'] == 1 ]
df2 = df2[['PUBCHEM_CID','rdkit_smiles_cln','clustid_0.750','clustpop_0.750']]
df2['clustid_0.750'] = df2['clustid_0.750'].astype('int')
df2['clustpop_0.750'] = df2['clustpop_0.750'].astype('int')

clst_id_list = df2['clustid_0.750'].sort_values().tolist()
ms = []
ms_titles = []

for clst_id in clst_id_list:
    clst_pop = df2.loc[ df2['clustid_0.750'] == clst_id, 'clustpop_0.750' ].iloc[0]
    cid = df2.loc[ df2['clustid_0.750'] == clst_id, 'PUBCHEM_CID' ].iloc[0]
    smiles = df2.loc[ df2['clustid_0.750'] == clst_id, 'rdkit_smiles_cln' ].iloc[0]
    ms.append( Chem.MolFromSmiles( smiles ) )
    wrapped_string = "clustID: "+str(clst_id)+"  size: "+str(clst_pop)+"  CID: "+str(cid)
    ms_titles.append( wrapped_string )
    print( "clust_id:{}, clust_pop:{}, pubchem_CID:{}, smiles:{}".format( clst_id, clst_pop, cid, smiles) )

#img = Draw.MolsToGridImage( ms, molsPerRow=7, subImgSize=(600,600), legends=ms_titles )
img = DrawMolsZoomed( ms, molsPerRow=26, subImgSize=(500,500), legends=ms_titles )
img.save( "cluster_medoids_all_0.750_landscape.png")

