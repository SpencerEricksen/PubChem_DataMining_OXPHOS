#!/home/ssericksen/anaconda2/envs/py36_chem/bin/python

import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import rdMolDraw2D

#from rdkit.Chem import AllChem
#import rdkit.Chem.Draw
#from rdkit.Chem.Draw import rdMolDraw2D
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

try:
    thresh = sys.argv[1]
    incsv = sys.argv[2]
except:
    print('')
    print('need thresh and incsv as arguments')
    print('')
    exit()


df1 = pd.read_csv( incsv, sep="|")
df2 = df1[df1['medoid']==1]

df3 = df1.drop_duplicates( subset=['PUBCHEM_CID'], keep='first' )
clst_id_list = df3['cluster_id'].value_counts().index.tolist()

ms = []
ms_titles = []
#for clst_id in sorted( df2['cluster_id'].unique().tolist() ):
for clst_id in clst_id_list:
    clst_pop  = len(df3[df3['cluster_id']==clst_id])
    pmid_list = df2[df2['cluster_id']==clst_id]['pmid'].unique().tolist()
    aid_list  = df2[df2['cluster_id']==clst_id]['AID'].unique().tolist()
    cid       = df2[df2['cluster_id']==clst_id]['PUBCHEM_CID'].iloc[0]
    smiles    = df2[df2['cluster_id']==clst_id]['smiles'].iloc[0]
    ms.append( Chem.MolFromSmiles( smiles ) )
    pmid_string = " ".join( [ str(i) for i in pmid_list ] )
    aid_string = " ".join( [str(i) for i in aid_list ] )
    wrapped_string = "cluster:"+str(clst_id)+"   members:"+str(clst_pop)+"   CID:"+str(cid)
    ms_titles.append( wrapped_string )
    print( "cluster:{}, members:{}, cid:{}, smiles:{}, PMIDs:{}, AIDs:{}".format( clst_id, clst_pop, cid, smiles, pmid_string, aid_string) )

#img = Draw.MolsToGridImage( ms, molsPerRow=7, subImgSize=(600,600), legends=ms_titles )
img = DrawMolsZoomed( ms, molsPerRow=20, subImgSize=(500,500), legends=ms_titles )
img.save( "cluster_reps_"+str(thresh)+".png")

