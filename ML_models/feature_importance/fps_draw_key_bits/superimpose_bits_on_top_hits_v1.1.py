
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import DataStructs
from rdkit.Chem import RDConfig
from rdkit.Chem import rdBase
# with CoordGen
from rdkit.Chem import rdCoordGen
print(rdBase.rdkitVersion)

# not sure if need this one
# https://iwatobipen.wordpress.com/2018/11/07/visualize-important-features-of-machine-leaning-rdkit/
def mol2fp(mol,nBits=1024):
    bitInfo={}
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, bitInfo=bitInfo)
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr, bitInfo

# https://rdkit.blogspot.com/2016/02/morgan-fingerprint-bit-statistics.html
# Functions for providing detailed descriptions of MFP bits from Nadine Schneider 
#  It's probably better to do this using the atomSymbols argument but this does work.

def includeRingMembership(s, n):
    r=';R]'
    d="]"
    return r.join([d.join(s.split(d)[:n]),d.join(s.split(d)[n:])])
 
def includeDegree(s, n, d):
    r=';D'+str(d)+']'
    d="]"
    return r.join([d.join(s.split(d)[:n]),d.join(s.split(d)[n:])])
 
def writePropsToSmiles(mol,smi,order):
    #finalsmi = copy.deepcopy(smi)
    finalsmi = smi
    for i,a in enumerate(order):
        atom = mol.GetAtomWithIdx(a)
        if atom.IsInRing():
            finalsmi = includeRingMembership(finalsmi, i+1)
        finalsmi = includeDegree(finalsmi, i+1, atom.GetDegree())
    return finalsmi
 
def getSubstructSmi(mol,atomID,radius):
    if radius>0:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol,radius,atomID)
        atomsToUse=[]
        for b in env:
            atomsToUse.append(mol.GetBondWithIdx(b).GetBeginAtomIdx())
            atomsToUse.append(mol.GetBondWithIdx(b).GetEndAtomIdx())
        atomsToUse = list(set(atomsToUse))
    else:
        atomsToUse = [atomID]
        env=None
    smi = Chem.MolFragmentToSmiles(mol,atomsToUse,bondsToUse=env,allHsExplicit=True, allBondsExplicit=True, rootedAtAtom=atomID)
    order = eval(mol.GetProp("_smilesAtomOutputOrder"))
    smi2 = writePropsToSmiles(mol,smi,order)
    return smi,smi2

# specificy feature type
feat_type = "fps"

# get important fps bit indices (top10?)
#df_fps_feat_imp = pd.read_csv('../MDI/fps/RF_nocalibr_fps_MDI.csv')
#df_fps_feat_imp = pd.read_csv('../permutation/fps/RF_calibr_fps_perm_avgprec.csv')
df_fps_feat_imp = pd.read_csv('../permutation/fps/RF_nocalibr_fps_perm_avgprec_5x.csv')
df_fps_feat_imp['feature'] = df_fps_feat_imp['feature'].astype('int')

# get top 20 most import fps bits by permutation method on RF no calibr model
important_bits = df_fps_feat_imp.head(40)['feature'].tolist()
# [1152, 1917, 41, 1816, 1060, 807, 464, 935, 1480, 1722]

# load fps data and RF scores on test set
# this test set was sorted by label and then RF score (best on top)
# note: there are but 255 true positives in this set of 23348 cpds
df_scores = pd.read_pickle('scores_RF_fps_test_set.pkl')
df_scores['score_percentile_rank'] = df_scores['score'].rank( pct=True )
df_scores['sum_key_bits'] = df_scores[important_bits].sum(axis=1)

# start with most important feature bit, look for first example
# in sorted positives for on-bit at that feature index
# draw molecule with that feature bit highlighted in molecule

from rdkit.Chem.Draw import rdMolDraw2D

#proto_molids = df_scores.loc[ (df_scores.label == 1) & (df_scores['sum_key_bits'] >= 5) ]['PUBCHEM_CID'].tolist()[:40]
proto_molids = df_scores.loc[ (df_scores.label == 0) & (df_scores['sum_key_bits'] >= 5) ]['PUBCHEM_CID'].tolist()[-40:]
#proto_molid = df_scores.loc[ (df_scores.label == 1) & (df_scores['sum_key_bits'] > 7) ]['PUBCHEM_CID'].iloc[0]

# loop through high priority cpds
#for molid in [ proto_molid ]:
for molid in proto_molids:
    match_mol = df_scores.loc[ df_scores['PUBCHEM_CID'] == molid ].iloc[0]
    mol = Chem.MolFromSmiles( match_mol['rdkit_smiles_cln']  )
    score = match_mol['score']
    label = match_mol['label']
    score_rank      = df_scores.loc[ df_scores['PUBCHEM_CID'] == molid ].index[0]
    score_perc_rank = df_scores.loc[ df_scores['PUBCHEM_CID'] == molid, 'score_percentile_rank' ].iloc[0]
    # get fingerprint and store bit info
    info = {}
    fp = AllChem.GetMorganFingerprintAsBitVect( mol, 2, 2048, bitInfo=info )
    rdCoordGen.AddCoords(mol)
    for rank, idx in enumerate(important_bits):
      try:
        # does compound have important bit on?
        if match_mol[idx] == 1:
            # if so, store atom id and radius
            match_frag in info[idx][0]:
            aid, rad = match_frag
            smi1, smi2 = getSubstructSmi( mol, aid, rad )
            patt = Chem.MolFromSmarts( smi2 )
            # highlightAtoms expects only one tuple, not tuple of tuples, so it needs to be flattened into single tuple 
            # store all occurrences of matched substructure as tuple of tuples (with atom indices) 
            matches = mol.GetSubstructMatches(patt)
            #matches = sum( mol.GetSubstructMatches(patt), () )
            hit_bonds = []
            for hit_ats in matches:
                #for bond in patt.GetBonds():
                    a1 = hit_ats[i][j]      
                    #a1 = hit_ats[ bond.GetBeginAtomIdx() ]
                    a2 = hit_ats[i][1]
                    #a2 = hit_ats[ bond.GetEndAtomIdx() ]
                    hit_bonds.append( mol.GetBondBetweenAtoms( a1, a2 ).GetIdx() )
            #        hit_bonds.append( mol.GetBondBetweenAtoms( a1, a2 ).GetIdx() )
            d = rdMolDraw2D.MolDraw2DCairo( 500, 500 )
            legend_text = "CID:{}, bitrnk:{} bit:{} scr:{:.3f}, scr_rank:{:.3f}, label:{}".format( molid, rank, idx, score, score_perc_rank, label)
            #rdCoordGen.AddCoords(mol)
            #rdMolDraw2D.PrepareAndDrawMolecule( d, mol, legend=legend_text, highlightAtoms=hit_ats, highlightBonds=hit_bonds )
            rdMolDraw2D.PrepareAndDrawMolecule( d, mol, legend=legend_text, highlightAtoms=matches )
            d.WriteDrawingText('CID'+str(molid)+'_rank'+str(rank).zfill(2)+'_bit'+str(idx)+'.png')
      except:
        print('{} rank:{} bit:{} failed'.format( str(molid), str(rank), str(idx) ) )
        pass

