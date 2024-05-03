
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools

try:
    thresh = str( sys.argv[1] )
except:
    print('')
    print('input threshold used as cutoff in clustering as only argument')
    print('')
    exit()


#thresh = "0.75"
inpkl = 'unique_pubchem_cids_all_wclusts.pkl'

df1 = pd.read_pickle( inpkl )
df1 = df1[['PUBCHEM_CID', 'clustid_'+thresh, 'clustmedoid_'+thresh, 'rdkit_smiles_cln']].dropna( subset=['clustid_'+thresh] )
#df1 = df1.drop_duplicates( subset=['PUBCHEM_CID'], keep='first' )

PandasTools.AddMoleculeColumnToFrame( df1, 'rdkit_smiles_cln', 'rdkit_mol', includeFingerprints=False )

for clustid in df1['clustid_'+thresh].unique():
    # for each cluster get list of mol objects, reorder with clustmedoid first
    ms_clustmedoid = df1.loc[ (df1['clustid_'+thresh] == clustid) & (df1['clustmedoid_'+thresh] == 1), 'rdkit_mol' ].to_list()
    ms        = df1.loc[ (df1['clustid_'+thresh] == clustid) & (df1['clustmedoid_'+thresh] == 0), 'rdkit_mol' ].to_list()
    ms = ms_clustmedoid + ms
    # get molecule identifiers
    m_names_clustmedoid = df1.loc[ (df1['clustid_'+thresh] == clustid) & (df1['clustmedoid_'+thresh] == 1), 'PUBCHEM_CID' ].to_list()
    m_names        = df1.loc[ (df1['clustid_'+thresh] == clustid) & (df1['clustmedoid_'+thresh] == 0), 'PUBCHEM_CID' ].to_list()
    m_names = m_names_clustmedoid + m_names
    cluster_pop = len(ms)
    print( "cluster_id:{} cluster_pop:{} len(ms):{} len(m_names):{}".format( clustid, cluster_pop, len(ms), len(m_names) ) )
    if cluster_pop == 1:
        img = Draw.MolsToGridImage( ms, molsPerRow=5, subImgSize=(400,400), legends=[ str(x) for x in m_names ] )
        img.save( "./clusters/"+thresh+"/cluster_pop"+str(cluster_pop).zfill(3)+"_index"+str(clustid).zfill(3)+".png")
    else:
        try:
            # if more than 25 cluster members, draw random sample
            if cluster_pop > 24:
                # obtain random draw of 24 members from cluster, add clustmedoid first in order (25 cpds drawn)
                rnd_idx = np.random.choice( cluster_pop - 1, 23, replace=False) + 1
                rnd_idx = [0] + list(rnd_idx)
                try:
                    ms = [ ms[i] for i in rnd_idx ]
                    m_names = [ m_names[i] for i in rnd_idx ]
                except:
                    print( "indices:{}".format( ",".join( [str(i) for i in rnd_idx] ) ) )
                    pass

            # find MCS in cluster (template), align all mols to this template
            mcs = rdFMCS.FindMCS(ms, ringMatchesRingOnly=True, completeRingsOnly=True, timeout=30)
            if len(mcs.smartsString) == 0:
                template = None
            else:
                template = Chem.MolFromSmarts( mcs.smartsString )
                AllChem.Compute2DCoords( template )
                for m in ms:
                    AllChem.GenerateDepictionMatching2DStructure( m, template )
            print( "template for cluster:"+str(clustid), mcs.smartsString )

            img = Draw.MolsToGridImage( [template] + ms, molsPerRow=5, 
                                        subImgSize=(400,400),
                                        legends=[ str(x) for x in ['mcs_template'] + m_names ]  
                                      )

            img.save( "./clusters/"+thresh+"/cluster_pop"+str(cluster_pop).zfill(3)+"_index"+str(clustid).zfill(3)+".png")
        except:
            m_smiles = df1.loc[ df1['clustid_'+thresh] == clustid, 'rdkit_smiles_cln' ].to_list()
            m_names =  df1.loc[ df1['clustid_'+thresh] == clustid, 'PUBCHEM_CID' ].to_list()
            m_names_smiles = zip( m_names, m_smiles )
            print( "cluster_id:{} bombed".format( clustid ) )
            for x in m_names_smiles:
                print( "\t{}:{}".format( str(x[0]), str(x[1]) ) )
            print('')
            pass

