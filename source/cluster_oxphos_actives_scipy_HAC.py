
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from itertools import combinations


def get_medoids( clusters, dist_mat, labels, dist_thresh, n_mols ):
    '''get cluster medoids using distance matrx, and
       return '''
    dist_dict = {}
    for cid in labels:
        member_mols = np.where( clusters == cid )[0]
        cid_dist_matrix = np.zeros( (n_mols, n_mols) )
        for pair in combinations( member_mols, 2 ):
            cid_dist_matrix[ pair[0], pair[1] ] = dist_mat[ pair[0], pair[1] ]
        cid_dist_sums = ( cid_dist_matrix.sum( axis=0 ) + cid_dist_matrix.sum( axis=1 ) ) / float(len(member_mols))
        for mol, dist_sum in zip( member_mols, cid_dist_sums[member_mols] ):
            # store clustering info for each mol
            if cid not in dist_dict:
                dist_dict[cid] = [ (mol, dist_sum) ]
            else:
                dist_dict[cid].append( (mol, dist_sum) )
        dist_dict[cid] = sorted( dist_dict[cid], key=lambda x:x[1] )
    # find representative cpd for each cid
    rep_mols = [ dist_dict[cid][0][0] for cid in dist_dict ]
    return rep_mols


def make_bitstring(fp):
    try:
        return "".join( [ str(i) for i in fp ] )
    except:
        return None


# get active cpds
df1 = pd.read_pickle('unique_pubchem_cids_all.pkl')
#df2 = df1.loc[ (df1.label == 1) & (df1.num_etc_linked_aids > 0), ['PUBCHEM_CID','ECFP6_bitstr'] ]
#df2 = df1.loc[ (df1.label == 1) & (df1.num_etc_linked_aids > 0) \
#                                & (df1.lipinski == 1), ['PUBCHEM_CID','ECFP6_bitstr'] ]
df2 = df1.loc[ (df1.label == 1) & (df1.num_etc_linked_aids > 0) \
                                & (df1.lipinski == 1) & (df1.pains == False), \
                                ['PUBCHEM_CID','ECFP6_bitstr'] ]
# get fps
df2['morgan'] = df2['ECFP6_bitstr'].apply( lambda x: [int(i) for i in list(x) ] )
X = np.stack( df2['morgan'].values )
# get PUBCHEM_CIDs
molid_list = df2['PUBCHEM_CID'].tolist()

# generate fingerprint matrix for top 2000 cpds
dists = pdist( X, metric='jaccard' )
# get distance matrix for medoid determination
dist_mat = squareform( dists )
n_mols = np.shape(X)[0]

# run HAC
Z = linkage( dists, method='average', optimal_ordering=True )

# run clustering slices through dendrogram at various thresholds
for dist_thresh in [ 0.700, 0.725, 0.750, 0.775, 0.800 ]:
    clusters = fcluster( Z, dist_thresh, criterion='distance' )
    # get set of cluster IDs
    labels = np.unique(clusters).tolist()
    d_thr_str = str( '{:.3f}'.format( dist_thresh ) )
    df2['clustid_'+d_thr_str] = clusters
    rep_mols = get_medoids( clusters, dist_mat, labels, dist_thresh, n_mols )
    rep_mols_pccids = [ molid_list[i] for i in rep_mols ]
    df2['clustmedoid_'+d_thr_str] = 0
    df2.loc[ df2['PUBCHEM_CID'].isin( rep_mols_pccids ), 'clustmedoid_'+d_thr_str ]  = 1
    df2[ 'clustpop_'+d_thr_str ] = df2.groupby( 'clustid_'+d_thr_str )['PUBCHEM_CID'].transform('count')

df2.drop( columns=['morgan','ECFP6_bitstr'], inplace=True )
df3 = df1.merge( df2, on='PUBCHEM_CID', how='left')
df3.to_pickle('unique_pubchem_cids_all_wclusts.pkl')

