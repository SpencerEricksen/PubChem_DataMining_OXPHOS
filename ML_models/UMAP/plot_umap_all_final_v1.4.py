
import time as time
import pandas as pd
import numpy as np
import umap
import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt
#import joblib
import pickle

print('loading fingerprint data...')
# get fingerprint data with activity labels
df1 = pd.read_pickle('../../ML_training_data/data_set_fps.pkl')

# get cluster info
print('loading cluster data...')
df_clust = pd.read_pickle('../../unique_pubchem_cids_all_wclusts.pkl')
df_clust = df_clust[['PUBCHEM_CID','clustid_0.750','clustpop_0.750','clustmedoid_0.750']]

# add cluster data to cpd features (fingerprints)
df2 = df1.merge( df_clust, how='left', on='PUBCHEM_CID' )

# fill in missing cluster information with 0
df2 = df2.fillna( 0 )

# get molids
molid_list = df2.PUBCHEM_CID.tolist()
clustid_list = df2['clustid_0.750'].tolist()
# top 135 clusters (triplets or larger)
label_list = df2['clustid_0.750'].value_counts().head(135).index.tolist()
# remove cpds with null clustid (=0)
label_list.remove(0.0)
num_clusters = len(set( label_list ))

# get cpds (points), represented by 2048 bit-length fingerprints
X_fps = df2[ df2.columns[2:2050] ].values


def build_umap( X_fps, df2, molid_list, nn, md ):
    '''embed coords in umap space (2D) for easier visualization'''
    mapper = umap.UMAP( n_components=2, n_neighbors=nn, min_dist=md, metric='jaccard', random_state=42 ).fit( X_fps )
    X_umap = mapper.transform( X_fps )
    umap_df = pd.DataFrame( X_umap, columns=['X','Y'], index=molid_list )
    umap_df = umap_df.reset_index()
    umap_df.rename( columns={'index':'PUBCHEM_CID'}, inplace=True )
    umap_df = umap_df.merge( df2[['PUBCHEM_CID', 'label', 'clustid_0.750', 'clustpop_0.750', 'clustmedoid_0.750']], on='PUBCHEM_CID', how='left' )
    # save umap transformer for later use
    fn = "umap_transformer_nn50_md0250_random42.pkl"
    pickle.dump( mapper, open(fn, 'wb') )
    #joblib.dump( mapper, fn )
    # loaded_umap = pickle.load( (open(fn, 'rb')) )
    return umap_df


def plot_umap( umap_df, label_list, nn, md ):
    ''' PLOT cpds and color according to cluster IDs '''
    fig = plt.figure()
    # define colormap
    cmap = plt.cm.hsv
    # extract all colors from the .jet map
    cmaplist = [ cmap(i) for i in range(cmap.N) ]
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = np.linspace( 0, num_clusters, num_clusters + 1 )
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax1 = fig.add_subplot(111)
    ax1.grid(False)
    ax1.set_title( 'UMAP chem space nn:{} md:{:.3f}'.format( nn, md) )
    df_labels = umap_df.loc[ umap_df['clustid_0.750'].isin( label_list ) ]
    df_medoids = df_labels[ df_labels['clustmedoid_0.750'] == 1.0 ]
    my_cmap = sns.color_palette("Blues_d", as_cmap=True )
    sns.kdeplot( x=umap_df['X'], y=umap_df['Y'], cmap=my_cmap, shade=True, bw_adjust=0.75, gridsize=500, levels=10, cut=0, thresh=0.01, ax=ax1 )
    ax1.scatter( umap_df.loc[ umap_df['label'] == 1, 'X' ], umap_df.loc[ umap_df['label'] == 1, 'Y' ], c='lightgrey', alpha=0.7, lw=0.1, s=0.1, edgecolors='black')
    ax1.scatter( df_labels['X'], df_labels['Y'], c=df_labels['clustid_0.750'], cmap=cmap, alpha=0.7, lw=0.1, s=0.2, edgecolors='black' )
    ax1.set_xlabel('umap-0')
    ax1.set_ylabel('umap-1')
    ax1.set_xlim([-13.5, 13.5])
    ax1.set_ylim([-13.5, 13.5])
    # label medoids from top 135 clusters
    for m in df_medoids['PUBCHEM_CID'].to_list():
        x = df_medoids.loc[ df_medoids['PUBCHEM_CID'] == m, 'X' ]
        y = df_medoids.loc[ df_medoids['PUBCHEM_CID'] == m, 'Y' ]
        cid = df_medoids.loc[ df_medoids['PUBCHEM_CID'] == m, 'clustid_0.750' ]
        cid = str(int(cid))
        ax1.annotate( cid, ( x, y ), size=0.5, rotation=10, ha='left', va='top' )
    plt.savefig( 'umap_nn{}_md{:.3f}_final.png'.format(nn, md), bbox_inches='tight', dpi=3200 )
    plt.close()
    return


# scan umap parameters, dump plots

#for nn in [ 5, 10, 15, 20, 25, 30, 40, 50, 100, 200 ]:
for nn in [ 50 ]:
    #for md in [ 0.001, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50 ]: 
    for md in [0.25]:
        umap_df = build_umap( X_fps, df2, molid_list, nn, md )
        plot_umap( umap_df, label_list, nn, md )
        umap_df.to_pickle( 'df_umap_nn{}_md{:.3f}_final.pkl'.format( nn, md) )
        print( 'dumping umap_nn{}_md{:.3f}_final42_2.png'.format(nn, md) )




