
import numpy as np
import pandas as pd
import pickle
import scipy.stats as st
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
#import matplotlib.cbook as cbook
import seaborn as sns

df = pd.read_pickle('../df_umap_nn50_md0.250_final.pkl')

df[['clustid_0.750','clustpop_0.750','clustmedoid_0.750']] = df[['clustid_0.750','clustpop_0.750','clustmedoid_0.750']].astype('int')

molid_list = df.PUBCHEM_CID.tolist()
clustid_list = df['clustid_0.750'].tolist()
# top 135 clusters (triplets or larger)
label_list =  df[ df['clustpop_0.750'] > 2 ]['clustid_0.750'].value_counts().index.tolist()
num_clusters = len(label_list)


# get 6 cpd candidate info
df_6 = pd.read_csv('hitlist_6_cpds_umap.csv')

# get the boundaries
xmin = -13.5
xmax = 13.5
ymin = -13.5
ymax = 13.5

'''
# compute KDE object for full point distribution in 2D UMAP space
xy_all = np.vstack( (df.X.values, df.Y.values) )
#kde_all = st.gaussian_kde( xy_all, bw_method=0.136 )
kde_all = st.gaussian_kde( xy_all, bw_method=0.100 )

# build a matrix to store bin values for each x,y bin (500 bins along each axis for 250,000 bins total)
gx, gy = np.mgrid[ xmin:xmax:500j, ymin:ymax:500j ]
gxy = np.dstack((gx, gy))

# get grid of values in 500x500 bins
z_all = np.apply_along_axis( kde_all, 2, gxy )
z_all = z_all.reshape(500,500)

# flip array vertically and horizontally
#z_all = np.flip( z_all, axis=0 )
#z_all = np.flip( z_all, axis=1 )
z_all = np.rot90( z_all )
z_all = np.flip( z_all, axis=0 )
'''
#np.save( "umap_kde_density_matrix_all.npy", z_all )
z_all = np.load( "umap_kde_density_matrix_all.npy" )

# plot heatmap
fig = plt.figure()

# define colormap for scatter points based on clusters
cmap = plt.cm.hsv
# extract all colors from the .jet map
cmaplist = [ cmap(i) for i in range(cmap.N) ]
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.linspace( 0, num_clusters, num_clusters + 1 )
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

ax1 = fig.add_subplot(111)

df_labels = df[ df['clustid_0.750'].isin( label_list ) ]

df_labels = df.loc[ df['clustid_0.750'].isin( label_list ) ]
df_medoids = df_labels[ df_labels['clustmedoid_0.750'] == 1 ]
#cmap='copper_r'
#kde_cmap = sns.color_palette("Blues_d", as_cmap=True )
#kde_cmap = 'Spectral_r'
kde_cmap = cm.get_cmap("Spectral_r").copy()
#cmap = cm.get_cmap("Blues").copy()
#cmap = cm.get_cmap("YlGnBu").copy()
#cmap='flare'
#pc = sns.heatmap(z_all, ax=ax, cmap=cmap, xticklabels=False, yticklabels=False)
levels = np.linspace( 0.000, 0.025, 20)
# need to provide x, y axis values for z_all
p = np.linspace(-13.5, 13.5, 500)
pc = ax1.contourf( p, p, z_all, levels=levels, cmap=kde_cmap, extend='max' )
#fig.colorbar(pc)
#pc.cmap.set_under( color='w')
#plt.axis('off')
plt.tick_params(
    axis='both',
    which='both',
    top=False,
    bottom=False,
    left=False,
    labelleft=False,
    labelbottom=False)
ax1.scatter( df.loc[ df['label'] == 1, 'X' ], df.loc[ df['label'] == 1, 'Y' ], c='lightgrey', alpha=1.0, lw=0.1, s=0.1, edgecolors='black')
ax1.scatter( df_labels['X'], df_labels['Y'], c=df_labels['clustid_0.750'], cmap=cmap, alpha=1.0, lw=0.1, s=0.2, edgecolors='black' )
ax1.scatter( df_6['X'], df_6['Y'], c=df_6['clustid_0.750'], cmap=cmap, alpha=1.0, lw=0.2, s=1.4, edgecolors='black' )
#ax1.scatter( df_6['X'], df_6['Y'], c=df_6['clustid_0.750'], facecolors='white', lw=0.2, s=1.4, edgecolors='black' )

#ax1.set_xlabel('umap-0')
#ax1.set_ylabel('umap-1')
#ax1.set_xlim([-13.5, 13.5])
#ax1.set_ylim([-13.5, 13.5])
for m in df_medoids['PUBCHEM_CID'].to_list():
    x = df_medoids.loc[ df_medoids['PUBCHEM_CID'] == m, 'X' ]
    y = df_medoids.loc[ df_medoids['PUBCHEM_CID'] == m, 'Y' ]
    cid = df_medoids.loc[ df_medoids['PUBCHEM_CID'] == m, 'clustid_0.750' ]
    cid = str(int(cid))
    ax1.annotate( cid, ( x, y ), size=0.5, rotation=10, ha='left', va='top' )
for m in df_6['PUBCHEM_CID'].to_list():
    x = df_6.loc[ df_6['PUBCHEM_CID'] == m, 'X' ]
    y = df_6.loc[ df_6['PUBCHEM_CID'] == m, 'Y' ]
    compound_name = df_6.loc[ df_6['PUBCHEM_CID'] == m, 'compound_name' ].iloc[0]
    #compound_number = df_6.loc[ df_6['PUBCHEM_CID'] == m, 'cpd_number' ].iloc[0]
    #compound_label = str(compound_name)+" - cpd "+str(compound_number) 
    ax1.annotate( compound_name, ( x, y ), size=2.2, rotation=0, ha='left', va='top', xytext = ( x + 0.15, y + 0.25 ) )
#plt.grid()
# could provide xticklabels as list ( np.arange(-13.5, 14.5, 1.5) )
plt.savefig( 'all_rot90_20levs_bw100_SpectralR_ngrid_wpts_6cpds.png', dpi=2400 )
plt.close()

