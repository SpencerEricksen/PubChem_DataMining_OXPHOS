
import numpy as np
import pandas as pd
import pickle
import scipy.stats as st
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors as colors
#import matplotlib.cbook as cbook
from matplotlib import cm
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

z_diff = np.load( "umap_kde_density_matrix_diff_bw0.100.npy" )
z_active = np.load( "umap_kde_density_matrix_active_bw0.100.npy" )
z_inactive = np.load( "umap_kde_density_matrix_inactive_bw0.100.npy" )

# plot diff map (active - inactive)
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
df_labels = df.loc[ df['clustid_0.750'].isin( label_list ) ]
df_medoids = df_labels[ df_labels['clustmedoid_0.750'] == 1 ]
kde_cmap = 'Spectral_r'
#cmap = cm.get_cmap("jet").copy()
#cmap = 'PiYG' #'seismic'
#levels = np.linspace( -0.040, 0.040, 15)
levels = np.linspace( -0.025, 0.050, 15)
umap_ticks = np.linspace( -13.5, 13.5, 500 )
pc = ax1.contourf( umap_ticks, umap_ticks, z_diff, cmap=kde_cmap, levels=levels, extend='max' )
fig.colorbar(pc)
pc.cmap.set_under( color='w')
# plot cluster medoid points for clusters with populations >= 3
ax1.scatter( df_medoids['X'], df_medoids['Y'], c=df_medoids['clustid_0.750'], cmap=cmap, alpha=0.7, lw=0.1, s=1.0, edgecolors='black' )
ax1.scatter( df_6['X'], df_6['Y'], c=df_6['clustid_0.750'], cmap=cmap, alpha=0.7, lw=0.2, s=3.0, edgecolors='black' )
#ax1.scatter( df_6['X'], df_6['Y'], 
'''
plt.tick_params(
    axis='both',
    which='both',
    top=False,
    bottom=False,
    left=False,
    labelleft=False,
    labelbottom=False)
'''
for m in df_medoids['PUBCHEM_CID'].to_list():
    x = df_medoids.loc[ df_medoids['PUBCHEM_CID'] == m, 'X' ]
    y = df_medoids.loc[ df_medoids['PUBCHEM_CID'] == m, 'Y' ]
    cid = df_medoids.loc[ df_medoids['PUBCHEM_CID'] == m, 'clustid_0.750' ]
    cid = str(int(cid))
    #ax1.annotate( cid, ( x, y ), size=4.0, rotation=10, ha='left', va='top' )
for m in df_6['PUBCHEM_CID'].to_list():
    x = df_6.loc[ df_6['PUBCHEM_CID'] == m, 'X' ]
    y = df_6.loc[ df_6['PUBCHEM_CID'] == m, 'Y' ]
    compound_name = df_6.loc[ df_6['PUBCHEM_CID'] == m, 'compound_name' ].iloc[0]
    #compound_number = df_6.loc[ df_6['PUBCHEM_CID'] == m, 'cpd_number' ].iloc[0]
    #compound_label = str(compound_name)+" - cpd "+str(compound_number) 
    ax1.annotate( compound_name, ( x, y ), size=6.0, rotation=0, ha='left', va='top', xytext = ( x + 2.75, y + 4.75 ), arrowprops=dict(arrowstyle='-', linewidth=0.5, connectionstyle='arc3,rad=0.0') )

plt.xlabel('umap-1')
plt.ylabel('umap-2')
plt.xlim(-7.0,8.0)
plt.ylim(-13.5,10.5)
#plt.grid()
plt.savefig( 'fig3_difference_map_bw100_lv15_spectralr_contour_wpts_6cpds_nogrid_wbar.png', dpi=2400 )
plt.close()

'''
# plot active KDE
fig, ax1 = plt.subplots(1,1)
#cmap = cm.get_cmap("pink_r").copy()
levels = np.linspace( 0.00000, 0.050, 15)
pc = ax1.contourf( umap_ticks, umap_ticks, z_active, levels=levels, cmap=kde_cmap, extend='max' )
ax1.scatter( df_medoids['X'], df_medoids['Y'], c=df_medoids['clustid_0.750'], cmap=cmap, alpha=0.7, lw=0.1, s=0.2, edgecolors='black' )
ax1.scatter( df_6['X'], df_6['Y'], c=df_6['clustid_0.750'], cmap=cmap, alpha=0.7, lw=0.2, s=1.4, edgecolors='black' )
plt.tick_params(
    axis='both',
    which='both',
    top=False,
    bottom=False,
    left=False,
    labelleft=False,
    labelbottom=False)
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
#pc.cmap.set_under( color='w')
#plt.grid()
#fig.colorbar(pc)
plt.savefig( 'actives_15levels_bw100_spectralr_ngrid_nbar.png', dpi=2400 )
plt.close()

# plot inactive KDE
fig, ax1 = plt.subplots(1,1)
#cmap = cm.get_cmap("pink_r").copy()
levels = np.linspace( 0.00000, 0.0250, 15)
pc = ax1.contourf( umap_ticks, umap_ticks, z_inactive, levels=levels, cmap=kde_cmap, extend='max' )
ax1.scatter( df_medoids['X'], df_medoids['Y'], c=df_medoids['clustid_0.750'], cmap=cmap, alpha=0.7, lw=0.1, s=0.2, edgecolors='black' )
ax1.scatter( df_6['X'], df_6['Y'], c=df_6['clustid_0.750'], cmap=cmap, alpha=0.7, lw=0.2, s=1.4, edgecolors='black' )
plt.tick_params(
    axis='both',
    which='both',
    top=False,
    bottom=False,
    left=False,
    labelleft=False,
    labelbottom=False)
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
#pc.cmap.set_under( color='w')
#plt.grid()
#fig.colorbar(pc)
plt.savefig( 'inactives_15levels_bw100_spectralr_ngrid_nbar.png', dpi=2400 )
plt.close()
'''
