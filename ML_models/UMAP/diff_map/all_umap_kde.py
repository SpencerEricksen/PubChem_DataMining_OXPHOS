
import numpy as np
import pandas as pd
import pickle
import scipy.stats as st
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
#import matplotlib.cbook as cbook
import seaborn as sns

df = pd.read_pickle('../df_umap_nn50_md0.250_final.pkl')

# get UMAP coords for actives and inactives
#df_active = df.loc[ df['label'] == 1 ][['X','Y']]
#df_inactive = df.loc[ df['label'] == 0 ][['X','Y']]

df_all = df[['X','Y']]

# get the boundaries
xmin = df['X'].min()
xmax = df['X'].max()
ymin = df['Y'].min()
ymax = df['Y'].max()

xmin = -13.5
xmax = 13.5
ymin = -13.5
ymax = 13.5

# compute KDE object for full point distribution in 2D UMAP space
xy_all = np.vstack( (df_all.X.values, df_all.Y.values) )
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

np.save( "umap_kde_density_matrix_all.npy", z_all )
#z_all = np.load( "umap_kde_density_matrix_all.npy" )

# plot heatmap
fig, ax = plt.subplots(1,1)
#cmap='copper_r'
cmap = sns.color_palette("Blues_d", as_cmap=True )
#cmap = cm.get_cmap("Blues").copy()
#cmap = cm.get_cmap("YlGnBu").copy()
#cmap='flare'
#pc = sns.heatmap(z_all, ax=ax, cmap=cmap, xticklabels=False, yticklabels=False)
levels = np.linspace( 0.00001, 0.025, 15)
pc = plt.contourf( z_all, levels=levels, cmap=cmap, extend='max' )
pc.cmap.set_under( color='w')
#plt.axis('off')
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
plt.grid()
fig.colorbar(pc, ax=ax)
# could provide xticklabels as list ( np.arange(-13.5, 14.5, 1.5) )
plt.savefig( 'pp_all_rot90_14levels_bw100_BluesD_wgrid.png', dpi=2400 )
plt.close()

