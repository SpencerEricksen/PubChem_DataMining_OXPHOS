
import numpy as np
import pandas as pd
import pickle
import scipy.stats as st
from matplotlib import pyplot as plt
import seaborn as sns

df = pd.read_pickle('../df_umap_nn50_md0.250_final.pkl')

# get UMAP coords for actives and inactives
#df_active = df.loc[ df['label'] == 1 ][['X','Y']]
df_inactive = df.loc[ df['label'] == 0 ][['X','Y']]

# get the boundaries
xmin = df['X'].min()
xmax = df['X'].max()
ymin = df['Y'].min()
ymax = df['Y'].max()

# compute KDE object for active point distribution in 2D UMAP space
#xy_active = np.vstack( (df_active.X.values, df_active.Y.values) )
#kde_active = st.gaussian_kde( xy_active )

# compute KDE object for inactive point distribution in 2D UMAP space
xy_inactive = np.vstack( (df_inactive.X.values, df_inactive.Y.values) )
kde_inactive = st.gaussian_kde( xy_inactive, bw_method=0.136 )

# build a matrix to store bin values for each x,y bin (500 bins along each axis for 250,000 bins total)
gx, gy = np.mgrid[ xmin:xmax:500j, ymin:ymax:500j ]
gxy = np.dstack((gx, gy))

# get grid of values in 500x500 bins
#z_active = np.apply_along_axis( kde_active, 2, gxy )
#z_active = z_active.reshape(500,500)

z_inactive = np.apply_along_axis( kde_inactive, 2, gxy )
z_inactive = z_inactive.reshape(500,500)

# get 2D difference map
#z_diff = z_active - z_inactive

# plot heatmap
import matplotlib.colors as colors
#import matplotlib.cbook as cbook
from matplotlib import cm

fig, ax = plt.subplots(1,1)
cmap='copper_r'
#cmap='flare'
ax1 = sns.heatmap( z_inactive, cmap=cmap, xticklabels=False, yticklabels=False )
plt.savefig( 'pp_inactive.png', dpi=2400 )

