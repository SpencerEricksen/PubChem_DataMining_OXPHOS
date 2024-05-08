
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, QuantileTransformer, RobustScaler
import matplotlib.pyplot as plt
import sys

try:
    scaler = sys.argv[1] # 'rb', 'qt', 'st'
    stat = sys.argv[2]  # 'mean' or 'median'
except:
    print('')
    print('provide scaler (rb, qt, or st) as first argument')
    print('provide statistic (mean or median) as 2nd argument')
    print('')
    exit()


df2 = pd.read_pickle('../../../ML_training_data/data_set_desc.pkl')

#df2 = df1.loc[ (df1.label == 1) & (df1.num_etc_linked_aids > 0) \
#               & (df1.lipinski == 1) & (df1.pains == False), \
#               ['PUBCHEM_CID','ECFP6_bitstr'] ]
#df3 = df1.loc[ (df1.label == 0) & (df1.num_etc_linked_aids > 0) & (df1.lipinski == 1) & (df1.pains == False) ]

# bad feature
df2.drop( columns=['desc_Ipc'], inplace=True )

# get only fragment based descritpors
df2.drop( columns=[ i for i in df2.columns[2:] if 'fr_' not in i ], inplace=True )

# choose feature scaler
if scaler == 'rb':
    transformer = RobustScaler()
elif scaler == 'qt':
    transformer = QuantileTransformer()
elif scaler == 'st':
    transformer = StandardScaler()

df2[ df2.columns[2:] ] = transformer.fit_transform(df2[ df2.columns[2:]] )

# sort descriptors by differences between active and inactive instances
if stat == 'mean':
    m_0 = df2.loc[ df2['label'] == 0, df2.columns[2:] ].mean()
    m_1 = df2.loc[ df2['label'] == 1, df2.columns[2:] ].mean()
elif stat == 'median':
    m_0 = df2.loc[ df2['label'] == 0, df2.columns[2:] ].median()
    m_1 = df2.loc[ df2['label'] == 1, df2.columns[2:] ].median()

df_diff = pd.concat( [m_0, m_1], axis=1 )
df_diff.columns = ['m_0', 'm_1']
df_diff['diff'] =  df_diff['m_1'] - df_diff['m_0']

df_diff['colors'] = [ 'red' if x < 0 else 'green' for x in df_diff['diff'] ]
df_diff.sort_values( by='diff', inplace=True )
df_diff.reset_index(inplace=True)
df_diff.rename( columns={'index':'frag'}, inplace=True )

# dump frag difference data to CSV
df_diff.to_csv('diff_data_desc_frags_mean_st.csv', index=False )

# get only most divergent features for plotting
df_diff2 = pd.concat( [ df_diff.iloc[0:20], df_diff.iloc[-20:] ] )
frag_list = [ i[8:] for i in df_diff2['frag'].tolist() ]
df_diff2 = df_diff2.reset_index()
df_diff2.drop( columns='index', inplace=True )

# Draw plot
plt.figure( figsize=(5,14.5) )
plt.hlines( y=df_diff2.index.tolist(), xmin=0, xmax=df_diff2['diff'].tolist(), color='black' )
for x, y, tex in zip( df_diff2['diff'].tolist(), df_diff2.index.tolist(), df_diff2['diff'].tolist() ):
    t = plt.text(x, y, str(round(tex, 2)), fontsize=12, horizontalalignment='right' if x < 0 else 'left',
                 verticalalignment='center', fontdict={'color':'red' if x < 0 else 'green', 'size':5})
# Decorations
plt.yticks( df_diff2.index.tolist(), frag_list, fontsize=12)
#plt.title('Feature Div, Mean/StandardScaler', fontdict={'size':12} )
plt.title('Relative Abundance of Substructures in Actives', fontdict={'size':12} )
plt.grid(linestyle='--', alpha=0.5)
plt.xlim(-1.25, 1.25)
plt.savefig( 'Fig2B_diff_bars_frags_'+scaler+'_'+stat+'.png', bbox_inches='tight', dpi=600 )

