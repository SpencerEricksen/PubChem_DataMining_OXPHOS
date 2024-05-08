
import pandas as pd

# load the feat divergence data
df_fps = pd.read_csv('./feature_divergence_plots/fps_bit_density_differences.csv')
#df_frag = pd.read_csv('diff_data_desc_frags_mean_st.csv')
#df_prop = pd.read_csv('diff_data_desc_props_mean_st.csv')
df_fps.rename( columns={'bit_index':'feature'}, inplace=True )

# load the MDI-based feat importance data
df_fps_mdi = pd.read_csv('./MDI/fps/RF_nocalibr_fps_MDI.csv')
df_fps_mdi.rename( columns={'stdev':'mdi_stdev'}, inplace=True )
df_fps_mdi.rename( columns={'mean_decrease_impurity':'mdi_mean'}, inplace=True )

# load the perm-based feat importance data
df_fps_perm = pd.read_csv('./permutation/fps/RF_nocalibr_fps_perm_avgprec_5x.csv')
df_fps_perm['feature'] = df_fps_perm['feature'].astype('int')
df_fps_perm.rename( columns={'perm_feat_importance':'perm_mean'}, inplace=True )
df_fps_perm.rename( columns={'stdev':'perm_stdev'}, inplace=True )

# now merge the div, MDI, and permutation-based importances
df_fps = df_fps.merge( df_fps_mdi, on='feature', how='left')
df_fps = df_fps.merge( df_fps_perm, on='feature', how='left')

# add percentile ranks for each feature type
df = df_fps.sort_values( by='perm_mean', ascending=False )
df['perm_rank'] = df['perm_mean'].rank(pct=True)
df['mdi_rank'] = df['mdi_mean'].rank(pct=True)
# sort divergence values by absolute value
df['diff_rank'] = df['diff'].apply( lambda x: abs(x) ).rank(pct=True)

# dump to CSV
df.to_csv('feature_importance_fps_div_mdi_perm.csv', index=False )

