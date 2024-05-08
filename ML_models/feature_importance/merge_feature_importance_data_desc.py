
import pandas as pd

# load the feat divergence data
df_frag = pd.read_csv('./feature_divergence_plots/diff_data_desc_frags_mean_st.csv')
df_prop = pd.read_csv('./feature_divergence_plots/diff_data_desc_props_mean_st.csv')
df_prop.rename( columns={'prop':'feature'}, inplace=True )
df_frag.rename( columns={'frag':'feature'}, inplace=True )
df_prop = df_prop.drop( columns='colors' )
df_frag = df_frag.drop( columns='colors' )
df_div = pd.concat( [df_prop, df_frag], axis=0 )

# load the MDI-based feat importance data
df_desc_mdi = pd.read_csv('./MDI/desc/RF_nocalibr_desc_MDI.csv')
df_desc_mdi.rename( columns={'stdev':'mdi_stdev'}, inplace=True )
df_desc_mdi.rename( columns={'mean_decrease_impurity':'mdi_mean'}, inplace=True )

# load the perm-based feat importance data
df_desc_perm = pd.read_csv('./permutation/desc/RF_nocalibr_desc_perm_feat_importances_test_avgprec.csv')
df_desc_perm.rename( columns={'perm_feat_importance':'perm_mean'}, inplace=True )
df_desc_perm.rename( columns={'stdev':'perm_stdev'}, inplace=True )

# now merge the div, MDI, and permutation-based importances
df_all = df_div.merge( df_desc_mdi, on='feature', how='left')
df_all = df_all.merge( df_desc_perm, on='feature', how='left')

# add percentile ranks for each feature type
df = df_all.sort_values( by='perm_mean', ascending=False )
df['perm_rank'] = df['perm_mean'].rank(pct=True)
df['mdi_rank'] = df['mdi_mean'].rank(pct=True)
# sort divergence values by absolute value
df['diff_rank'] = df['diff'].apply( lambda x: abs(x) ).rank(pct=True)

# dump to CSV
df.to_csv('feature_importance_desc_div_mdi_perm.csv', index=False )
