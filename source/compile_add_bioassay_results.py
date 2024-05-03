
import pandas as pd

# load master dataframe
df1 = pd.read_pickle('oxphos_merged_aids_cids_clnsmi_assaydesc_ETCassay_pmid.pkl')

# pivot master to build unique compound dataframe -- each row a unique CID with assays and results listed
pt = pd.pivot_table( df1, values=['AID','PUBCHEM_ACTIVITY_OUTCOME','PUBCHEM_ACTIVITY_SCORE','ETC_linked_AID','ETC_linked_PMID'], index='PUBCHEM_CID', aggfunc={ 'AID':list, 'PUBCHEM_ACTIVITY_OUTCOME':list, 'PUBCHEM_ACTIVITY_SCORE':list, 'ETC_linked_AID':list, 'ETC_linked_PMID':list } )

# add columns with counts of outcomes
pt['num_aids'] = pt['AID'].apply( lambda x: len(x) )
pt['num_active']       = pt['PUBCHEM_ACTIVITY_OUTCOME'].apply( lambda x: sum( i == "Active" for i in x ) )
pt['num_inactive']     = pt['PUBCHEM_ACTIVITY_OUTCOME'].apply( lambda x: sum( i == "Inactive" for i in x ) )
pt['num_inconclusive'] = pt['PUBCHEM_ACTIVITY_OUTCOME'].apply( lambda x: sum( i == 'Inconclusive' for i in x ) )
pt['num_unspecified']  = pt['PUBCHEM_ACTIVITY_OUTCOME'].apply( lambda x: sum( i == 'Unspecified' for i in x ) )
pt['num_etc_linked_aids'] = pt['ETC_linked_AID'].apply( lambda x: sum(x) )
pt['num_etc_linked_pmids'] = pt['ETC_linked_PMID'].apply( lambda x: sum(x) )

# rename some cols
pt.rename( columns={'AID':'aids','PUBCHEM_ACTIVITY_OUTCOME':'pc_act_outcomes','PUBCHEM_ACTIVITY_SCORE':'pc_act_scores'}, inplace=True )
pt = pt.reset_index()

# get unique set of CIDs (w rdkit_smiles_cln)
df2 = pd.read_pickle('unique_pubchem_cids.pkl')

# merge the complete assay info for each unique CID
df3 = df2.merge( pt, on='PUBCHEM_CID', how='left')
# dump to pickle
df3.to_pickle('unique_pubchem_cids_complete_assay_results.pkl')


