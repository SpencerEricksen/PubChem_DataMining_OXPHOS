
import pandas as pd
import numpy as np

# load the original dataframe of all cpd records in OXPHOS-related assays
df1 = pd.read_pickle('oxphos_merged_aids_cids_clnsmi_assaydesc_ETCassay.pkl')
df1['pmid'] = df1['pmid'].replace( "None", 0 )
df1['pmid'] = df1['pmid'].astype( 'int')

# get list of PMIDs that Liping confirmed to be associated with complexes I-IV
df2 = pd.read_csv( './source/data/ETC_inhibitor_PMIDs.csv')
df2['PMID'] = df2['PMID'].astype('int')


pmid_list = df2['PMID'].tolist()

df1['ETC_linked_PMID'] = 0
df1.loc[ df1['pmid'].isin( pmid_list ), 'ETC_linked_PMID'] = 1

# dump with new column indicating ETC_linked_PMID
df1.to_pickle('oxphos_merged_aids_cids_clnsmi_assaydesc_ETCassay_pmid.pkl')


# extras...

# cpd records in original dataframe
# len(df1) = 361124 (cpd records)

# number of cpd records with ETC_linked_PMIDs
#df1['ETC_linked_PMID'].sum()
# 2263 

# number of cpd records with that are active and with ETC_linked_PMIDs
#df1.loc[ (df1['ETC_linked_PMID'] == True) & (df1['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active' ) ]['PUBCHEM_CID']
# 894  

# list of unique compounds (CIDs) that are active and with ETC_linked_PMIDs, and their number of occurrences
#df1.loc[ (df1['ETC_linked_PMID'] == True) & (df1['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active' ) ]['PUBCHEM_CID'].value_counts()
# 316 (number of unique cpds--PUBCHEM_CIDs)



