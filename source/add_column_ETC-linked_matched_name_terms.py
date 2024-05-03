
import pandas as pd
import numpy as np

#### get counts of unique cpds, aids, pmids, for all cpds in original unfiltered data set
df0 = pd.read_csv('all_oxphos_aids_cids_assaydesc.csv', sep="|", low_memory=False)
#df0 = pd.read_csv('ultimate_oxphos_pcassay_dataset_npscores.csv', sep="|", low_memory=False)
orig_aids = set(df0['AID'])
orig_cids = set(df0['PUBCHEM_CID'])
orig_pmids = set(df0['pmid'])
orig_names = set(df0['name'])

#### get counts of the unique active cpds in original unfiltered data set
df0_act = df0.loc[ df0['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active' ]
orig_cids_act  = set(df0_act['PUBCHEM_CID'])
orig_aids_act  = set(df0_act['AID'])
orig_pmids_act = set(df0_act['pmid'])
orig_names_act = set(df0_act['name'])

print('')
print('original data set...')
print( "{:<40}{:<}".format( 'number of unique PUBCHEM_CIDs: ', len(orig_cids) ) )
print( "{:<40}{:<}".format( 'number of unique active PUBCHEM_CIDs: ', len(orig_cids_act) ) )
print('')
print( "{:<40}{:<}".format( 'number of unique PMIDs: ', len(orig_pmids) ) )
print( "{:<40}{:<}".format( 'number of unique active PMIDs: ', len(orig_pmids_act) ) )
print('')
print( "{:<40}{:<}".format( 'number of unique AIDs: ', len(orig_aids) ) )
print( "{:<40}{:<}".format( 'number of unique active AIDs: ', len(orig_aids_act) ) )
print('')
print(f"{'number of assay names: ':<40}{len(orig_names):<}") 
print(f"{'number of active assay names: ':<40}{len(orig_names_act):<}")
print('')

# The original pcassay query (PubChem) returns 4645 assays (match Assay Description):

# ("electron transport chain"[Assay Description] OR 
#  "mitochondrial complex"[Assay Description] OR 
#  "mitochondrial respiratory chain"[Assay Description] OR 
#  "mitochondrial membrane potential"[Assay Description])


#df1 = df0.loc[ df0['PUBCHEM_ACTIVITY_OUTCOME']=='Active' ]

term_list = ['mitochondrial', 
             'mitochondria', 
             'ROS', 
             'electron transport',
             'respiration',
             'respiratory chain',
             'ETC', 
             'OXPHOS', 
             'NADH oxidase', 
             'ubiquinone', 
             'oxidoreductase', 
             'respiration', 
             'ataxia', 
             'hypoxia', 
             'NAD1', 
             'bc1 complex', 
             'b-c1 complex',
             'bc-1 complex',
             'ytochrome c',
             'ytochrome bc1',
             'ytochrome b-c1',
             'ytochrome bc-1',
             'ATP synthase', 
             'oxygen consumption', 
             'superoxide', 
             'radical', 
             'omplex I', 
             'omplex II', 
             'omplex III', 
             'UQCR', 
             'omplex IV', 
             'omplex V', 
             'respiratory chain', 
             'NADH-CoQ reductase', 
             'NADH dehydrogenase',
             'succinate dehydrogenase',
             'succinate-CoQ reductase', 
             'CoQH2', 
             'CoQ', 
             'AA3']

#removed ['reduction', 'sulforhodamine', 'SRB', 'MTT', 'glutathione']

print('adding assay names containing good query terms...')
name_list = []
for i in term_list:
    temp = df0[ df0['name'].str.contains( i )]['name'].unique().tolist()
    print("{:<40}{:<50}".format( 'term: '+i, 'unique name matches: '+str(len(temp)) ) )
    name_list.extend( temp )
print('')

new_names = list(set(name_list))

# remove some mismatch strings
print('removing assay names with "bad" terms...')
bad_terms = ['hloroplast','thylakoid']
for bad_term in bad_terms:
    temp = [ name for name in new_names if not (bad_term in name) ]
    n_matches = len(new_names) - len(temp)
    print("{:<40}{:<}".format( 'term: '+bad_term, 'unique name matches: '+str(n_matches) ) )
    new_names = temp
print('')

#### get counts of unique active cpds, aids, pmids with new query term-matched assay names
df4 = df0[ df0['name'].isin( new_names) ]
new_cids = set(df4['PUBCHEM_CID'])
new_pmids = set(df4['pmid'])
new_aids = set(df4['AID'])
new_names = set(new_names)
print('matches with new "name" list search terms')
print(f"{'number of unique PUBCHEM_CIDs: ':<40}{len(new_cids):<}")
print(f"{'number of unique PMIDs: ':<40}{len(new_pmids):<}")
print(f"{'number of unique AIDs: ':<40}{len(new_aids):<}")
print(f"{'number of assay names: ':<40}{len(new_names):<}")
#print( "{:<40}{:<}".format( 'number of assay names: ', len(new_names) ) )

# adding column to indicate ETC-linked AIDs
print('')
print('adding column to indicate ETC-linked AIDs that match new "name" list....')
df0['ETC_linked_AID'] = False
df0['ETC_linked_AID'].loc[ df0['AID'].isin( new_aids ) ] = True
print('writing new dataframe to CSV...')
df0.to_csv('all_oxphos_aids_cids_assaydesc_ETC.csv', sep="|", index=False)

# df0[ (df0['pmid'].isin( old_pmids - new_pmids )) & (df0['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active') ]['name'].unique() 
print('')
print('list unique cpds that are active in at least 1 ETC-linked AID and number of ETC-linked assays it hits on')
print( df0.loc[ (df0['ETC_linked_AID'] == True) & (df0['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active'), 'PUBCHEM_CID'].value_counts() )

# number of ETC-linked assays each molecule hits on
num_assay_hits = df0.loc[ (df0['ETC_linked_AID'] == True) & (df0['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active'), 'PUBCHEM_CID'].value_counts()


len( df0[ df0['PUBCHEM_CID'].isin( list( num_assay_hits[ num_assay_hits >= 2 ].index) ) ]['smiles'].unique() )

# etc_aid_hits  #_cpds
# 1             1959
# 2             1142
# 3             292
# 4             138
# 5             52
# 6             18
# 7             13
# 8             9
# 9             5
# 10            3
