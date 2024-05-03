
import pandas as pd
import numpy as np


df = pd.read_pickle('oxphos_merged_aids_cids_clnsmi_assaydesc.pkl')

# alternatively could try converting name, Title, and Abstract to all 
# lower case to simplify terms search

# df[colname] = df[colname].str.lower()

#### get counts of unique cpds, aids, pmids, for all cpds in original unfiltered data set
#### get counts of the unique active cpds in original unfiltered data set

good_terms = ['mitochondrial', 
             'mitochondria', 
             'ROS', 
             'electron transport',
             'respiration',
             'respiratory',
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
             'NADH-CoQ reductase', 
             'NADH dehydrogenase',
             'succinate dehydrogenase',
             'succinate-CoQ reductase', 
             'CoQH2', 
             'CoQ', 
             'AA3']

bad_terms = ['hloroplast','thylakoid']

#removed ['reduction', 'sulforhodamine', 'SRB', 'MTT', 'glutathione']

df.loc[ df['Title'].isnull(), 'Title' ] = ""
df.loc[ df['Abstract'].isnull(), 'Abstract' ] = ""

df['ETC_linked_AID'] = 0
for term in good_terms:
    df.loc[ df['name'].str.contains( term ), 'ETC_linked_AID' ] = 1
    df.loc[ df['Title'].str.contains( term ), 'ETC_linked_AID' ] = 1
    df.loc[ df['Abstract'].str.contains( term ), 'ETC_linked_AID' ] = 1

for term in bad_terms:
    df.loc[ df['name'].str.contains( term ), 'ETC_linked_AID' ] = 0
    df.loc[ df['Title'].str.contains( term ), 'ETC_linked_AID' ] = 0
    df.loc[ df['Abstract'].str.contains( term ), 'ETC_linked_AID' ] = 0

df.to_pickle('oxphos_merged_aids_cids_clnsmi_assaydesc_ETCassay.pkl')

#print('list unique cpds that are active in at least 1 ETC-linked AID and number of ETC-linked assays it hits on')
#print( df0.loc[ (df0['ETC_linked_AID'] == True) & (df0['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active'), 'PUBCHEM_CID'].value_counts() )

# number of ETC-linked assays each molecule hits on
#num_assay_hits = df0.loc[ (df0['ETC_linked_AID'] == True) & (df0['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active'), 'PUBCHEM_CID'].value_counts()



