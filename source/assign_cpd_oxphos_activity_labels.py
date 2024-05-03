
import pandas as pd

df = pd.read_pickle('unique_pubchem_cids_complete_assay_results.pkl')

df['label'] = 'unlabeled'

# assign simple labels

df['etc_linked_bool_map'] = df['ETC_linked_AID'].apply( lambda x: [ i==1 for i in x ] )
# keep only ETC-linked AID pc_act_outcomes
#df['pc_act_outcomes_etc'] = [b for a, b in zip( df['etc_linked_bool_map'], df['pc_act_outcomes']) if a]

# label as 'active' cpds with at least one active record and also number of active records outnumbers inactive records 
df.loc[ (df['num_active'] > 0) & (df['num_active'] >= df['num_inactive']), 'label' ] = 1

# label as 'inactive' cpds with at least one inactive record and no 'active', 'inconclusive', or 'unspecified' records
df.loc[ (df['num_inactive'] > 0) & (df['num_active'] == 0) & (df['num_inconclusive'] == 0) & (df['num_unspecified'] == 0), 'label' ] = 0

# dump file
df.to_pickle('unique_pubchem_cids_complete_assay_results_w_labels.pkl')


