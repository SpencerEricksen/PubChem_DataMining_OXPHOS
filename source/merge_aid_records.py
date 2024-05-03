
from glob import glob
import pandas as pd
import numpy as np

aids = glob("AIDs/pcba-aid*.csv")

df_list = []

for f in aids:
    aid = f.split('/')[1].split('.')[0][8:]
    temp_df = pd.read_csv(f)
    temp_df = temp_df[ temp_df.columns[0:6] ]
    temp_df['AID'] = aid
    df_list.append( temp_df )

df = pd.concat( df_list, axis=0 )


df.dropna( subset=['PUBCHEM_CID'], inplace=True )

df2 = df.loc[ df['PUBCHEM_ACTIVITY_SCORE'] != 'PUBCHEM_ACTIVITY_SCORE' ]
df2['PUBCHEM_ACTIVITY_SCORE'] = df2['PUBCHEM_ACTIVITY_SCORE'].astype('float')
df2['PUBCHEM_CID'] = df2['PUBCHEM_CID'].astype('int')
df2['PUBCHEM_SID'] = df2['PUBCHEM_SID'].astype('int')
df2['AID'] = df2['AID'].astype('int')


df2.to_pickle('oxphos_merged_aids_records.pkl')
# df = pd.read_pickle('oxphos_merged_aids_records.pkl')

