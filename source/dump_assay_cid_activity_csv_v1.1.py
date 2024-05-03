
import pandas as pd
import numpy as np
import sys

# pcba-aid1055297.csv
input_aid = sys.argv[1]

aid_number = str( input_aid.split('-')[-1].split('.')[0][3:] )

df1 = pd.read_csv( input_aid, low_memory=False )
df2 = df1[ pd.to_numeric( df1['PUBCHEM_CID'], errors='coerce' ).notnull() ]
df3 = df2[['PUBCHEM_CID', 'PUBCHEM_ACTIVITY_OUTCOME', 'PUBCHEM_ACTIVITY_SCORE' ] ]
df3 = df3.astype( {'PUBCHEM_CID':'int64','PUBCHEM_ACTIVITY_SCORE':'float64'} )
df4 = df3.sort_values( by='PUBCHEM_ACTIVITY_SCORE', ascending=False ).drop_duplicates( subset='PUBCHEM_CID', keep='first')
df4.set_index( 'PUBCHEM_CID', inplace=True)
df4.to_csv('./CID_lists/pcba-aid'+aid_number+'_activities.csv', index_label='PUBCHEM_CID')

