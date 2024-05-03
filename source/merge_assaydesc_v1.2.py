import pandas as pd
import numpy as np

df1 = pd.read_pickle('oxphos_merged_aids_cids_clnsmi.pkl')
df2 = pd.read_csv('./assay_descriptions/all_assays_desc.csv', sep="|")
df3 = df1.merge( df2, how='left', on='AID')
df3.to_pickle('oxphos_merged_aids_cids_clnsmi_assaydesc.pkl')

