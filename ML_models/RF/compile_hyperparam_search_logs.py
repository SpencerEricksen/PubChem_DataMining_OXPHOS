



import pandas as pd


df1 = pd.read_csv('rf_opt1-5_desc_log.csv')

# new log file to add into the mix
df2 = pd.read_csv('rf_opt6_desc.log')

# fix column headers in new log file
df2.columns = [ i.strip() for i in df2.columns ]

# merge data sets, get average metric, sort descending
df3 = pd.concat([df1, df2], axis=0 )
df3['avg_metrics'] = df3[['bal_acc','avg_prec','roc_auc','f1']].mean(axis=1)
df3.sort_values( by='avg_metrics', ascending=False, inplace=True )
df3.to_csv('rf_opt1-6_desc_log.csv', index=False )


df4 = df3.dropna()

