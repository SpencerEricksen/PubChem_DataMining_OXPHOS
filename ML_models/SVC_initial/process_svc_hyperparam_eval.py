


import pandas as pd
import numpy as np
import sys


incsv = sys.argv[1]  # svc_hyperparams_search.csv
outcsv = sys.argv[2] # svc_hyperparams_search_sorted.csv


# vim (remove header lines after line 1
# :2,$g/c_kernel/d
df = pd.read_csv( incsv )

# make fields for hyperparams
df['kernel'] = df['c_kernel_gamma_feats'].apply( lambda x: x.split('_')[1] )

# specify feature types used 'fps', 'desc', or 'both'
df['feats'] = df['c_kernel_gamma_feats'].apply( lambda x: x.split('_')[-1] )

# all models have C
df['C'] = df['c_kernel_gamma_feats'].apply( lambda x: x.split('_')[0] )

# RBF, poly, and sigmoid kernels use 'gamma'
df['gamma'] = df['c_kernel_gamma_feats'].apply( lambda x: x.split('_')[2] )
df = df.replace('nogamma',np.nan)

# poly kernel has a 'degree'
df['degree'] = np.nan
df['degree'] = df.loc[ df['c_kernel_gamma_feats'].str.contains('poly'), 'c_kernel_gamma_feats' ].apply( lambda x: x.split('_')[-2] )

# both 'poly' and 'sigmoid' kernels have 'coef0'
df['coef0'] = np.nan
df.loc[ df['c_kernel_gamma_feats'].str.contains('poly'), 'coef0'] = df['c_kernel_gamma_feats'].apply( lambda x: x.split('_')[-3] )
df.loc[ df['c_kernel_gamma_feats'].str.contains('sigmoid'), 'coef0'] = df['c_kernel_gamma_feats'].apply( lambda x: x.split('_')[-2] )

# add avg_metric (mean of avg_prec and auc_roc)
df['avg_metric'] = df[['avg_prec','roc_auc']].mean(axis=1)

df2 = df.sort_values(by=['kernel','feats','avg_metric'], ascending=[True, True, False] )
df2.to_csv( outcsv, index=False)

