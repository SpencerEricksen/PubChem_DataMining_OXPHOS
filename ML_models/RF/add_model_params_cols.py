

import pandas as pd

# read in hyperparam sweep data
df = pd.read_csv('df_opt1-4_desc.log')

# add individual cols for each model parameter
df['model_ntrees'] = df['ntrees_nsamples_ndepth_nfeats'].apply( lambda x: int( x.split('_')[0])  )
df['model_nsamples'] = df['ntrees_nsamples_ndepth_nfeats'].apply( lambda x: float( x.split('_')[1])  )
df['model_ndepth'] = df['ntrees_nsamples_ndepth_nfeats'].apply( lambda x: str( x.split('_')[2])  )
df['model_nfeats'] = df['ntrees_nsamples_ndepth_nfeats'].apply( lambda x: float( x.split('_')[3])  )

# reorder columns
df2 = df[ df.columns[-4:].tolist() + df.columns[0:-4].tolist() ]

