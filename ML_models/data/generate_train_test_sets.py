

import pandas as pd
from sklearn.model_selection import train_test_split

# load the featurized compound data from OXPHOS assays
df_fps = pd.read_pickle('../../ML_training_data/data_set_fps.pkl')
df_dsc = pd.read_pickle('../../ML_training_data/data_set_desc.pkl')

# merge fps and desc features back into single dataframe before split
df = df_fps.merge( df_dsc.drop( columns='label' ), on='PUBCHEM_CID', how='inner')

# split 85% training, 15% for test
df_train, df_test, y_train, y_test = train_test_split( df, df.label, test_size=0.15, random_state=42 )

# cal_ratio = calibrate_size / ( 1.0 - test_size ) = 17.65%
calibrate_ratio = 0.15 / (1.0 - 0.15)

# split using calibrate ratio (so test and calibrate sizes are the same)
df_train, df_calibrate, y_train, y_calibrate = train_test_split( df_train, df_train.label, test_size=calibrate_ratio, random_state=42 )

# dump splits
df_train.to_pickle('training_data.pkl')
df_test.to_pickle('test_data.pkl')
df_calibrate.to_pickle('calibrate_data.pkl')


