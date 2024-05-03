import pandas as pd

df = pd.read_pickle('unique_pubchem_cids_desc.pkl')

# set all cpds to fail (0)
df['lipinski'] = 0

# update cpds to that pass Lipinski as passing (1)
df.loc[ (df['desc_MolWt'] > 200.0 ) &
        (df['desc_MolLogP'] < 5.8 ) &
        (df['desc_TPSA'] < 150 ) &
        (df['desc_HeavyAtomCount'] > 20) &
        (df['desc_NOCount'] > 0) &
        (df['desc_RingCount'] > 0), 'lipinski' ] = 1

# actual Lipinski
# MW <= 500 Da
# MlogP <= 4.15
# N/O <= 10
# NH/OH <= 5

# Ghose
# 160 <= MW <= 480 Da
# -0.4 <= WlogP <= 5.6
# 40 <= MR <= 130
# 20 <= Heavy Atoms <= 70  ****

# Veber
# RotBonds <= 10
# TPSA <= 140

# Egan
# WlogP <= 5.88         ****
# TPSA <= 131.6

# Muegge
# 200 <= MW <= 600 Da   ****
# -2 <= XlogP <= 5
# TPSA <= 150           ****
# # rings <= 7
# # C atoms > 4
# # heteroatoms > 1     ****
# # RotBnds <= 15
# # HBA <= 10
# # HBD <- 5

# dump to pickle
df[['PUBCHEM_CID','lipinski']].to_pickle('unique_pubchem_cids_lipinski.pkl')

