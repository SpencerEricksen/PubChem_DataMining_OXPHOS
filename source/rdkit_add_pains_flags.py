#!/home/ssericksen/anaconda2/envs/py36_chem/bin/python3.6

# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import sys
from rdkit import Chem

def pains_flags_from_smi( smi ):
    try:
        mol = Chem.MolFromSmiles( smi )
        for k,v in dic.items():
            subs = Chem.MolFromSmarts( k )
            if subs != None:
                if mol.HasSubstructMatch( subs ):
                    mol.SetProp(v,k)
        props = [ prop for prop in mol.GetPropNames() ]
        if len(props) == 0:
            props = False
    except:
        props = False 
        pass
    return props

# load CSV file into dataframe
inpkl = sys.argv[1]
outpkl = sys.argv[2]

#df = pd.read_csv( incsv, sep="|" )
df = pd.read_pickle( inpkl )

# get the pains definitions, load into dict
inf = open("./source/data/pains.txt", "r")
sub_strct = [ line.rstrip().split(" ") for line in inf ]

smarts = [ line[0] for line in sub_strct]
desc = [ line[1] for line in sub_strct]
dic = dict(zip(smarts, desc))

# add pains flags
df['pains'] = df['rdkit_smiles_cln'].apply( pains_flags_from_smi )

df.to_pickle( outpkl )

