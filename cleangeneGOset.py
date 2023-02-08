import pandas as pd
import numpy as np

# DATASET PREPARATION
## GO annotation from iprscan and eggnog
setOfor = pd.read_csv("/Users/terezafrankova/Documents/GeneSetStart/OforconcatGO.tsv",sep='\t')
print(setOfor)

# Cleaning dataset from blank values
setOfor2 = setOfor[setOfor[setOfor.columns[1]] != '-']
print(setOfor2) # cleaned blank spaces
setOfor = setOfor2[setOfor2[setOfor2.columns[1]].notna()]
print(setOfor) # cleaned NaN rows

# Separate GO terms into separate rows
setOfor = setOfor.set_index(['query']).apply(lambda x: x.str.split(',').explode()).reset_index()
setOfor = setOfor.set_index(['query']).apply(lambda x: x.str.split('|').explode()).reset_index()

# Drop duplicate lines
setOfor = setOfor.drop_duplicates()

#Save file and move to R
setOfor3.to_csv('/Users/terezafrankova/Documents/GeneSetStart/OforGOforR.tsv',sep='\t', header=True, index=False)
