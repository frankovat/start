import numpy as np
import pandas as pd

TbedDec=pd.read_csv("/Users/terezafrankova/Desktop/Cluster/Cafe/maDec.tsv",sep="\t") #import file with gene and function (olfactory receptor)
TbedDec=TbedDec.drop('gene_function',axis=1) #drop the second column, 'gene_function' is the name of the second column, use your name!
TbedDec=TbedDec.values
TbedDec=TbedDec.flatten()

Orthogroups=pd.read_csv("/Users/terezafrankova/Downloads/Orthogroups.csv",sep='\t') #upload the peptidic file of your species
#Tbed=Orthogroups[['Orthogroups','Taur_pep']].copy() - you will not need this
#Tbed=Tbed.dropna()
#Tbed=Tbed.set_index(['Orthogroups']).apply(lambda x: x.str.split(', ').explode()).reset_index()
#HERE YOU WILL NEED TO CREATE TWO COLUMNS, ONE FROM THE GENES AND SECOND FROM THE SEQUENCES!!! example:
#     q54224.2              AOFHEBSMDLFJSNSBCJCDKJSBDD
#     q3736.13              URHSXCKFPSJFHGGCZSJAKLWOHF
#then you will name the column one as gene and column two as sequence

tbedMinus=Tbed[Tbed['Orthogroups'].isin(TbedDec)] #here you will ask if the column gene is in the created TbedDec
tbedMinus=tbedMinus.drop('Orthogroups',axis=1)
tbedMinus=tbedMinus.values
tbedMinus=tbedMinus.flatten()
#tbedMinus=[s.replace('ta','g') for s in tbedMinus] #changes ta1.12 to g1.12, you don't need to use it!

### AND THIS SHOULD BE ENOUGHT TO GET THE SEQUENCES FROM IT!!!
==========================================================================
tbedGF=pd.read_csv("/Users/terezafrankova/Downloads/TaurGF.csv", sep='\t')
tbedGF=tbedGF.dropna()
tbedGF.drop(tbedGF.loc[tbedGF['Function']=='-'].index, inplace=True)


tbedGF=tbedGF[tbedGF['Gene'].isin(tbedMinus)]

aggregation_functions={'Function':'sum'}
tbedGF=tbedGF.groupby(tbedGF['Gene']).aggregate(aggregation_functions)
print(macuGF)

tbedGF.to_csv('/Users/terezafrankova/Downloads/taurfinal.csv',sep='\t',header=False,index=False)
