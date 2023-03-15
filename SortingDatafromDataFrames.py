import numpy as np
import pandas as pd

TbedDec=pd.read_csv("/Users/terezafrankova/Desktop/Cluster/Cafe/maDec.tsv",sep="\t")
TbedDec=TbedDec.drop('First',axis=1)
TbedDec=TbedDec.values
TbedDec=TbedDec.flatten()

Orthogroups=pd.read_csv("/Users/terezafrankova/Downloads/Orthogroups.csv",sep='\t')
Tbed=Orthogroups[['Orthogroups','Taur_pep']].copy()
Tbed=Tbed.dropna()
Tbed=Tbed.set_index(['Orthogroups']).apply(lambda x: x.str.split(', ').explode()).reset_index()

tbedMinus=Tbed[Tbed['Orthogroups'].isin(TbedDec)]
tbedMinus=tbedMinus.drop('Orthogroups',axis=1)
tbedMinus=tbedMinus.values
tbedMinus=tbedMinus.flatten()
tbedMinus=[s.replace('ta','g') for s in tbedMinus]

tbedGF=pd.read_csv("/Users/terezafrankova/Downloads/TaurGF.csv", sep='\t')
tbedGF=tbedGF.dropna()
tbedGF.drop(tbedGF.loc[tbedGF['Function']=='-'].index, inplace=True)


tbedGF=tbedGF[tbedGF['Gene'].isin(tbedMinus)]

aggregation_functions={'Function':'sum'}
tbedGF=tbedGF.groupby(tbedGF['Gene']).aggregate(aggregation_functions)
print(macuGF)

tbedGF.to_csv('/Users/terezafrankova/Downloads/taurfinal.csv',sep='\t',header=False,index=False)
