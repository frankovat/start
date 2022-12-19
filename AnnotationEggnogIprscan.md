Eggnog mapper
https://wiki.metacentrum.cz/wiki/Conda_-_modules
Interproscan
https://wiki.metacentrum.cz/wiki/Interproscan_


Interproscan
```
/mnt/gpfsA/home/frankova/software/my_interproscan/interproscan-5.59-91.0/interproscan.sh -i $HOMEDIR/Ofo_0.1_longestiso.fa 
\--seqtype p -goterms -iprlookup -b Ofo_longestiso.fa.iprscan -pa -cpu 16
```

Eggnog

```
emapper.py -m diamond --sensmode more-sensitive --no_annot -i Mnat_0.1_longestiso.fa -o Mnat_diamond
emapper.py -m no_search --annotate_hits_table Mnat_diamond1.emapper.seed_orthologs -o Mnat_diamond --dbmem
```
