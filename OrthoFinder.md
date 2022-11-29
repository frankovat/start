MANUAL: https://github.com/davidemms/OrthoFinder/blob/master/OrthoFinder-manual.pdf

OFFICIAL PAPER: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y

If you do not understand something (like a terminology or the pipeline itself), try finding the answer! Alternatively, you can post your questions here: https://github.com/frankovat/start/discussions/2


Orthofinder can be found on metacentrum https://wiki.metacentrum.cz/wiki/Orthofinder
For me it ran about 2 hours I think? With the data here. So run it now with about 6 hours and prolong it if needed.
This is the code:
'''
#!/bin/bash
#PBS -q SMP
#PBS -l ncpus=8
#PBS -l mem=100gb

HOMEDIR="/mnt/gpfsA/home/frankova/A-Compare/Orthofinder/data"

source /mnt/gpfsA/home/frankova/miniconda3/bin/activate orthofinder

/mnt/gpfsA/home/frankova/software/OrthoFinder_source/orthofinder.py -f $HOMEDIR -t 8
'''
