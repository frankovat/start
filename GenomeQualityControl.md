After we finish the assembly, we have to asses the quality of the assembly. For this, we can use 3 programs. BUSCO, fCat and taXaminer.

# BUSCO
You can use the already preinstalled python environment in Metacentrum (https://wiki.metacentrum.cz/wiki/BUSCO). It is the easiest way to asses the quality of the genome assembly: BUSCO is fast and reliable.
To understand what BUSCO is telling us:

![This is an image](https://training.galaxyproject.org/archive/2019-12-01/topics/genome-annotation/images/busco_genome_summary.png) \n
So you want to have Complete and Single Copy BUSCOs as close to the Complete BUSCOs, and very little of duplicated BUSCOs, fragmented and missing BUSCOs. The closer you get to the Complete BUSCOs the better the assembly is. 

```
 #!/bin/bash
 #PBS -l select=1:ncpus=12:ompthreads=12:mem=15gb:scratch_local=100gb
 #PBS -l walltime=24:00:00 ##24 is a bit too much tbh, you can have just around 3-5 hours. It was enough for me
 #PBS -m abe

 trap 'clean_scratch' TERM EXIT

 DATADIR="/your/assembled/genome/"

 cp $DATADIR/yourAssembledGenome.fasta $SCRATCHDIR || exit 1
 cd $SCRATCHDIR || exit 2

 module add conda-modules-py37 ##this adds the module of anaconda, where the busco environment is

 conda activate busco ##we activate the environment with this command

 mkdir $SCRATCHDIR/augustus
 cp -r /software/augustus/3.4.0 $SCRATCHDIR/augustus/3.4.0
 export AUGUSTUS_CONFIG_PATH=$SCRATCHDIR/augustus/3.4.0/config ##this creates a writable directory of augustus configuration directory, otherwise it doesn´t work

 busco -i yourAssembledGenome.fasta -o outputDirectory -l configurationData -m geno -c 12 --long --augustus_species augustusSpecies
 ## -l is a configurationData you will choose from BUSCO, you can find the list bellow this box
 ## --augustus_species is another configuration file you can choose, you can find the list bellow this box
 conda deactivate ##we deactivate the environment

 cp -r outputDirectory $DATADIR || export CLEAN_SCRATCH=false
```
