# BRAKER
[Braker](https://github.com/Gaius-Augustus/BRAKER#what-is-braker) is a genome annotation pipeline implementing several programs. It combines [GeneMark-ET](https://github.com/gatech-genemark/GeneMark-EP-plus) and [Augustus](https://github.com/Gaius-Augustus/Augustus). BRAKER pipeline implements RNA-seq data to predict gene structure annotation in a novel genome.
By now, you should be familiar with Augustus, because we use it while running BUSCO.
GeneMark uses a [Hidden Markov Model](https://www.sciencedirect.com/topics/medicine-and-dentistry/hidden-markov-model), which is quite an interesting thing, so go ahead and read it if you want. Augustus on the other hand uses a [Generalized Hidden Markov Model](https://pubmed.ncbi.nlm.nih.gov/8877513/). It will help you to understand how the annotations are carried out, if you are interested :)

```
#!/bin/bash
#PBS -l select=1:ncpus=10:mem=250gb:scratch_local=200gb
#PBS -l walltime=20:00:00
#PBS -m abe

OUTDIR="/storage/plzen1/home/frankovat/BRAKER/Cchal"

test -n "$SCRATCHDIR" || exit 1
cp /auto/plzen1/home/frankovat/RepeatMasker/Cchal/Cchal_rm_masked/Cchal_gapcloser.fasta.masked $SCRATCHDIR || exit 1
cp /storage/plzen1/home/frankovat/BRAKER/Cchal/job_11802907.meta-pbs.metacentrum.cz/Cchal/Cchal_sorted.bam $SCRATCHDIR || exit 2
cd $SCRATCHDIR || exit 3

module load braker2-2.1.6
module load exonerate-2.2.0

export AUGUSTUS_CONFIG_PATH=$SCRATCHDIR/augustus_config
cp -r /software/augustus/3.4.0/config $AUGUSTUS_CONFIG_PATH || exit 3
chmod -R u+rwX $AUGUSTUS_CONFIG_PATH

braker.pl --genome=$SCRATCHDIR/Cchal_gapcloser.fasta.masked --bam=$SCRATCHDIR/Cchal_sorted.bam --species=CchalBraker --gff3 --cores=10 --softmasking --overwrite --AUGUSTUS_ab_initio

cp -r $SCRATCHDIR $OUTDIR || exit 4
clean_scratch

#a nebo můžete vykopírovat přímo tyto dva soubory
cp /annotation/braker/name/braker/name_braker/augustus.hints.gff3 $HOMEDIR

cp /annotation/braker/name/braker/name_braker/augustus.ab_initio.gff3 $HOMEDIR
```
First, we have to index the DNA again with `hisat2-build` as we did before and `hisat` consequentially. If you forgot, what these two do, just go back a couple of files.
Afterwards, we finally implement the [braker.pl](https://github.com/Gaius-Augustus/BRAKER/blob/master/scripts/braker.pl) file, which you can read yourself where everything is explained. I will still make a tab for you here.
|braker.pl command|explanation|
|---|---|
|--genome|of course this is a path to the **masked** genome from RepeatMasker|
|--bam|path to the bam file we created with hisat2|
|--species|here you assign the name of the species, so Ccyp for C. cypriaca and so on|
|--gff3|we want an output in [gff3](https://learn.gencore.bio.nyu.edu/ngs-file-formats/gff3-format/) format|
|--cores|again, you assign the same number of CPUs used in PBS|
|--AUGUSTUS_CONFIG_PATH|you should be familiar with this from BUSCO|
|--GENEMARK_PATH|path to gmes_petap.pl script, I forgot what this one is about and can´t find it online, just many errors with this perl script, so hopefully you will not have that many problems with it|
|--softmasking|we softmasked our genome before, so we are letting braker know that|
|--overwrite|overwrites existing files so we do not take much space|
|--AUGUSTUS_ab_initio|it gives us ab initio predictions of genes in the genome|

**DON´T FORGET TO NAME ALL THE FILES PROPERLY SO YOU DON´T GET LOST! IN THIS PHASE, WE ARE GENERATING A HUGE AMOUNT OF FILES AND IT FEELS LIKE A LABYRINTH WHEN YOU HAVE A MESS IN IT!!!**
