#BRAKER
[Braker](https://github.com/Gaius-Augustus/BRAKER#what-is-braker) is a genome annotation pipeline implementing several programs. It combines [GeneMark-ET](https://github.com/gatech-genemark/GeneMark-EP-plus) and [Augustus](https://github.com/Gaius-Augustus/Augustus). BRAKER pipeline implements RNA-seq data to predict gene structure annotation in a novel genome.
By now, you should be familiar with Augustus, because we use it while running BUSCO.
GeneMark uses a [Hidden Markov Model](https://www.sciencedirect.com/topics/medicine-and-dentistry/hidden-markov-model), which is quite an interesting thing, so go ahead and read it if you want. Augustus on the other hand uses a [Generalized Hidden Markov Model] (https://pubmed.ncbi.nlm.nih.gov/8877513/). It will help you to understand how the annotations are carried out, if you are interested :)

```
#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=100000

module add BRAKER1

module add exonerate

module add ncbi-blast

hisat2-build /path/to/the/masked/genome.fasta.masked /annotation/from/braker/name/name_build # don´t forget to copy this directory

cd $SCRATCHDIR/annotation/from/braker/name

hisat2 --max-intronlen 100000 -p 10 -x /annotation/from/braker/name/name_build --phred33 -1 /path/to/transcriptome/+strand -2 /path/to/transcriptome/-strand | samtools view -bS - | samtools sort -o /annotation/from/braker/name/name_sorted.bam

braker.pl --genome=/path/to/the/masked/genome.fasta.masked --bam=/annotation/from/braker/name/name_sorted.bam --species=name --gff3 --cores=10 --AUGUSTUS_CONFIG_PATH=/path/to/augustus/config/ --GENEMARK_PATH=/path/to/gm_et_linux_64/gmes_petap/ --softmasking --overwrite --AUGUSTUS_ab_initio

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
|AUGUSTUS_CONFIG_PATH|you should be familiar with this from BUSCO|
|--GENEMARK_PATH|path to gmes_petap.pl script, I forgot what this one is about and can´t find it online, just many errors with this perl script, so hopefully you will not have that many problems with it|
|--softmasking|we softmasked our genome before, so we are letting braker know that|
|--overwrite|overwrites existing files so we do not take much space|
|--AUGUSTUS_ab_initio|it gives us ab initio predictions of genes in the genome|
