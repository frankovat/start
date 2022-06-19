# Trinity GG = genome guided transcriptome assembly
Now we run trinity again, but this time the genome-guided version. However, this is a bit more complicated because we have to prepare our genome first with hisat2.
On metacentrum, you can easily find only [hisat](https://wiki.metacentrum.cz/wiki/Hisat) and not hisat2. It is said, that `HISAT2 is distributed under the GPLv3 license` which for me right now doesn´t make any sense.
But when you have a look at the available modules, hisat2 is there, so it shouldn´t be a problem to run it. *Don´t forget to add the modules into the script!!!*

```
#!/bin/bash
#PBS -l select=1:ncpus=10:mem=250gb:scratch_local=200gb
#PBS -l walltime=72:00:00
#PBS -m abe

trap 'clean_scratch' TERM EXIT

OUTDIR="/storage/plzen1/home/frankovat/TRINITY/trinity_GG/"

cp /auto/plzen1/home/frankovat/GenomesTranscriptomes/Transcriptomes/CchalT/CchalT_1.fq.gz $SCRATCHDIR || exit 1
cp /auto/plzen1/home/frankovat/GenomesTranscriptomes/Transcriptomes/CchalT/CchalT_2.fq.gz $SCRATCHDIR || exit 2
cp /auto/plzen1/home/frankovat/Postprocessing/Gapcloser/Cchal/Cchal_gapcloser.fasta $SCRATCHDIR || exit 3

cd $SCRATCHDIR

mkdir -p trinity_GG/Cchal/Cchal_build

module load samtools-1.9
module load hisat2-2.2.1
#module load trinity-2.9.1

hisat2-build $SCRATCHDIR/Cchal_gapcloser.fasta $SCRATCHDIR/trinity_GG/Cchal/Cchal_build

hisat2 --max-intronlen 100000 -p 10 -x $SCRATCHDIR/trinity_GG/Cchal/Cchal_build --phred33 -1 CchalT_1.fq.gz -2 CchalT_2.fq.gz | samtools view -bS - | samtools sort -o $SCRATCHDIR/trinity_GG/Cchal/Cchal_sorted.bam

Trinity --no_version_check --genome_guided_bam $SCRATCHDIR/Cchal_sorted.bam --max_memory 200G --CPU 10 --SS_lib_type RF --genome_guided_max_intron 100000 --output $SCRATCHDIR/trinity_GG/Cchal/Cchal_trinity_GG

cp -r $SCRATCHDIR $OUTDIR || export CLEAN_SCRATCH=false
```
## Hisat2
[Here](http://daehwankimlab.github.io/hisat2/manual/) is a manual to hisat2. 
First, we will run `hisat2-build` command, to create indexes for our genome. Later, we will be using the command `hisat2` which creates the alignment.
|hisat2 command specification|explanation|
|---|---|
|--max-intronlen|this was tested on halictides and I would not change the number, the maximum length of introns|
|-p|parallelization, aka number of CPUs used, same as in the PBS|
|-x|here, you have to guide the program to the dirrectory you copied after running hisat_build(dont forget you have to copy the subdir NAME_build|
|--phred33|the ASCII qualities in fastq file|
|-1|the + string of the RNA data|
|-2|path to - string of the RNA data|
|samtools view and sort| here we make a pipeline with the straight slash ad piping the sorting with samtools - don´t forget the samtools module!|

## Trinity
Previous file mentiones Trinity
|Trinity command specification|explanation|
|---|---|
|--genome_guided_bam|here we used the output file that is name_sorted.bam|
|--max_memory|the maximum memory allocated for this task should be the same as in PBS|
|--CPU|again, number of CPUs, same as in PBS|
|--SS_lib_type|:Strand-specific RNA-Seq read orientation. if paired: RF or FR, if single: F or R.   (dUTP method = RF) See web documentation.(copied from trinity manual)|
|--genome_guided_max_intron|again, this was tested on halictids, leave it like that|
|--output||


```
/Genomics/kocherlab/berubin/local/src/hisat2-2.0.5/hisat2-build /Genomics/kocherlab/berubin/assembly/hic/LLEU/LLEU_genome_v2.0.fasta /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_build

cd /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU

/Genomics/kocherlab/berubin/local/src/hisat2-2.0.5/hisat2 --max-intronlen 100000 -p 10 -x /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_build --phred33 -1 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_head_R1.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_abd_R1.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_ant_R1.fastq.gz -2 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_head_R2.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_abd_R2.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_ant_R2.fastq.gz | samtools view -bS - | samtools sort -o /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted.bam

/Genomics/kocherlab/berubin/local/src/Trinity/Trinity --genome_guided_bam /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted.bam --max_memory 100G --CPU 10 --SS_lib_type RF --genome_guided_max_intron 100000 --output /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG

cp /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG_jaccard/Trinity-GG.fasta /Genomics/kocherlab/berubin/transcriptomes/trinity_GG/LLEU/LLEU_Trinity-GG.fasta
cp /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted.bam /Genomics/kocherlab/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted.bam
rm -r /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG/

/Genomics/kocherlab/berubin/local/src/hisat2-2.0.5/hisat2 --max-intronlen 100000 -p 10 -x /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_build --phred33 -1 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/S81A_R1.fastq.gz -2 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/S81A_R2.fastq.gz | samtools view -bS - | samtools sort -o /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted_callum.bam

/Genomics/kocherlab/berubin/local/src/Trinity/Trinity --genome_guided_bam /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted_callum.bam --max_memory 100G --CPU 10 --genome_guided_max_intron 100000 --output /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG_callum

cp /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG_callum/Trinity-GG.fasta /Genomics/kocherlab/berubin/transcriptomes/trinity_GG/LLEU/LLEU_Trinity-GG_callum.fasta
cp /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted_callum.bam /Genomics/kocherlab/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted_callum.bam
rm -r /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG_callum/

rm -r /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/
```
