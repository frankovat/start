# Trinity GG = genome guided transcriptome assembly
Now we run trinity again, but this time the genome-guided version. However, this is a bit more complicated because we have to prepare our genome first with hisat2.
On metacentrum, you can easily find only [hisat](https://wiki.metacentrum.cz/wiki/Hisat) and not hisat2. It is said, that `HISAT2 is distributed under the GPLv3 license` which for me right now doesn´t make any sense.
But when you have a look at the available modules, hisat2 is there, so it shouldn´t be a problem to run it. *Don´t forget to add the modules into the script!!!*

```
#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=100000

module add hisat2-2.0.5 #I added it here to remind you to not forget
#here create a dirrectory in the scratchdir for the outputs
hisat2-build /path/to/your/assembled/genome.fasta /here/is/the/output/directory

cd /here/is/the/output/directory

hisat2 --max-intronlen 100000 -p 10 -x /theDirectory/you/copied/from/hisat2-build/subdir
_build --phred33 -1 /path/to/rna+.fastq.gz -2 /path/to/rna-.fastq.gz | samtools view -bS - | samtools sort -o /path/to/output/dir/name_sorted.bam

Trinity --genome_guided_bam /path/to/output/dir/name_sorted.bam --max_memory 100G --CPU 10 --SS_lib_type RF --genome_guided_max_intron 100000 --output /path/to/output/dir/name_trinity_GG

cp /trinity_GG/name/name_trinity_GG_jaccard/Trinity-GG.fasta $HOMEDIR/trinity_GG/name/name_Trinity-GG.fasta
cp /trinity_GG/LLEU/LLEU_sorted.bam $HOMEDIR/trinity_GG/name/name_sorted.bam

hisat2 --max-intronlen 100000 -p 10 -x /theDirectory/you/copied/from/hisat2-build/subdir
_build --phred33 -1 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/S81A_R1.fastq.gz -2 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/S81A_R2.fastq.gz | samtools view -bS - | samtools sort -o $SCRATCHDIR

/Genomics/kocherlab/berubin/local/src/Trinity/Trinity --genome_guided_bam /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted_callum.bam --max_memory 100G --CPU 10 --genome_guided_max_intron 100000 --output /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG_callum

cp /trinity_GG/name/name_trinity_GG_callum/Trinity-GG.fasta $HOMEDIR/trinity_GG/name/name_Trinity-GG_callum.fasta
cp /trinity_GG/name/name_sorted_callum.bam $HOMEDIR/trinity_GG/name/name_sorted_callum.bam
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

|Trinity command specification|explanation|
|---|---|
|--genome_guided_bam|here we used the output file that is name_sorted.bam|
|--max_memory|the maximum memory allocated for this task should be the same as in PBS|
|--CPU|again, number of CPUs, same as in PBS|
|--SS_lib_type|:Strand-specific RNA-Seq read orientation. if paired: RF or FR, if single: F or R.   (dUTP method = RF) See web documentation.(copied from trinity manual)|
|--genome_guided_max_intron|again, this was tested on halictids, leave it like that|
|--output||
