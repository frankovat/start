# Trinity - transcriptome preparation
Trinity, a RNA-seq data assembler, is available at [Metacentrum](https://wiki.metacentrum.cz/wiki/Trinity).
We will need the assembled transcriptomes for all gene prediction programs, that is why we are preparing them now.
The [trinity manual](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running%20Trinity) gives you enough answers tu successfully run the program. 
However, I will be sharing with you a script from Priceton, which is easy to adjust to our computational resources.

```
#!/bin/bash
#PBS -l select=1:ncpus=10:mem=200gb:scratch_local=150gb
#PBS -l walltime=74:00:00
#PBS -m abe 

trap 'clean_scratch' TERM EXIT

DATADIR="/storage/plzen1/home/sharu/sharu/rna/AclaT"
OUTDIR="/storage/plzen1/home/sharu/sharu/rna/AclaT/Acla_trinity"

cp $DATADIR/AclaT_1.fq.gz $SCRATCHDIR || exit 1
cp $DATADIR/AclaT_2.fq.gz $SCRATCHDIR || exit 2

cd $SCRATCHDIR || exit 3

module load trinity-2.9.1

Trinity --no_version_check --seqType fq --left /storage/plzen1/home/sharu/sharu/rna/AclaT/AclaT_1.fq.gz --right /storage/plzen1/home/sharu/sharu/rna/AclaT/AclaT_2.fq.gz --max_memory 200G --CPU 10 --trimmomatic --jaccard_clip --output /storage/plzen1/home/sharu/sharu/rna/AclaT/Acla_trinity

export CLEAN_SCRATCH=false

```
|Command | Explanation|
|---|---|
|--seqType fq|the type of the sequence is fastq|
|--left|here we pass the + string (or the first sequence)|
|--right| here we pass the - string (or the second sequence)|
|--max-memory|the maximum memory that can be used by trinity|
|--CPU|the number of CPUs used for the job (should be the same as in the PBS)|
|--trimmomatic| this will trim the adapter sequences left|
|--normalize_reads|uses program Silico to normalise the reads|
|--jaccard_clip|minimizes the falsely fused transcriptes in UTR positions |
|--output|there you put a path to the output dirrectory|
