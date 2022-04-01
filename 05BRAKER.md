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
