# Trinity - transcriptome preparation
Trinity, a RNA-seq data assembler, is available at [Metacentrum](https://wiki.metacentrum.cz/wiki/Trinity).
We will need the assembled transcriptomes for all gene prediction programs, that is why we are preparing them now.
The [trinity manual](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running%20Trinity) gives you enough answers tu successfully run the program. 
However, I will be sharing with you a script from Priceton, which is easy to adjust to our computational resources.

```
#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=150000

Trinity --seqType fq --left RNA+string.fastq.gz --right RNA-string.fastq.gz --max_memory 100G --CPU 10 --trimmomatic --normalize_reads --jaccard_clip --output outputDir
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
