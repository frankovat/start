# Making our sequence neat!
## What is RepeatModeler
RepeatModeler identifies transposable elements in your genome. You should always clean your genome and find transposable elements and repetitive elements. They can cause many problems when finding genes and gene areas in your genome. This is always the first step!
RepeatModeler is quite easy to use. You can find the official documentation below.

## What is cd-hit:
Cluster Database at High Identity with Tolerance
Produces non-redundant representative sequences as output, reduces overal size of the database without removing any sequence information by only removing redundant or highly similar sequences. It basically cleans the database to make it easier to predict genes.

## Job example:
```bash
#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=20000
#SBATCH --mail-user=frankova@princeton.edu
#SBATCH --job-name=Cnig_repeats
#SBATCH --output=/Genomics/grid/users/frankova/RepeatModeler/Cnig/repeats-%j.out
#SBATCH --error=/Genomics/grid/users/frankova/RepeatModeler/Cnig/repeats-%j.err

module add RepeatModeler
export PATH=$PATH:/Genomics/grid/users/frankova/software/ncbi-blast-2.5.0+/bin/

BuildDatabase -name Cnig_db -engine ncbi /Genomics/grid/users/frankova/Genome/Cnig.fasta 

RepeatModeler -database Cnig_db -pa 3 -engine ncbi >& Cnig.out

/Genomics/grid/users/frankova/software/cdhit/cd-hit-est -i ./RM_*/consensi.fa -o rm_noredun.fa -c 0.8 -M 20000 -n 5 -aS 0.8 -r 1 -T 10

blastx -query rm_noredun.fa -db /Genomics/grid/users/frankova/Reference-genomes/uniprot/uniprot_sprot_blastdb -outfmt '6 qseqid sseqid pident length evalue bitscore qlen slen' -out rm_noredun_against_uniprot.txt -num_threads 10 

blastx -query rm_noredun.fa -db /Genomics/grid/users/frankova/Reference-genomes/drosophila/dmel_blastdb -outfmt '6 qseqid sseqid pident length evalue bitscore qlen slen' -out rm_noredun_against_dmel.txt -num_threads 10

python /Genomics/grid/users/frankova/Python/filter_repeats.py rm_noredun.fa rm_noredun_noprot.fa rm_noredun_against_uniprot.txt rm_noredun_against_dmel.txt
```
`#!/bin/bash`
is a mandatory line when using bash. It acts as a language interpreter, it means, that you are telling the system that your code uses the bash shell scripting language. <br/>

`#SBATCH -n 1` with `#SBATCH`, we are adding parameters for slurm (cluster management and job scheduling system Princeton is using). With `n 1`, we are telling the system that we want to use 1 node for this computation. <br/>

`#SBATCH --cpus-per-task=10` specifies how many CPUs we want to use. Basically the more CPUs the faster the computation. <br/>

`#SBATCH --time=6-23:00 --qos=1wk` specifies time, 6 means six days -23:00 plus 23 hours, qos= <br/>

`#SBATCH --mem=20000` how much memory we want to allocate for this computation in megabytes, this means we want 20000 megabytes <br/>

`#SBATCH --mail-user=frankova@princeton.edu` here you specify your email address. You can add `--mail-type=BEGIN/FAIL/END/REQUEUE or ALL` if you want to get notified about the progress your job is making.<br/>

`#SBATCH --job-name=Cnig_repeats` assigns a name for your job. You can name it accordingly.<br/>

`#SBATCH --output=/Genomics/grid/users/frankova/RepeatModeler/Cnig/repeats-%j.out` creates an output file. `%j` assigns the output file the number of your job in which it runs.<br/>

`#SBATCH --error=/Genomics/grid/users/frankova/RepeatModeler/Cnig/repeats-%j.err` creates an error file. You will damn need this when you want to debug your script. <br/>

`module add RepeatModeler` since RepeatModeler is installed in gen-comp2, you call this software like this. <br/>

`export PATH=$PATH:/Genomics/grid/users/frankova/software/ncbi-blast-2.5.0+/bin/` I installed my own ncbi-blast, if you want to use your own, you export the path to it like this. Always call bin since this is where your exes are.<br/>

`BuildDatabase -name Cnig_db -engine ncbi /Genomics/grid/users/frankova/Genome/Cnig.fasta` for the repeatmodeler to work, you need to create a database. This means that you will use your organism's fasta file which will be used to create the database for RepeatModeler<br/>
|Parameter|Function|
|---------|--------|
|-name|Names your database|
|-engine|Which engine you are going to use, see documentation|<br/>

`RepeatModeler -database Cnig_db -pa 3 -engine ncbi >& Cnig.out` Here we call the program itself which will identify the TEs.<br/>
|Parameter|Function|
|---------|--------|
|-database|name of the database we designed|
|-pa|number of cores we want to use for the search of repetetive elements with our engine|
|-engine|choose desired engine (f.e. NCBI)|
|>&|creates a new file, in our case Cnig.out|<br/>

`/Genomics/grid/users/frankova/software/cdhit/cd-hit-est -i ./RM_*/consensi.fa -o rm_noredun.fa -c 0.8 -M 20000 -n 5 -aS 0.8 -r 1 -T 10` In this line, we call the executable of cd-hit for nucleotide sequences (=est). To use the output from RepeatModeler, we have to clean it a bit more.<br/>
|Parameter|Function|
|---------|--------|
|-i consensi.fa|input file from RepeatModeler, always use the full path|
|-o rm_noredun.fa|output file that is created by the cd-hit-est, consists of a fasta file of representative sequences and also a text file of the list of clusters|
|-c 0.8|sequence identity threshold, default is 0.9=number of identical AAs in alignment divided by the full length of the shorter sequence, number determined by Ben for Halictids|
|-M 20000|max available memory in megabytes, default is 400|
|-n 5|word length, default is 5, see user's guide|
|-aS 0.8|shorter sequence, default is 0.0, if set to 0.9, the alignment must cover 90% of the sequence. Is set as percentage.|
|-r 1|0=default, 1=comparing both strands, ++ and +-|
|-T 10|default 1, 0= all cpus will be used.|<br/>

`blastx -query rm_noredun.fa -db /Genomics/grid/users/frankova/Reference-genomes/uniprot/uniprot_sprot_blastdb -outfmt '6 qseqid sseqid pident length evalue bitscore qlen slen' -out rm_noredun_against_uniprot.txt -num_threads 10` blastx=searches protein database using translated nucleotide query. So we use our translated sequences (now protein sequences) and find protein alignments. </br>
|Parameter|Function|
|---------|--------|
|-query|output file from cd-hit|
|-database|use desired database of protein sequences|
|-outfmt|set parameters for your output, see http://www.metagenomics.wiki/tools/blast/blastn-output-format-6|
|-out|output file|
|-num_threads| number of cpus we want to use|</br>

`blastx -query rm_noredun.fa -db /Genomics/grid/users/frankova/Reference-genomes/drosophila/dmel_blastdb -outfmt '6 qseqid sseqid pident length evalue bitscore qlen slen' -out rm_noredun_against_dmel.txt -num_threads 10` is very similar to the previous line. However! **-database** is going to be an outgroup to our species. 

`python /Genomics/grid/users/frankova/Python/filter_repeats.py rm_noredun.fa rm_noredun_noprot.fa rm_noredun_against_uniprot.txt rm_noredun_against_dmel.txt` Here we call a python script that Ben wrote.</br>

### Python script for filtering repeats by Ben Rubin
```python
from Bio import SeqIO
#nimport sys library
import sys 
#Sys.Argv creates a string of arguments used in your script, 0 is the name of the script, 1-n is the arguments
rm_file = sys.argv[1] 
outfile = sys.argv[2]
uniprot_blast = sys.argv[3]
dmel_blast = sys.argv[4]
#name of your variable=open(file, "r=open for reading, U=universal newlines, expects all types of line endings)
blast_file = open(uniprot_blast, 'rU')
#create a blank list in which you'll assign values
sim_seq_list = []
#for each line in our blast_file
for line in blast_file:
    #variable sequence id = first line [0] that ends with space will be split from the file ()
    seqid = line.split()[0]
    #variable pident, float=gives you decimal value, can transform even string into numbers, you get a number from line 2
    pident = float(line.split()[2])
    #the same for length
    length = float(line.split()[3])
    #the same for bitscore
    bitscore = float(line.split()[5])
    #same for qlength
    qlen = float(line.split()[6])
    #same for slen
    slen = float(line.split()[7])
    #if the bitscore is more or equals to 100 and we do not have the sequence id in the list of sequences
    if bitscore >= 100 and seqid not in sim_seq_list:
        #we add our sequence id into the list of sequences
        sim_seq_list.append(seqid)
     #or if the length/qlength is smaller or equals to 0.50 and pident is more than 50 and seqid is is not in the list of sequences
    elif length / qlen >= 0.50 and pident > 50 and seqid not in sim_seq_list:
        #add the sequence id into the list of sequences
        sim_seq_list.append(seqid)
    #or if the length/slength is bigger or equals to 0.50 and pident is bigger than 50 and seqid is not in the list of sequences
    elif length / slen >= 0.50 and pident > 50 and seqid not in sim_seq_list:
        #add the seqid into the list of sequences
        sim_seq_list.append(seqid)
#the if function went through each line and the blast file is done, the program closes the blast_file
blast_file.close()
#here we do exactly the same for the D. melanogaster blast file
blast_file = open(dmel_blast, 'rU')
for line in blast_file:
    seqid = line.split()[0]
    pident = float(line.split()[2])
    length = float(line.split()[3])
    bitscore = float(line.split()[5])
    qlen = float(line.split()[6])
    slen = float(line.split()[7])
    if bitscore >= 100 and seqid not in sim_seq_list:
        sim_seq_list.append(seqid)
    elif length / qlen >= 0.50 and pident > 50 and seqid not in sim_seq_list:
        sim_seq_list.append(seqid)
    elif length / slen >= 0.50 and pident > 50 and seqid not in sim_seq_list:
        sim_seq_list.append(seqid)

blast_file.close()

seq_file = SeqIO.parse(rm_file, format = 'fasta')
#open a file called outfile for writing = 'w' (reading='r'))
outfile = open(outfile, 'w')
#for each record in sequence file, 
for rec in seq_file:
    if rec.id not in sim_seq_list and len(rec.seq.tostring()) >= 80:
        outfile.write(">%s\n%s\n" % (rec.id, rec.seq.tostring()))

outfile.close()
```
## Documentations
### RepeatModeler
https://blaxter-lab-documentation.readthedocs.io/en/latest/repeatmodeler.html
```BuildDatabase:
NAME
    BuildDatabase - Format FASTA files for use with RepeatModeler

SYNOPSIS
      BuildDatabase [-options] -name "mydb" <seqfile(s) in fasta format>
     or
      BuildDatabase [-options] -name "mydb"
                                  -dir <dir containing fasta files *.fa, *.fasta,
                                         *.fast, *.FA, *.FASTA, *.FAST, *.dna,
                                         and  *.DNA >
     or
      BuildDatabase [-options] -name "mydb"
                                  -batch <file containing a list of fasta files>

DESCRIPTION
      This is basically a wrapper around AB-Blast's and NCBI Blast's
      DB formating programs.  It assists in aggregating files for processing
      into a single database.  Source files can be specified by:

          - Placing the names of the FASTA files on the command
            line.
          - Providing the name of a directory containing FASTA files
            with the file suffixes *.fa or *.fasta.
          - Providing the name of a manifest file which contains the
            names of FASTA files ( fully qualified ) one per line.

      NOTE: Sequence identifiers are not preserved in this database. Each
            sequence is assigned a new GI ( starting from 1 ).  The
            translation back to the original sequence is preserved in the
            *.translation file.

    The options are:

    -h(elp)
        Detailed help

    -name <database name>
        The name of the database to create.

    -engine <engine name>
        The name of the search engine we are using. I.e abblast/wublast or
        ncbi (rmblast version).

    -dir <directory>
        The name of a directory containing fasta files to be processed. The
        files are recognized by their suffix. Only *.fa and *.fasta files
        are processed.

    -batch <file>
        The name of a file which contains the names of fasta files to
        process. The files names are listed one per line and should be fully
        qualified.

SEE ALSO
        RepeatModeler, ABBlast, NCBIBlast

COPYRIGHT
    Copyright 2004-2010 Institute for Systems Biology

AUTHOR
    Robert Hubley <rhubley@systemsbiology.org>


RepeatModeler
NAME
    RepeatModeler - Model repetitive DNA

SYNOPSIS
      RepeatModeler [-options] -database <XDF Database>

DESCRIPTION
    The options are:

    -h(elp)
        Detailed help

    -database
        The prefix name of a XDF formatted sequence database containing the
        genomic sequence to use when building repeat models. The database
        may be created with the WUBlast "xdformat" utility or with the
        RepeatModeler wrapper script "BuildXDFDatabase".

    -engine <abblast|wublast|ncbi>
        The name of the search engine we are using. I.e abblast/wublast or
        ncbi (rmblast version).

    -pa #
        Specify the number of shared-memory processors available to this
        program. RepeatModeler will use the processors to run BLAST searches
        in parallel. i.e on a machine with 10 cores one might use 1 core for
        the script and 9 cores for the BLAST searches by running with "-pa
        9".

    -recoverDir <Previous Output Directory>
        If a run fails in the middle of processing, it may be possible
        recover some results and continue where the previous run left off.
        Simply supply the output directory where the results of the failed
        run were saved and the program will attempt to recover and continue
        the run.

SEE ALSO
        RepeatMasker, WUBlast

COPYRIGHT
     Copyright 2005-2014 Institute for Systems Biology

AUTHOR
     Robert Hubley <rhubley@systemsbiology.org>
     Arian Smit asmit@systemsbiology.org

Interpret the results
RepeatModeler produces a voluminous amount of temporary files stored in a directory created at runtime named like:
       RM_<PID>.<DATE> ie. "RM_5098.MonMar141305172005" 
and remains after each run for debugging purposes or for the purpose of resuming runs if a failure occures. At the succesful completion of a run, two files are generated:
       <database_name>-families.fa  : Consensus sequences
       <database_name>-families.stk : Seed alignments
The seed alignment file is in a Dfam compatible Stockholm format and may be uploaded to the Dfam database by submiting the data to help@dfam.org. In the near future we will provide a tool for uploading families directly to the database.
The fasta format is useful for running quick custom library searches using RepeatMasker. Ie.:
       <RepeatMaskerPath>/RepeatMasker -lib <database_name>-families.fa mySequence.fa 
Other files produced in the working directory include:
       RM_<PID>.<DATE>/
          consensi.fa
          families.stk
          round-1/
               sampleDB-#.fa       : The genomic sample used in this round
               sampleDB-#.fa.lfreq : The RepeatScout lmer table
               sampleDB-#.fa.rscons: The RepeatScout generated consensi
               sampleDB-#.fa.rscons.filtered : The simple repeat/low 
                                               complexity filtered
                                               version of *.rscons
               consensi.fa         : The final consensi db for this round
               family-#-cons.html  : A visualization of the model
                                     refinement process.  This can be opened
                                     in web browsers that support zooming.
                                     ( such as firefox ).
                                     This is used to track down problems 
                                     with the Refiner.pl
               index.html          : A HTML index to all the family-#-cons.html
                                     files.
          round-2/
               sampleDB-#.fa       : The genomic sample used in this round
               msps.out            : The output of the sample all-vs-all 
                                     comparison
               summary/            : The RECON output directory
                    eles           : The RECON family output
               consensi.fa         : Same as above
               family-#-cons.html  : Same as above
               index.html          : Same as above
          round-3/
               Same as round-2
           ..
          round-n/
Recover from a failure
If for some reason RepeatModeler fails, you may restart an analysis starting from the last round it was working on. The -recoverDir [ResultDir] option allows you to specify a diretory ( i.e RM_./ ) where a previous run of RepeatModeler was working and it will automatically determine how to continue the analysis
```

### blastx
http://nebc.nox.ac.uk/bioinformatics/docs/blastx.html

### CD-hit
http://www.bioinformatics.org/cd-hit/cd-hit-user-guide.pdf
