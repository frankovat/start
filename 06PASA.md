# PASA pipeline

First, format transcriptomes from Trinity. You have to run it twice, since the script didn't work combined

```
#For denovo, you have to name the input file z Trinity YourSpecies_denovo.fasta
from Bio import SeqIO
from glob import glob
import sys

species = sys.argv[1]

denovos = open("%s_denovo.fasta" % (species), 'w')
#ggs = open("%s_gg.fasta" % (species), 'w')
denovo_counter = 0
#gg_counter = 0
for transcripts in glob("*Trinity*.fasta"):
    print(transcripts)
    #if "_GG" in transcripts:
    reader = SeqIO.parse(transcripts, format = 'fasta')
    for rec in reader:
    	denovos.write(">%s\n%s\n" % (rec.id.replace("TRINITY", "TRINITY" + str(denovo_counter)), str(rec.seq)))
    denovo_counter += 1
    #else:
     #   reader = SeqIO.parse(transcripts, format = 'fasta')
      #  for rec in reader:
       #     ggs.write(">%s\n%s\n" % (rec.id.replace("TRINITY_GG", "TRINITY_GG" + str(gg_counter)), str(rec.seq)))
       # gg_counter += 1
denovos.close()
#ggs.close()
```

```
#this is for the genome guided trinity, name your input file YourSpecies_gg.fasta
from Bio import SeqIO
from glob import glob
import sys

species = sys.argv[1]

#denovos = open("%s_denovo.fasta" % (species), 'w')
ggs = open("%s_gg.fasta" % (species), 'w')
#denovo_counter = 0
gg_counter = 0
for transcripts in glob("*Trinity*.fasta"):
    print(transcripts)
    #if "_GG" in transcripts:
    reader = SeqIO.parse(transcripts, format = 'fasta')
    for rec in reader:
    	ggs.write(">%s\n%s\n" % (rec.id.replace("TRINITY", "TRINITY" + str(gg_counter)), str(rec.seq)))
    gg_counter += 1
    #else:
     #   reader = SeqIO.parse(transcripts, format = 'fasta')
      #  for rec in reader:
       #     ggs.write(">%s\n%s\n" % (rec.id.replace("TRINITY_GG", "TRINITY_GG" + str(gg_counter)), str(rec.seq)))
       # gg_counter += 1
#denovos.close()
ggs.close()
```

Ok, now we prepared our transcriptomes (made the name of sequences unique), we can continue with actual PASA
```
#copy all files you will need
cp alignAssembly.config #can be found further down 
cd $SCRATCHDIR

export PATH=$PATH:/software/blat/35/bin
export PATH=$PATH:/software/fasta-36.3.5c/bin
export PERL5LIB=/software/bioperl-1.6.1/lib/perl/5.10.1:$PERL5LIB
export PATH=$PATH:/software/transdecoder/3.0.1

sed -i 's;\$SCRATCHDIR;'"$SCRATCHDIR;" alignAssembly.config
export PASACONF=/afs/ics.muni.cz/software/pasa_2.3.3/PASApipeline-v2.3.3/pasa_conf/sample_test.conf
module load pasa-2.3.3
module load samtools-1.6

species=YourSpecies #the name must be the same as the name of your files

cat  ${species}_denovo.fasta ${species}_gg.fasta > ${species}_transcripts.fasta

software/pasa_2.3.3/pathTOpasaScriptSeqClean/seqclean ${species}_transcripts.fasta

$PASA_HOME/misc_utilities/accession_extractor.pl < $SCRATCHDIR/Cchal_denovo.fasta > tdn.accs

/software/pasa_2.3.3/PASApipeline-v2.3.3/Launch_PASA_pipeline.pl -c alignAssembly.config -C -r -R -g ${species}_gapcloser.fasta.masked -t -u ${species}_transcripts.fasta --ALIGNERS blat,gmap --CPU 16 -I 100000 --TDN tdn.accs --TRANSDECODER
```

alignAssembly.config
```
## templated variables to be replaced exist as <__var_name__>

# database settings
DATABASE=$SCRATCHDIR/YourSpecies.sqlite #this creates the file, the database, you just name it


#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter" 
# assign a value as done above.

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=95
validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY=0

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50

```
