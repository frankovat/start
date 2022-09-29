## 01Maker_reformatgff
```
#!/bin/bash
#PBS -l select=1:ncpus=8:ompthreads=8:mem=300gb:scratch_local=150gb
#PBS -l walltime=5:00:00
#PBS -m abe

cp /auto/plzen1/home/frankovat/PASA/job_12326530.meta-pbs.metacentrum.cz/job_12299021.meta-pbs.metacentrum.cz/Cchal.sqlite.pasa_assemblies.gff3 $SCRATCHDIR || exit  1

module add python-3.6.2-gcc

python /auto/plzen1/home/frankovat/Maker/reformat_gff.py /auto/plzen1/home/frankovat/PASA/job_12326530.meta-pbs.metacentrum.cz/job_12299021.meta-pbs.metacentrum.cz/Cchal.sqlite.pasa_assemblies.gff3 /auto/plzen1/home/frankovat/Maker/Cchal.sqlite.pasa_assemblies_form.gff3 Cchal #transforming previous gff file into new one
```

## 02Maker_ifloop
```
#!/bin/bash
#PBS -l select=1:ncpus=14:mem=300gb:scratch_local=150gb
#PBS -l walltime=70:00:00
#PBS -m abe

cp /auto/plzen1/home/frankovat/Maker/maker_bopts.ctl $SCRATCHDIR || exit 1
cp /auto/plzen1/home/frankovat/Maker/maker_opts.ctl $SCRATCHDIR || exit 2
cp /auto/plzen1/home/frankovat/Maker/maker_exe.ctl $SCRATCHDIR || exit 3
cp /auto/plzen1/home/frankovat/Maker/Cchal.sqlite.pasa_assemblies_form.gff3 $SCRATCHDIR || exit 4
cp -r /auto/plzen1/home/frankovat/Maker/02Output/Cchal.maker.output $SCRATCHDIR || exit 5
cd $SCRATCHDIR || exit 6


module load maker-2.31.10

mkdir tmp

mpirun maker -fix_nucleotides -base Cchal -TMP "$SCRATCHDIR/tmp" &
sleep 20m

mpirun maker -fix_nucleotides -base Cchal -TMP "$SCRATCHDIR/tmp"&

rm -r -f $SCRATCHDIR/tmp
#force removal of a directory

mkdir round_1_output_brakerpred
cd round_1_output_brakerpred

fasta_merge -d $SCRATCHDIR/Cchal.maker.output/Cchal_master_datastore_index.log -o Cchal
gff3_merge -d $SCRATCHDIR/Cchal.maker.output/Cchal_master_datastore_index.log -o Cchal_all.gff3
gff3_merge -d $SCRATCHDIR/Cchal.maker.output/Cchal_master_datastore_index.log -o Cchal_genes.gff3 -g

cp -r $SCRATCHDIR /storage/plzen1/home/frankovat/Maker || export CLEAN_SCRATCH=false
```

## 03Maker_busco1
```
#!/bin/bash
#PBS -l select=1:ncpus=14:mem=300gb:scratch_local=250gb
#PBS -l walltime=15:00:00
#PBS -m abe

cp -r /storage/plzen1/home/frankovat/Maker/02Output/round_1_output_brakerpred $SCRATCHDIR || exit1
cp -r /storage/plzen1/home/frankovat/GenomesTranscriptomes/MakerFiles/hymenoptera_odb9 $SCRATCHDIR || exit 4
cd $SCRATCHDIR || exit 2

#module load maker-2.31.10
#module load mpich-3.0.2-pgi
module load busco-3.0.2

cd round_1_output_brakerpred

run_BUSCO.py -i $SCRATCHDIR/round_1_output_brakerpred/Cchal.all.maker.proteins.fasta -o busco_prots -l $SCRATCHDIR/hymenoptera_odb9 -m proteins -c 14 -sp Cchal_braker

cp -r $SCRATCHDIR /storage/plzen1/home/frankovat/Maker/ || export CLEAN_SCRATCH=false
```

## 04Maker_interproscan
```
#!/bin/bash
#PBS -l select=1:ncpus=9:mem=300gb:scratch_local=250gb
#PBS -l walltime=15:00:00
#PBS -m abe

cp -r /auto/plzen1/home/frankovat/Maker/02Output/round_1_output_brakerpred $SCRATCHDIR || exit 1
cd $SCRATCHDIR || exit 2

module load interproscan-5.55-88.0

cd round_1_output_brakerpred

interproscan.sh -T $SCRATCHDIR -i Cchal.all.maker.non_overlapping_ab_initio.proteins.fasta --seqtype p -goterms -iprlookup -b Cchal.all.maker.non_overlapping_ab_initio.proteins.fasta.iprscan -pa

grep "IPR" Cchal.all.maker.non_overlapping_ab_initio.proteins.fasta.iprscan.tsv | awk '{print $1}' | sort | uniq > ipr_hits.txt

cp -r $SCRATCHDIR /storage/plzen1/home/frankovat/Maker || export CLEAN_SCRATCH=false
```

## 05Maker_gff3-select
```
#!/bin/bash
#PBS -l select=1:ncpus=9:mem=300gb:scratch_local=250gb
#PBS -l walltime=15:00:00
#PBS -m abe

cp /auto/plzen1/home/frankovat/Maker/02Output/round_1_output_brakerpred/Cchal_all.gff3 $SCRATCHDIR || exit 1
cp /auto/plzen1/home/frankovat/Maker/maker_bopts.ctl $SCRATCHDIR || exit 2
cp /auto/plzen1/home/frankovat/Maker/maker_opts.ctl $SCRATCHDIR || exit 3
cp /auto/plzen1/home/frankovat/Maker/maker_exe.ctl $SCRATCHDIR || exit 4
cp /auto/plzen1/home/frankovat/Maker/maker_opts_ipradd.ctl $SCRATCHDIR || exit 5
cp /auto/plzen1/home/frankovat/Maker/gff3_select $SCRATCHDIR || exit 6
cd $SCRATCHDIR || exit 7

module load maker-2.31.10

perl gff3_select Cchal_all.gff3 ipr_hits.txt > ipr_hits.gff3

mkdir "/tmp/05Maker"

mpirun maker -fix_nucleotides -base Cchal_ipradd -TMP "/tmp/05Maker" maker_opts_ipradd.ctl

rm -r /tmp/05Maker

cp -r $SCRATCHDIR /storage/plzen1/home/frankovat/Maker || export CLEAN_SCRATCH=false
```

## 06Maker_gff3merge2
```
#!/bin/bash
#PBS -l select=1:ncpus=9:mem=300gb:scratch_local=250gb
#PBS -l walltime=15:00:00
#PBS -m abe

#cp /auto/plzen1/home/frankovat/Maker/02Output/round_1_output_brakerpred/Cchal_all.gff3 $SCRATCHDIR || exit 1
cp -r  /auto/plzen1/home/frankovat/Maker/02Output/Cchal_ipradd.maker.output $SCRATCHDIR || exit 1
cp /auto/plzen1/home/frankovat/Maker/maker_bopts.ctl $SCRATCHDIR || exit 2
#cp /auto/plzen1/home/frankovat/Maker/maker_opts.ctl $SCRATCHDIR || exit 3
cp /auto/plzen1/home/frankovat/Maker/maker_exe.ctl $SCRATCHDIR || exit 3
cp /auto/plzen1/home/frankovat/Maker/maker_opts_ipradd.ctl $SCRATCHDIR || exit 4
cd $SCRATCHDIR || exit 5

module load maker-2.31.10
#module load bioperl-1.6.9-gcc
#module load boost-1.60-gcc
#module load hmmer-3.3.2
#module load blast+-2.8.0a-src

mkdir tmp

mpirun maker -fix_nucleotides -base Cchal_ipradd -TMP "$SCRATCHDIR/tmp" &
sleep 20m
 
mpirun maker -fix_nucleotides -base Cnig -TMP "$SCRATCHDIR/tmp"&

rm -r -f $SCRATCHDIR/tmp
#force removal of a directory

mkdir round_2_output_ipradd
cd round_2_output_ipradd

fasta_merge -d $SCRATCHDIR/Cchal_ipradd.maker.output/Cchal_ipradd_master_datastore_index.log -o cchal

gff3_merge -d $SCRATCHDIR/Cchal_ipradd.maker.output/Cchal_ipradd_master_datastore_index.log -o cchal_all.gff3

gff3_merge -d $SCRATCHDIR/Cchal_ipradd.maker.output/Cchal_ipradd_master_datastore_index.log -o cchal_genes.gff3 -g

maker_map_ids --prefix cchal_ --justify 5 cchal_genes.gff3 > Cchal_genes.ids

cp -r $SCRATCHDIR/round_2_output_ipradd /auto/plzen1/home/frankovat/Maker/02Output || CLEAN_SCRATCH=false
```


# CTL files

### maker_opts
```
#-----Genome (these are always required)
genome=/auto/plzen1/home/frankovat/RepeatMasker/Cchal/Cchal_rm_masked/Cchal_gapcloser.fasta.masked #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=/auto/plzen1/home/frankovat/Maker/Cchal.sqlite.pasa_assemblies_form.gff3 #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=/auto/plzen1/home/frankovat/GenomesTranscriptomes/MakerFiles/amel_OGSv3.2_pep.fa,/auto/plzen1/home/frankovat/GenomesTranscriptomes/MakerFiles/uniprot_trembl.fasta, /storage/plzen1/home/frankovat/GenomesTranscriptomes/uniprot_amel.fasta #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff= #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= arthropoda #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/afs/ics.muni.cz/software/maker/2.31.10/maker/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff=/auto/plzen1/home/frankovat/BRAKER/Cchal/job_11820497.meta-pbs.metacentrum.cz/braker/augustus.ab_initio.gff3 #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
#run_evm=1 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
#trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
#snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
#snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=300000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=3000 #skip genome contigs below this length (under 10kb are often useless)
pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=1 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=100000 #length for the splitting of hits (expected max intron size for evidence alignments)
#min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=1 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```

### maker_exe.ctl
```
#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=/software/blastPlus-2.2.29/bin/makeblastdb #location of NCBI+ makeblastdb executable
blastn=/software/blastPlus-2.2.29/bin/blastn #location of NCBI+ blastn executable
blastx=/software/blastPlus-2.2.29/bin/blastx #location of NCBI+ blastx executable
tblastx=/software/blastPlus-2.2.29/bin/tblastx #location of NCBI+ tblastx executable
formatdb= #location of NCBI formatdb executable
blastall= #location of NCBI blastall executable
xdformat= #location of WUBLAST xdformat executable
blasta= #location of WUBLAST blasta executable
#prerapsearch= #location of prerapsearch executable
#rapsearch= #location of rapsearch executable
RepeatMasker=/software/RepeatMasker/RepeatMasker-4-0-0a/RepeatMasker #location of RepeatMasker executable
exonerate=/software/exonerate/2.2.0/bin/exonerate #location of exonerate executable

#-----Ab-initio Gene Prediction Algorithms
snap=/software/snap-1.0.16 #location of snap executable
gmhmme3= #location of eukaryotic genemark executable
gmhmmp= #location of prokaryotic genemark executable
augustus=/software/augustus/3.4.0/bin/augustus #location of augustus executable
fgenesh= #location of fgenesh executable
#evm=/storage/plzen1/home/frankovat/EVidenceModeler-1.1.1/evidence_modeler.pl #location of EvidenceModeler executable
#tRNAscan-SE= #location of trnascan executable
#snoscan= #location of snoscan executable

#-----Other Algorithms
probuild= #location of probuild executable (required for genemark)
```

### maker_bopts.ctl
```
#-----BLAST and Exonerate Statistics Thresholds
blast_type=ncbi+ #set to 'ncbi+', 'ncbi' or 'wublast'
#use_rapsearch=0 #use rapsearch instead of blastx, 1 = yes, 0 = no
pcov_blastn=0.8 #Blastn Percent Coverage Threhold EST-Genome Alignments
pid_blastn=0.85 #Blastn Percent Identity Threshold EST-Genome Aligments
eval_blastn=1e-10 #Blastn eval cutoff
bit_blastn=40 #Blastn bit cutoff
depth_blastn=0 #Blastn depth cutoff (0 to disable cutoff)

pcov_blastx=0.5 #Blastx Percent Coverage Threhold Protein-Genome Alignments
pid_blastx=0.4 #Blastx Percent Identity Threshold Protein-Genome Aligments
eval_blastx=1e-06 #Blastx eval cutoff
bit_blastx=30 #Blastx bit cutoff
depth_blastx=0 #Blastx depth cutoff (0 to disable cutoff)

pcov_tblastx=0.8 #tBlastx Percent Coverage Threhold alt-EST-Genome Alignments
pid_tblastx=0.85 #tBlastx Percent Identity Threshold alt-EST-Genome Aligments
eval_tblastx=1e-10 #tBlastx eval cutoff
bit_tblastx=40 #tBlastx bit cutoff
depth_tblastx=0 #tBlastx depth cutoff (0 to disable cutoff)

pcov_rm_blastx=0.5 #Blastx Percent Coverage Threhold For Transposable Element Masking
pid_rm_blastx=0.4 #Blastx Percent Identity Threshold For Transposbale Element Masking
eval_rm_blastx=1e-06 #Blastx eval cutoff for transposable element masking
bit_rm_blastx=30 #Blastx bit cutoff for transposable element masking

ep_score_limit=20 #Exonerate protein percent of maximal score threshold
en_score_limit=20 #Exonerate nucleotide percent of maximal score threshold
```

### maker_evm.ctl
```
#-----Transcript weights
evmtrans=10 #default weight for source unspecified est/alt_est alignments
evmtrans:blastn=0 #weight for blastn sourced alignments
evmtrans:est2genome=10 #weight for est2genome sourced alignments
evmtrans:tblastx=0 #weight for tblastx sourced alignments
evmtrans:cdna2genome=7 #weight for cdna2genome sourced alignments

#-----Protein weights
evmprot=10 #default weight for source unspecified protein alignments
evmprot:blastx=2 #weight for blastx sourced alignments
evmprot:protein2genome=10 #weight for protein2genome sourced alignments

#-----Abinitio Prediction weights
evmab=30 #default weight for source unspecified ab initio predictions
evmab:snap=10 #weight for snap sourced predictions
evmab:augustus=10 #weight for augustus sourced predictions
evmab:fgenesh=10 #weight for fgenesh sourced predictions
evmab:genemark=7 #weight for genemark sourced predictions
```
### maker_opts_ipradd.ctl
```
#-----Genome (these are always required)
genome=/auto/plzen1/home/frankovat/Postprocessing/Gapcloser/Cchal/Cchal_gapcloser.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=#protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff= #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= $SCRATCHDIR/02Output/round_1_output_brakerpred/ipr_hits.gff3 #ab-initio predictions from an external GFF3 file
model_gff=$SCRATCHDIR/02Output/round_1_output_brakerpred/Cchal_genes.gff3 #annotated gene models from an external GFF3 file (annotation pass-through)
#run_evm=1 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
#trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
#snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
#snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=300000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=3000 #skip genome contigs below this length (under 10kb are often useless)
pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=1 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=100000 #length for the splitting of hits (expected max intron size for evidence alignments)
#min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=1 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```

### Reformat_gff
```
import sys

ingff = sys.argv[1]
outgff = sys.argv[2]
prefix = sys.argv[3]

reader = open(ingff, 'r')

outfile = open(outgff, 'w')

id_list = []
prev_id = ""
matches = 0
for line in reader:
    cur_line = line.split()
    cur_id = cur_line[8].split(";")[0].split("ID=align_")[1]
    if True: #cur_id == prev_id:
        matches += 1
        new_line = line.replace("cDNA_match", "match_part")
        if matches < 10: 
            new_line = new_line.replace("ID=align_%s" % cur_id, "ID=match0000%s;Parent=align_%s" % (matches, cur_id)) #if the gene is a parent to the exon, i
            #making a match id, making it easier to sort it and read it
        else:
            new_line = new_line.replace(cur_id, "ID=match000%s;Parent=align_%s" % (matches, cur_id))
#        outfile.write(new_line.replace(prefix, "%s_scaff" % prefix))
            #I don't think this is necessary anymore so eliminating
        outfile.write(new_line.replace(prefix, "%s" % prefix))
    else:
        matches = 0
        prev_id = cur_id
#        outfile.write(line.replace(prefix, "%s_scaff" % prefix))
        outfile.write(line.replace(prefix, "%s" % prefix))

outfile.close()
#intention: correctly label children and parents

# reader = open(outgff, 'w')
# reader = open(outgff + "_temp", 'rU')

# for line in reader:
#     cur_line = line.split()
#     if cur_line[2] == "cDNA_match":
#         cur_id = cur_line[8].split(";")[0].split("ID=align_")[1]
```

### gff3_select
```
#!/usr/bin/perl
use strict;
use warnings;

my $usage = "
USAGE:
     gff3_select <maker.gff3> <list_of_ids.txt> <filter>

     This script will select a subset set of features from a maker produced
     GFF3 file and return them with inheritance resolved. You can add an
     optional filter to be used in searching column 3 of top level features.

";

my $gff_file = shift;
my $id_file = shift; #file to get IDs from
my $filter = shift;


if(!$gff_file || ! -f $gff_file || ! $id_file || ! -f $id_file){
    print $usage;
    exit();
}

my $all = parse_gff($gff_file);

my $alias = {};
foreach my $f (values %$all){
    if(my $name = $f->{attributes}{Name}){
	$alias->{$name} = $f;
    }
}

my @selected;
open(IN, "<  $id_file");
while(my $line = <IN>){
    chomp $line;
    $line =~ s/^\s+|\s+$//g;
    
    my $f = $all->{$line};
    $f = $alias->{$line} if(!$f);

    if(! $f){
	warn "**WARNING: No top level feature found for ID $line\n";
	next;
    }

    next if($filter && $f->{type} ne $filter);

    dump_feature($f);
}
close(IN);


sub dump_feature{
    my $f = shift;

    if(! $f->{dumped}){
        print $f->{seqid}."\t";
        print $f->{source}."\t";
        print $f->{type}."\t";
        print $f->{start}."\t";
        print $f->{end}."\t";
        print $f->{score}."\t";
        print $f->{strand}."\t";
        print $f->{phase}."\t";

        my @atts;
        while(my $key = each %{$f->{attributes}}){
            next if(! $key);
            my $att = "$key=";
            if(ref($f->{attributes}{$key}) eq 'ARRAY'){
                $att .= join(',', @{$f->{attributes}{$key}});
            }
            else{
                $att .= $f->{attributes}{$key};
            }
            push(@atts, $att);
        }
        print join(';', @atts)."\n";
        $f->{dumped}++;
    }
    else{
	return;
    }

    foreach my $s (@{$f->{children}}){
	dump_feature($s);
    }
}

sub parse_gff {
    my $ann_file = shift;

    my %top_index;
    my %id_index;
    my %famtree;
    open(IN, "< $ann_file");
    while(my $line = <IN>){
	last if($line =~ /^>/);
	next if($line =~ /^\#/);
	chomp($line);
	next if(!$line);

	my @F = split(/\t/, $line);
	next if(@F != 9);

	my %att;
	foreach my $a (split(/\;/, $F[8])){
	    $a =~ s/^\s+//;
	    $a =~ s/\s+$//;
	    next unless($a);
	    my ($key, $value) = split(/=/, $a);

	    if($key =~ /^(ID|Name|Target|Gap|Derives_from|Is_circular)$/){
		$att{$key} = $value;
	    }
	    else{
		$att{$key} = [split(/\,/, $value)];
	    }
	}
	
	my %f = (seqid => $F[0],
		 source => $F[1],
		 type => $F[2],
		 start => $F[3],
		 end => $F[4],
		 score => $F[5],
		 strand => $F[6],
		 phase => $F[7],
		 attributes => \%att);

	if($f{attributes}{ID} && !$f{attributes}{Parent}){
	    $top_index{$f{attributes}{ID}}  = \%f;
	}
	
	if($f{attributes}{ID}){
	    my $id = $f{attributes}{ID};
	    if($id_index{$id} && $id_index{$id}{children}){
		my $children = $id_index{$id}{children};
		$id_index{$id} = \%f;
		$id_index{$id}{children} = $children;
	    }
	    else{
		$id_index{$id} = \%f;
	    }
	}
	
	if($f{attributes}{Parent}){
	    foreach my $p (@{$f{attributes}{Parent}}){
		push(@{$id_index{$p}{children}}, \%f);
	    }
	}
    }
    close(IN);

    return \%top_index;
}
'''
