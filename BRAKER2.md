#!/bin/bash 
#PBS -l select=1:ncpus=32:mem=250gb:scratch_local=200gb 
#PBS -l walltime=23:00:00 
#PBS -m abe

HOMEDIR="/auto/brno2/home/frankovat/C_cucurbitina/braker2"


cp /auto/brno2/home/frankovat/C_cucurbitina/RepeatMasker/Ccuc_rm_masked/Ccuc_mira.fa.masked $SCRATCHDIR || exit 1
cp -r /auto/brno2/home/frankovat/Orthofinder/ceratinas $SCRATCHDIR || exit 2
cd $SCRATCHDIR || exit 3


source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module add braker/2.1.6-gcc-10.2.1-reuw2tp
test -r ~/.gm_key || ln -s /software/genemark/4.68/genemark4.68_key_64 ~/.gm_key

mkdir augustus; lndir $(dirname $AUGUSTUS_CONFIG_PATH) augustus >/dev/null
export AUGUSTUS_CONFIG_PATH=$SCRATCHDIR/augustus/config

export GENEMARK_PATH=$(dirname $(which gmes_petap.pl) )

export PROTHINT_PATH=$SCRATCHDIR/ProtHint/bin
cp -r /storage/brno2/home/frankovat/software/ProtHint $SCRATCHDIR || exit 6
chmod -R u=rwX $PROTHINT_PATH

export ALIGNMENT_TOOL_PATH=$SCRATCHDIR/spaln
cp -r /storage/brno2/home/frankovat/software/download/spaln-master/bin $ALIGNMENT_TOOL_PATH || exit 7
chmod -R u=rwX $ALIGNMENT_TOOL_PATH

braker.pl --genome=$SCRATCHDIR/Ccuc_mira.fa.masked --prot_seq=$SCRATCHDIR/ceratinas/Cchal_pep.fa,$SCRATCHDIR/ceratinas/Ccuc_pep.fa,$SCRATCHDIR/ceratinas/Ccyp_pep.fa,$SCRATCHDIR/ceratinas/Cmor_pep.fa \
 --species=Ccuc_braker2 --gff3 --cores=32 --softmasking --overwrite --epmode --AUGUSTUS_ab_initio

cp -r $SCRATCHDIR/braker $HOMEDIR/${PBS_JOBID}-braker && clean_scratch
