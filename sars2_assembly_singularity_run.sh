#!/bin/bash
#================================================================
# HEADER
#================================================================
#% USAGE
#+    bash sars2_assembly_singularity_run.sh reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa iam_sarscov:0.0.4
#%
#% DESCRIPTION
#%      This script run the sars2_assembly_singularity inside the singularity enviroment
#%    
#================================================================
# END_OF_HEADER
#================================================================

#================================================================
# ARGUMENTS
#================================================================

#================================================================
# ARGUMENTS
#================================================================

FASTA=$1 #reference genome
FASTQ1=$2 #foward reads
FASTQ2=$3 #reverse reads
PREFIXOUT=$4 #prefix for output
THREADS=$5 #number of threads
DEPTH=$6 #minum depth to mask regions
MIN_LEN=$7 #minimum length to trimm reads
ADAPTERS=$8 #fasta file with adapters
IMAGE=$9

singularity run --env REFERENCE=$FASTA --env FASTQ1=$FASTQ1 --env FASTQ2=$FASTQ2 --env PREFIXOUT=$PREFIXOUT --env THREADS=$THREADS --env DEPTH=$DEPTH --env MIN_LEN=$MIN_LEN --env ADAPTERS=$ADAPTERS --writable $IMAGE