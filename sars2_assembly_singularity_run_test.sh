#!/bin/bash
#================================================================
# HEADER
#================================================================
#% USAGE
#%    bash sars2_assembly_singularity_run.sh reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa iam_sarscov:0.0.5
#%    bash sars2_assembly_singularity_run.sh reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa iam_sarscov:0.0.5 50 5
#%
#% DESCRIPTION
#%    This script run the sars2_assembly_singularity inside the singularity enviroment
#%
#================================================================
# END_OF_HEADER
#================================================================

#================================================================
# ARGUMENTS
#================================================================
INPUT_DIR=$1 # input dir path
# WARNING: All remaining arguments must be define in relation to INPUT_DIR
# The container assumes all input files are inside the specified dir.
# ex: if /my/input/dir/fastq_R1.fastq.gz, then FASTQ1 should be fastq_R1.fastq.gz
#     if /my/input/dir/reference/ref.fasta, then FASTA must be referece/ref.fasta
FASTA=$2 #reference genome

#FASTQ1=$3 #foward reads
#FASTQ2=$4 #reverse reads
#PREFIXOUT=$5 #prefix for output

ADAPTERS=$3 #fasta file with adapters
THREADS=$4 #number of threads
DEPTH=$5 #minum depth to mask regions
MIN_LEN=$6 #minimum length to trimm reads
IMAGE=$7 # container image file
DP_INTRAHOST=${8-100} #minimum dp for intrahost analysis
TRIMM_LEN=${9-0} #length to trim front and tail reads with fastp

singularity run --bind $INPUT_DIR:/data/ --env REFERENCE=$FASTA --env THREADS=$THREADS --env DEPTH=$DEPTH --env MIN_LEN=$MIN_LEN --env ADAPTERS=$ADAPTERS --env DP_INTRAHOST=$DP_INTRAHOST --env TRIMM_LEN=$TRIMM_LEN --writable-tmpfs $IMAGE