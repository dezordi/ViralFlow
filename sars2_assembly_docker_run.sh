#!/bin/bash
#================================================================
# HEADER
#================================================================
#% USAGE
#+    $bash sars2_assembly_docker_run.sh <REFERENCEGENOME> <001.fastq.gz> <002.fastq.gz> <PREFIX> <NUM_THREADS> <DEPTH> <MIN_LEN> <ADAPTERS>
#%
#% DESCRIPTION
#%    This script automates the running of sars2_assembly.sh via docker.
#%    
#% OPTIONS
#%    <REFERENCEGENOME>    -   Fasta file with reference genome
#%    <001.fastq.gz>       -   Fasqt file with positive sense reads (R1)
#%    <002.fastq.gz>       -   Fastq file with negative sense reads (R2)
#%    <PREFIX>             -   Prefix string to store results and to rename consensus genome
#%    <NUM_THREADS>        -   Number of threads
#%    <DEPTH>              -   Minimum depth to mask unanssembled regions
#%    <MIN_LEN>            -   Minimum length to trimm sequences
#%    <ADAPTERS>           -   Fasta file with adapters used in the sequencing analysis
#%    <IMAGE>              -   image:tag
#%
#% EXAMPLES
#%    $bash  sars2_assembly_docker_run.sh reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa sars2_assembly:0.0.3
#%
#% DEPENDENCIES
#%    BWA Version: 0.7.17-r1198-dirty
#%    samtools 1.7 Using htslib 1.7-2
#%    fastp 0.20.1
#%    iVar version 1.3.1
#%    bam-readcount version: 0.8.0-unstable-7-625eea2
#%    Python 3.7.10
#%    mafft v7.310 (2017/Mar/17)    
#%    nextclade 0.14.2
#%    pangolin 2.3.9
#%    bedtools v2.26.0
#%    bamdst 1.0.6
#%
#================================================================
#- IMPLEMENTATION
#-    version         $sars2_assembly_docker_run 0.0.1
#-    authors         Filipe Dezordi and Gabriel Wallau
#-    maintainer      Filipe Dezordi (zimmer.filipe@gmail.com)
#-    username        dezordi
#-    license         GPL
#-    information     dezordi.github.io
#================================================================
#  HISTORY
#     2021/06/04 : dezordi : Script creation
# 
#================================================================
# END_OF_HEADER
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

docker run --env REFERENCE=$FASTA --env FASTQ1=$FASTQ1 --env FASTQ2=$FASTQ2 --env PREFIXOUT=$PREFIXOUT --env THREADS=$THREADS --env DEPTH=$DEPTH --env MIN_LEN=$MIN_LEN --env ADAPTERS=$ADAPTERS -v $(pwd)/:/home/ -w /home/ --rm $IMAGE