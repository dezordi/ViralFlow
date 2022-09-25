#!/bin/bash


# get bash script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR

# user input
organism_name=$1 # Dengue
organism_refseq_code=$2 # NC_001474.2

# hardcoded paths
SNPEFF_CTNR="singularity_snpeff.sif"
SNPEFF_PATH="/opt/conda/share/snpeff-5.0-1/"
EFETCH_CTNR="singularity_edirect.sif"

# add line to snpeff config
echo -e "# $organism_name, version $organism_refseq_code\n$organism_refseq_code.genome: $organism_name" >> $SNPEFF_CTNR/$SNPEFF_PATH/snpEff.config

# build the directory at DB
mkdir -p $SNPEFF_CTNR/$SNPEFF_PATH/data/$organism_refseq_code

# download fasta
singularity exec $EFETCH_CTNR efetch -db nucleotide -id $organism_refseq_code -format gb > $SNPEFF_CTNR/$SNPEFF_PATH/data/$organism_refseq_code/genes.gbk

# build database
singularity exec $SNPEFF_CTNR snpEff build -genbank -v $organism_refseq_code
