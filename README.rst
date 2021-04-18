IAM_SARS-CoV-2
=========

This repository contains a set of scripts to performs a reference guided genome assembly of SARS-CoV-2. Python scripts were developed based on the wuhan SARS-CoV-2 reference genome NC_045512.2. The workflow was developed to work with Illumina paired-end reads. Tests with other technologies should be performed.

=====
Dependencies
=====

* BWA Version: 0.7.17-r1198-dirty
* samtools 1.7 Using htslib 1.7-2
* fastp 0.20.1
* iVar version 1.3.1
* bam-readcount version: 0.8.0-unstable-7-625eea2
* Python 3.7.10
    * argparse 1.1
    * pandas 1.1.3
    * numpy 1.19.2
    * biopython 1.78
* mafft v7.310 (2017/Mar/17)    
* nextclade 0.14.2
* pangolin 2.3.9
* bedtools v2.26.0
* bamdst 1.0.6

=====
Usage
=====

bash sars2_assembly.sh <REFERENCEGENOME> <001.fastq.gz> <002.fastq.gz> <PREFIX> <NUM_THREADS> <DEPTH> <MIN_LEN> <ADAPTERS>

* Arguments
    * <REFERENCEGENOME> -   Fasta file with reference genome
    * <001.fastq.gz>    -   Fasqt file with positive sense reads (R1)
    * <002.fastq.gz>    -   Fastq file with negative sense reads (R2)
    * <PREFIX>          -   Prefix string to store results and to rename consensus genome
    * <NUM_THREADS>     -   Number of threads
    * <DEPTH>           -   Minimum depth to mask unanssembled regions
    * <MIN_LEN>         -   Minimum length to trimm sequences
    * <ADAPTERS>        -   Fasta file with adapters used in the sequencing analysis

**Sugestion to paired-end reads with 150 of length:**

.. code:: bash
    
    bash sars2_assembly.sh reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa

**Sugestion to paired-end reads with 75 of length:**

.. code:: bash

    bash sars2_assembly.sh reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 35 adapters.fa

=====
Disclaimer
=====

* If you use this workflow for academic pourposes, please cite this repository;
* More information `Here <https://dezordi.github.io/>`_;