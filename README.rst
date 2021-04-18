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

.. code-block:: text

    current_directory/
    ├-sars2_assembly.sh
    ├-major_minor.py
    ├-minor_finder.py
    ├-put_minor.py
    ├-reference.fasta
    ├-code_R1.fastq.gz
    ├-code_R2.fastq.gz
    └-adapters.fasta

* Arguments
    * <REFERENCEGENOME> -   Fasta file with reference genome
    * <001.fastq.gz>    -   Fasqt file with positive sense reads (R1)
    * <002.fastq.gz>    -   Fastq file with negative sense reads (R2)
    * <PREFIX>          -   Prefix string to store results and to rename consensus genome
    * <NUM_THREADS>     -   Number of threads
    * <DEPTH>           -   Minimum depth to mask unanssembled regions
    * <MIN_LEN>         -   Minimum length to trimm sequences
    * <ADAPTERS>        -   Fasta file with adapters used in the sequencing analysis

**Suggestion to paired-end reads with 150 of length:**

.. code:: bash
    
    bash sars2_assembly.sh reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa

**Suggestion to paired-end reads with 75 of length:**

.. code:: bash

    bash sars2_assembly.sh reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 35 adapters.fa

Both of those examples will generate the following results:


.. code-block:: text

    current_directory/
    ├-sars2_assembly.sh
    ├-major_minor.py
    ├-minor_finder.py
    ├-put_minor.py
    ├-reference.fasta
    ├-code_R1.fastq.gz
    ├-code_R2.fastq.gz
    ├-adapters.fasta
    └-prefix_name.results/
     ├-chromosomes.report #tsv file with genomic metrics
     ├-coverage.report #txt file with all assembly metrics
     ├-prefix_name.<R1/R2>.fq.gz #trimmed fastq files
     ├-prefix_name.depthX.fa #consensus defined with iVar
     ├-prefix_name.depthX.all.fa #in case of minor variant detection, this file contain the 2 genome versions (major and minor consensus)
     ├-prefix_name.depthX.fa.nextclade.csv or prefix_name.depthX.all.fa.nextclade.csv #in case of minor variant detection #nextclade csv output
     ├-prefix_name.depthX.fa.pango.csv or prefix_name.depthX.all.fa.pango.csv #in case of minor variant detection #pangolin lineages information
     ├-prefix_name.depthX.fa.bc #bamreadcount output, with all nucleotide frequencies by genomic position
     ├-prefix_name.depth5.fa.bc.fmt.minors.tsv.fmt #tsv file with minor variant informations
     ├-prefix_name.depth5.minor.fa #fasta file with minor consensus genome
     ├-prefix_name.quality.html #html file with quality controll informations
     ├-prefix_name.sorted.bam #sorted bam file
     ├-prefix_name.sorted.bam.bai #index of sorted bam file
     ├-prefix_name.time.txt #time in minutes of each step of analysis.
     └-prefix_name.tsv #tsv output from iVar with the frequencies of iSNVs

=====
Disclaimer
=====
* The fastq files should have in the same directory of sars2_assembly.sh and the python scripts.
* The minor consensus version is only based on minor allele frequencies throught the genome mapping, without any statistical method to reconstruct quasispecies genomic populations.
* If you use this workflow for academic  purposes, please cite this repository;
* More information `Here <https://dezordi.github.io/>`_;