IAM_SARS-CoV-2
=========

This repository contains a set of scripts to performs a reference guided genome assembly of SARS-CoV-2. Python scripts were developed based on the wuhan SARS-CoV-2 reference genome NC_045512.2. The workflow was developed to work with Illumina paired-end reads. Tests with other technologies should be performed.

.. image:: images/workflow.png
   :width: 600

=====
Dependencies
=====

* BWA Version: 0.7.17-r1188
* samtools 1.9 Using htslib 1.9
* fastp 0.20.1
* iVar version 1.3.1
* bam-readcount version: 0.8.0-unstable-7-625eea2
* Python 3.8.1
    * argparse 1.4
    * pandas 1.0.1
    * numpy 1.20.3
    * biopython 1.74
* mafft v7.310 (2017/Mar/17)    
* nextclade 0.14.2
* pangolin v3.1.4
* bedtools v2.27.0
* bamdst 1.0.6

=====
Files info
=====

.. code-block:: text

    IAM_SARSCOV2/
    ├-Dockerfile                            ### Recipe to build local docker image
    ├-sars2_assembly_docker                 ### Script called into ENTRYPOINT of local docker image
    ├-sars2_assembly_docker_run.sh          ### Script for users unfamiliar with docker run sintaxe 
    ├-Singularityfile                       ### Recipe to build local singularity sandbox
    ├-sars2_assembly_singularity            ### Script called into ENTRYPOINT of local singularity sandbox
    ├-sars2_assembly_singularity_run.sh     ### Script for users unfamiliar with singularity run sintaxe 
    ├-pango_update                          ### Script to activate conda and update pangolin, run automatically during docker or singularity build
    └-python_scripts:                       
      ├-assembly_metrics.py                 ### Run bamdst 
      ├-bwa_index.py                        ### Run bwa index
      ├-bwa_mem.py                          ### Run bwa mem
      ├-fastp.py                            ### Run fastp
      ├-get_mvs.py                          ### Perform intrahost variant analysis with bam-readcount and intrahost.py
      ├-intrahost.py                        ### Identify genomic positions with multi-allele frequencies
      ├-ivar.py                             ### Run ivar variant and ivar consensus
      └-pango_nextclade.py                  ### Run pangolin and nextclade


=====
Docker
=====

A docker image with all tools and libraries can be found `here <https://hub.docker.com/repository/docker/dezordi/iam_sarscov2/>`_.
The last update of the pangolin in the docker images was carried out on June 22, 2021 to the version v3.1.4.
You can create a container and run as an interactive session the sars2_assembly following:

.. code:: bash
    
    docker run -tdi --name iam_sarscov2 --cpus <number> --memory <number> dezordi/iam_sarscov2:0.0.4 /bin/bash
    docker cp  <REFERENCEGENOME> <001.fastq.gz> <002.fastq.gz> <ADAPTERS_FILE> iam_sarscov2:home
    docker attach iam_sarscov2
    cd home
    bash sars2_assembly <REFERENCEGENOME> <001.fastq.gz> <002.fastq.gz> <PREFIX> <NUM_THREADS> <DEPTH> <MIN_LEN> <ADAPTERS_FILE>


* Arguments docker run
    * tdi     -   t and i create an interactive environment similar to terminal connection session, d run the container in background;
    * name    -   container name;
    * cpus    -   number maximum of threads;
    * memory  -   ram memory limit;

Or you can use the Dockerfile and sars2_assembly_docker_run.sh to run the docker without the interactive mode:

.. code:: bash
    
    docker build -t <image>:<tag> .
    bash sars2_assembly_docker_run.sh <REFERENCEGENOME> <001.fastq.gz> <002.fastq.gz> <PREFIX> <NUM_THREADS> <DEPTH> <MIN_LEN> <ADAPTERS_FILE> <image>:<tag>

Using the Dockerfile and sars2_assembly_docker_run.sh a directory named 'prefix.results' will be created in the current directory storing the results.

**Suggestion to paired-end reads with 150 of length using Dockerfile:**

.. code:: bash
    
    docker build -t iam_sarscov2:0.0.4 .
    bash sars2_assembly_docker_run.sh reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa iam_sarscov2:0.0.4

=====
Singularity
=====

For environments with non-root privileges, you can run the analysis using singularity. A recipe file was create using the same docker image.
The recipe file and following steps were tested for singularity version 3.7.1.

.. code:: bash
    
    singularity build --fakeroot <imagename> Singularityfile
    bash sars2_assembly_singularity_run.sh <REFERENCEGENOME> <001.fastq.gz> <002.fastq.gz> <PREFIX> <NUM_THREADS> <DEPTH> <MIN_LEN> <ADAPTERS_FILE> <imagename>

**Suggestion to paired-end reads with 150 of length using Singularity:**

.. code:: bash
    
    singularity build --fakeroot iam_sarscov2.0.0.4 Singularityfile
    bash sars2_assembly_singularity_run.sh reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa iam_sarscov2:0.0.4

=====
Explained Usage
=====

**Into interactive docker container**

bash sars2_assembly <REFERENCEGENOME> <001.fastq.gz> <002.fastq.gz> <PREFIX> <NUM_THREADS> <DEPTH> <MIN_LEN> <ADAPTERS_FILE>

* Arguments
    * <REFERENCEGENOME> -   Fasta file with reference genome.
    * <001.fastq.gz>    -   Fasqt file with positive sense reads (R1).
    * <002.fastq.gz>    -   Fastq file with negative sense reads (R2).
    * <PREFIX>          -   Prefix string to store results and to rename consensus genome.
        * The user can set the gisaid format genome name, and the workflow will automatically format the consensus name, as the prefix will be used to create the directory output, the slash '/' should be replaced by '__' and the pipe '|' should be replaced by '--'.
        * e.g. prefix:       hCoV-19__Brazil__PE-FIOCRUZ-IAM1234__2020--2020-06-01.
        * e.g. outdir:       hCoV-19__Brazil__PE-FIOCRUZ-IAM1234__2020--2020-06-01.results.
        * e.g. consensus:    hCoV-19/Brazil/PE-FIOCRUZ-IAM1234/2020|2020-06-01.
    * <NUM_THREADS>     -   Number of threads.
    * <DEPTH>           -   Minimum depth to mask unanssembled regions.
    * <MIN_LEN>         -   Minimum length to trimm sequences.
    * <ADAPTERS_FILE>   -   Fasta file with adapters used in the sequencing analysis.

**Suggestion to paired-end reads with 150 of length:**

.. code:: bash
    
    bash sars2_assembly reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa

**Suggestion to paired-end reads with 75 of length:**

.. code:: bash

    bash sars2_assembly reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 35 adapters.fa

Both of those examples will generate the following results:


.. code-block:: text


    current_directory/
    ├-sars2_assembly
    ├-reference.fasta
    ├-code_R1.fastq.gz
    ├-code_R2.fastq.gz
    ├-adapters.fasta
    ├-python_scripts/
    └-prefix_name.results/
     ├-chromosomes.report                            ### tsv file with genomic metrics
     ├-coverage.report                               ### txt file with all assembly metrics
     ├-prefix_name.<R1/R2>.fq.gz                     ### trimmed fastq files
     ├-prefix_name.depthX.fa                         ### consensus defined with iVar
     ├-prefix_name.depthX.amb.fa                     ### consensus defined with iVar with ambiguous nucleotideos on positions where major allele frequencies correspond at least 60% of depth.
     ├-prefix_name.depthX.all.fa                     ### in case of minor variant detection, this file contain the 2 genome versions (major and minor consensus)
     ├-prefix_name.depthX.fa.nextclade.csv           ### or prefix_name.depthX.all.fa.nextclade.csv in case of minor variant detection, nextclade csv output
     ├-prefix_name.depthX.fa.pango.csv               ### or prefix_name.depthX.all.fa.pango.csv in case of minor variant detection, pangolin lineages information
     ├-prefix_name.depthX.fa.bc                      ### bamreadcount output, with all nucleotide frequencies by genomic position
     ├-prefix_name.depthX.fa.bc.intrahost.tsv        ### tsv file with minor variant informations
     ├-prefix_name.depthX.fa.bc.intrahost.short.tsv  ### short tsv file with minor variant informations
     ├-prefix_name.depthX.minor.fa                   ### fasta file with minor consensus genome
     ├-prefix_name.quality.html                      ### html file with quality controll informations
     ├-prefix_name.sorted.bam                        ### sorted bam file
     ├-prefix_name.sorted.bam.bai                    ### index of sorted bam file
     ├-prefix_name.time.txt                          ### time in minutes of each step of analysis.
     └-prefix_name.tsv                               ### tsv output from iVar with the frequencies of iSNVs

=====
Disclaimer
=====
* The fastq files should be in the same directory of sars2_assembly and the python scripts.
* The minor consensus version is based only on replacing the nucleotide from the consensus (majority consensus) with the minor allele (supported by 5 to 49% of the reads), without any statistical method to reconstruct quasispecies genomic populations. For minor variants with percentage near of 50%, the results of this step should be curated mannualy owing the possibility of different frequencies from ivar and bamreadcount analysis.
* In the interactive container with Docker, a pangolin update is strongly recommended (pangolin --update);
* Using Dockerfile or Singularity a pangolin update will be performed automatically, but periodical updates are recommended (re-building the docker image);
* If you use this workflow for academic  purposes, please cite this repository;
* More information `Here <https://dezordi.github.io/>`_;