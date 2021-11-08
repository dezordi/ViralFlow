ViralFlow
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
* pangolin v3.1.11
* bedtools v2.27.0
* bamdst 1.0.6
* seqtk 1.3-r106

=====
How to install
=====

You can install viralflow via pip

.. code-block:: text

  git clone https://github.com/dezordi/ViralFlow.git
  cd ViralFlow/
  git checkout new_cli
  conda env create -f envs/viralflow.yml
  conda activate viralflow
  pip install -e ./

The recommended way to run ViralFlow is via **Singularity container**, be sure `Singularity is installed <https://hub.docker.com/repository/docker/dezordi/iam_sarscov2/>`_.
If you plan to run on your local environment, be sure all requirements are met.

=====
Files info
=====

.. code-block:: text

    IAM_SARSCOV2/
    ├-Singularityfile                       ### Recipe to build local singularity sandbox
    ├-sars2_assembly_singularity            ### Script called into ENTRYPOINT of local singularity sandbox
    ├-sars2_assembly_singularity_run.sh     ### Script for users unfamiliar with singularity run sintaxe
    ├-pango_update                          ### Script to activate conda and update pangolin, run automatically during docker or singularity build
    ├-setup.py                              ### install instructions for pip
    ├-viralflow
    | ├-__init__.py                         ### viralflow python library definition
    | ├-calls.py                            ### command calls module
    | ├-containers.py                       ### containers handling functions module
    | ├-intrahost.py                        ### intrahost bam processing functions module
    | └-pipeline.py                         ### wrapper functions for running pipeline
    |
    ├-scripts:
    | └-viralflow                           ### CLI ViralFlow interface
    └-images:
      └-workflow.png                        ### image of workflow

====
Quick guide
====

Building and running a ViralFlow singularity container

.. code:: bash

  viralflow --build -singFilePath /path/to/ViralFlow/Singularityfile_test
  viralflow --runContainer -inputDir path/to/input/  \
                           -referenceGenome reference_genome.fasta \
                           -adaptersFile adapters.fasta -totalCpus 4 \
                           -depth 5 -minLen 75 \
                           -containerImg /path/to/viralflow_container \
                           -minDpIntrahost 100 -trimLen 0

Run locally (Be sure all requirements are met on your machine)

.. code:: bash

  viralflow --run -inputDir path/to/input/data/ \
                  -referenceGenome $FASTA \
                  -adaptersFile adapters.fasta -totalCpus 4 -depth 5 \
                  -minLen 75 -minDpIntrahost 100 -trimLen 75 \
                  -nxtBin /path/to/nextclade \
                  -nxtDtset /path/to/nextclade/dataset/sars-cov-2/ -v


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

    singularity build --fakeroot viralflow.0.0.5 Singularityfile
    bash sars2_assembly_singularity_run.sh /my/input/dir/ reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa viralflow:0.0.5

For Singularity > 3.7.1 versions, follow:

.. code:: bash

    singularity build --fakeroot --sandbox <imagename> Singularityfile
    bash sars2_assembly_singularity_run.sh <PATH_TO_INPUT_DIR> <REFERENCEGENOME> <001.fastq.gz> <002.fastq.gz> <PREFIX> <NUM_THREADS> <DEPTH> <MIN_LEN> <ADAPTERS_FILE> <imagename>

This method will create a sandbox, and all files to analysis should be in the same directory of the sandbox.
The input directory will be mounted on the container directory /data/ and the ViralFlow repository will be available inside de container at /app/

=====
Explained Usage
=====

**Into interactive docker container**

.. code:: bash

    bash sars2_assembly <REFERENCEGENOME> <001.fastq.gz> <002.fastq.gz> <PREFIX> <NUM_THREADS> <DEPTH> <MIN_LEN> <ADAPTERS_FILE>

* Arguments
    * <REFERENCEGENOME> -   Fasta file with reference genome.
    * <001.fastq.gz>    -   Fasqt file with positive sense reads (R1).
    * <002.fastq.gz>    -   Fastq file with negative sense reads (R2).
    * <PREFIX>          -   Prefix string to store results and to rename consensus genome. The user can set the gisaid format genome name, and the workflow will automatically format the consensus name, as the prefix will be used to create the directory output, the slash '/' should be replaced by '__' and the pipe '|' should be replaced by '--'.
        * e.g. prefix:       hCoV-19__Brazil__PE-FIOCRUZ-IAM1234__2020--2020-06-01.
        * e.g. outdir:       hCoV-19__Brazil__PE-FIOCRUZ-IAM1234__2020--2020-06-01.results.
        * e.g. cons.:    hCoV-19/Brazil/PE-FIOCRUZ-IAM1234/2020|2020-06-01.
    * <NUM_THREADS>     -   Number of threads.
    * <DEPTH>           -   Minimum depth to mask unanssembled regions.
    * <MIN_LEN>         -   Minimum length to trimm sequences.
    * <ADAPTERS_FILE>   -   Fasta file with adapters used in the sequencing analysis.
    * <DP_INTRAHOST>    -   Argument created on workflow v.0.0.5. Minimum depth value to consider intrahost minor allele, optional, default = 100.
    * <TRIMM_LEN>       -   Argument created on workflow v.0.0.5. Length to trimm front and tail of reads on fastp analysis,optional, default = 0.

**Suggestion to paired-end reads with 150 of length:**

.. code:: bash

    bash sars2_assembly reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa

**Suggestion to paired-end reads with 150 of length, considering 50 of depth threshold for intrahost minor alleles:**

.. code:: bash

    bash sars2_assembly reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa 50

**Suggestion to paired-end reads with 150 of length, considering 50 of depth threshold for intrahost minor alleles and trimming 10 bases of front and tail of reads:**

.. code:: bash

    bash sars2_assembly reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa 50 10

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
     ├-prefix_name.depthX.fa.algn.minor.fa           ### fasta file with minor consensus genome
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
