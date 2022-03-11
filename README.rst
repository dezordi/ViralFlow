ViralFlow
=========

This repository contains the code  of Viralflow, a workflow to performs a reference guided genome assembly of SARS-CoV-2. Python scripts were developed based on the wuhan SARS-CoV-2 reference genome NC_045512.2. The workflow was developed to work with Illumina paired-end reads. Tests with other technologies should be performed.

If you use this workflow for academic purposes, please cite: `ViralFlow: A Versatile Automated Workflow for SARS-CoV-2 Genome Assembly, Lineage Assignment, Mutations and Intrahost Variant Detection <https://www.mdpi.com/1999-4915/14/2/217>`_.

.. image:: images/workflow_develop.png
   :width: 600
   :align: center

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
* mafft v7.453 (2019/Nov/8)
* nextclade 1.4.5
* pangolin v3.1.11
* bedtools v2.27.1
* bamdst 1.0.9
* seqtk 1.3-r106

=====
How to install
=====

You can install viralflow via pip

.. code-block:: text

  git clone https://github.com/dezordi/ViralFlow.git
  cd ViralFlow/
  conda env create -f envs/viralflow.yml
  conda activate viralflow
  pip install -e ./

The recommended way to run ViralFlow is via **Singularity container**, be sure `Singularity is installed <https://hub.docker.com/repository/docker/dezordi/iam_sarscov2/>`_. But you can also run via Docker and with local environment. If you plan to run on your local environment, be sure all requirements are met.

====
Quick guide
====

Run locally (Be sure all requirements are met on your machine)

.. code:: bash

  viralflow --run -inputDir path/to/input/data/ \
                  -referenceGenome $FASTA \
                  -adaptersFile adapters.fasta -totalCpus 4 -depth 5 \
                  -minLen 75 -minDpIntrahost 100 -trimLen 0 \
                  -nxtBin /path/to/nextclade \
                  -nxtDtset /path/to/nextclade/dataset/sars-cov-2/ -v

Building and running ViralFlow with singularity container

.. code:: bash

  viralflow --build -singFilePath ./Singularityfile
  viralflow --runContainer -inArgsFile ./test_files/test_args.conf

Building and running ViralFlow with docker container

.. code:: bash

  viralflow --build -containerService docker
  viralflow --runContainer -containerService docker -inArgsFile ./test_files/test_args_docker.conf


Compile the outputs

.. code:: bash

  viralflow --compileOutput -inputDir <path/to/directory/with/results> -outDir <path/to/store/compiled/results>
  #example
  viralflow --compileOutput -inputDir ./test_files/ -outDir ./test_files/

Check negative controls

.. code:: bash

  viralflow --checkNegControls -negControlLabels <negative_control_sample_code> -pangoCSV <path/to/compiled/pango.csv>
  #example
  viralflow --checkNegControls -negControlLabels Cneg_R1 -pangoCSV ./test_files/RESULTS/pango.csv

Get lineage summary

.. code:: bash

  viralflow --getLineageSummary -pangoCSV <path/to/compiled/pango.csv> -chromCSV <path/to/compiled/chromossomes.csv> -outDir <path/to/store/summaries>
  #example
  viralflow --getLineageSummary -pangoCSV ./test_files/RESULTS/pango.csv -chromCSV ./test_files/RESULTS/chromossomes.csv -multifasta ./test_files/RESULTS/seqbatch.fa -outDir ./test_files/RESULTS/

=====
Explained Usage
=====

Reference Genome
-------

We recommend the use of wuhan SARS-CoV-2 reference genome NC_045512.2, which is included in the ./test_files/reference.fasta.

viralflow --run
-------

This option can be used if the user compile all the dependencies into the local machine.

viralflow --runContainer
-------

This option can be used if the user have singularity or docker pre installed into the local machine. A file with paremeters should be parsed in the argument -inArgisfile (see Quick guide). The following arguments can be parsed:

.. code-block:: text

  inputDir                   ### The path to directory with fastq.gz files -- e.g. ./test_files;
  referenceGenome            ### The name of fasta file with reference genome, this file should be inside the directory with fastq.gz files (inputDir) -- e.g. reference.fasta;
  adaptersFile               ### The name of fasta file with adapters, this file should be inside the directory with fastq.gz files (inputDir) -- e.g. ART_adapters.fa;
  totalCpus                  ### Total CPU number used in workflow -- e.g. 2;
  depth 5                    ### Minimum depth to consider a sequenced nucleotide -- e.g. 5;
  minLen 75                  ### Minimum length to maintain reads after fastp processing -- e.g. 75;
  containerImg               ### Container image -- e.g. ./viralflow_container with Singularity or viralflow_container:latest with Docker; 
  minDpIntrahost 100         ### Minimum depth to consider an iSNV -- e.g. 100;
  trimLen 0                  ### Length to trim read extremities -- e.g. 0;
  cpusPerSample 1            ### Number of CPUs per sample during analysis.

viralflow --compileOutput
-------

This option compile all the outpus generated in a batch of samples, resulting in the following directory structure:

.. code-block:: text

    inputDir/
    ├-reference.fasta
    ├-code1_R1.fastq.gz
    ├-code1_R2.fastq.gz
    ├-code2_R1.fastq.gz
    ├-code2_R2.fastq.gz
    ├-adapters.fasta
    ├-code1.results/
    ├-code2.results/
    └-RESULTS/
     ├-chromosomes.csv       ### csv file with compiled bamdst information;
     ├-erross_detected.csv   ### csv file with errors detected by sample;
     ├-mutations.csv         ### csv file with nucleotidic mutations by sample;
     ├-nextclade.csv         ### csv file with compiled nextclade information;
     ├-pango.csv             ### csv file with compiled pangolin information;
     └-seqbatch.fa           ### fasta file with consensus sequences;

viralflow --checkNegControls
-------

This option check if some sample have the same linage of negative control.

viralflow --getLineageSummary
-------

This option summarize the lineage information;

.. code-block:: text

    inputDir/
    ├-reference.fasta
    ├-code1_R1.fastq.gz
    ├-code1_R2.fastq.gz
    ├-code2_R1.fastq.gz
    ├-code2_R2.fastq.gz
    ├-adapters.fasta
    ├-code1.results/
    ├-code2.results/
    └-RESULTS/
     ├-chromosomes.csv       ### csv file with compiled bamdst information;
     ├-erross_detected.csv   ### csv file with errors detected by sample;
     ├-mutations.csv         ### csv file with nucleotidic mutations by sample;
     ├-nextclade.csv         ### csv file with compiled nextclade information;
     ├-pango.csv             ### csv file with compiled pangolin information;
     ├-seqbatch.fa           ### fasta file with consensus sequences;
     ├-major_summary.csv     ### csv file with depth, coverage, and lineage of major consensus genomes;
     ├-minor_summary.csv     ### csv file with depth, coverage, and lineage of minor consensus genomes;
     └-lineage_summary.csv   ### csv file lineage count of batch analysis.

=====
Files info
=====

After running, a directory .results will be created for each sample.

Repository directory structure

.. code-block:: text

    ViralFlow/
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

Results directory structure for each sample

.. code-block:: text


    inputDir/
    ├-reference.fasta
    ├-code_R1.fastq.gz
    ├-code_R2.fastq.gz
    ├-adapters.fasta
    └-code.results/
     ├-chromosomes.report                                  ### tsv file with genomic metrics
     ├-coverage.report                                     ### txt file with all assembly metrics
     ├-code_<fastp/mafft/nextclade/pangolin/bwa/sam>.log   ### txt file with log of tool
     ├-code.<R1/R2>.fq.gz                                  ### trimmed fastq files
     ├-code.depthX.fa                                      ### consensus defined with iVar
     ├-code.depthX.amb.fa                                  ### consensus defined with iVar with ambiguous nucleotideos on positions where major allele frequencies correspond at least 60% of depth.
     ├-code.depthX.all.fa                                  ### in case of minor variant detection, this file contain the 2 genome versions (major and minor consensus)
     ├-code.depthX.fa.nextclade.csv                        ### or code.depthX.all.fa.nextclade.csv in case of minor variant detection, nextclade csv output
     ├-code.depthX.fa.gene<SC2 genes>.fasta                ### or code.depthX.all.fa.gene<SC2 genes>.fasta in case of minor variant detection, fasta with aminoacid sequence of each gene, generated with nextclade
     ├-code.depthX.fa.pango.csv                            ### or code.depthX.all.fa.pango.csv in case of minor variant detection, pangolin lineages information
     ├-code.depthX.fa.bc                                   ### bamreadcount output, with all nucleotide frequencies by genomic position
     ├-code.depthX.fa.bc.intrahost.tsv                     ### tsv file with minor variant informations
     ├-code.depthX.fa.bc.intrahost.short.tsv               ### short tsv file with minor variant informations
     ├-code.depthX.fa.algn.minor.fa                        ### fasta file with minor consensus genome
     ├-code.fastp.html                                     ### html file with fastp quality controll informations
     ├-code.fastp.json                                     ### json file with fastp quality controll informations
     ├-code.sorted.bam                                     ### sorted bam file
     ├-code.sorted.bam.bai                                 ### index of sorted bam file
     ├-code.time.txt                                       ### time in minutes of each step of analysis.
     └-code.tsv                                            ### tsv output from iVar with the frequencies of iSNVs

=====
Disclaimer
=====
* The adapters and reference file should be in the same directory of fastq files.
* The minor consensus version is based only on replacing the nucleotide from the consensus (majority consensus) with the minor allele (supported by 5 to 49% of the reads), without any statistical method to reconstruct quasispecies genomic populations. For minor variants with percentage near of 50%, the results of this step should be curated mannualy owing the possibility of different frequencies from ivar and bamreadcount analysis.
* Using Dockerfile or Singularity a pangolin update will be performed automatically, but periodical updates are recommended (re-building the docker image);
* More information `Here <https://dezordi.github.io/>`_;
