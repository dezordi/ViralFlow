ViralFlow
=========

This repository contains the code of Viralflow, a workflow to performs a reference guided genome assembly of SARS-CoV-2 written in Nextflow. The workflow was developed to work with Illumina paired-end reads. Tests with other technologies should be performed.

If you use this workflow for academic purposes, please cite: `ViralFlow: A Versatile Automated Workflow for SARS-CoV-2 Genome Assembly, Lineage Assignment, Mutations and Intrahost Variant Detection <https://www.mdpi.com/1999-4915/14/2/217>`_.


=====
How to install (quick and dirty)
=====

Viralflow is a nextflow pipeline and the recommended usage is as such. A CLI wrapper was implemented to make it more accessible for users not familiar with nextflow.

You can install viralflow wrapper via pip

.. code:: bash

  git clone https://github.com/dezordi/ViralFlow.git
  cd ViralFlow/
  pip install -e ./


Install conda, singularity and nextflow:

.. code:: bash

  viralflow_dev -setup_dependencies


Build containers

.. code:: bash

  viralflow_dev -build_containers


=====
How to run sars-cov-2 (quick and dirty)
=====

.. code:: bash

  viralflow_dev -run --params_file test_files/sars-cov-2.params

=====
How to run denv (custom) (quick and dirty)
=====

.. code:: bash

  viralflow_dev -run --params_file test_files/denv.params

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

====
Quick guide
====

=====
Explained Usage
=====


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
