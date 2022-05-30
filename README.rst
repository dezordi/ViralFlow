ViralFlow
=========

This repository contains the ViralFlow code, a workflow that performs reference guided genome assembly of SARS-CoV-2 and a number of additional analyses. Python scripts were developed based on the Wuhan SARS-CoV-2 reference genome NC_045512.2. The workflow was developed to work with Illumina paired-end reads. Other sequencing technologies were not accessed so far and any attempt to run ViralFlow with reads from other platforms should be carefully evaluated. The ViralFlow authors can not provide support for other sequencing platforms for now.

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
* pangolin v4.0.5
* bedtools v2.27.1
* bamdst 1.0.9
* seqtk 1.3-r106

To container modes:

* singularity-ce version ≥ 3.8.0
* Docker version ≥ 20.10.12

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

=====
How to install Singularity to run ViralFlow using singularity container
=====

You can check the complete documentation `Here <https://sylabs.io/guides/3.8/admin-guide/installation.html/>`_.

.. code-block:: text
   
   #install dependencies
   
   sudo apt-get update && sudo apt-get install -y \
    build-essential \
    uuid-dev \
    libgpgme-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup-bin

   #install GO
   
   export VERSION=1.14.12 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz
    
   echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc
    
   #download singularity
   
   git clone https://github.com/sylabs/singularity.git && \
    cd SingularityCE && \
    git checkout v3.8.0
   
   #compile singularity
   
   ./mconfig && \
    make -C ./builddir && \
    sudo make -C ./builddir install

=====
How to install Docker to run ViralFlow using docker container
=====

You can check the complete documentation `Here <https://docs.docker.com/engine/install/ubuntu/>`_.

.. code-block:: text
   
   #install
   
   curl -fsSL https://get.docker.com | bash 
   sudo usermod -aG docker <your_username>
   newgrp docker
   docker version
   docker container ls
   systemctl enable docker
   
   #check installation
   
   docker container ls
   docker container run -ti hello-world

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
  viralflow --compileOutput -inputDir ./test_files/ -outDir ./test_files/RESULTS/

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

We recommend the use of the Wuhan SARS-CoV-2 reference genome NC_045512.2, which is included in the ./test_files/reference.fasta.

viralflow --run
-------

This option can be used if the user compiles all the dependencies into the local machine.

viralflow --runContainer
-------

This option can be used if the user has singularity or docker pre installed into the local machine. A file with parameters should be parsed in the argument -inArgisfile (see Quick guide). The following arguments can be parsed:

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

This option compile all the output generated in a batch of samples, resulting in the following directory structure:

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

This option summarize the lineage information.

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

After running, a directory *.results will be created for each sample.

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
Frequently Asked Questions
=====

* Should I build the docker/singularity containers every single run?

Answer: No, ViralFlow --build command should be performed only one time. Obs. in each new run using the once built ViralFlow image pangolin tool will be updated automatically.

* Should I necessarily provide as input the adapters file to ViralFlow?

Answer: Yes, in the current version of ViralFlow the argument of adapters file (primers used in the PCR amplification step) is mandatory. More flexibility of this parameter will be implemented in the next versions of ViralFlow. If you are not working with amplicon sequencing, you can pass an empty file as an argument.

* My ViralFlow run froze after the consensus generation step, why?

Answer: Check if you have directory results (prefix.results) in the output directory. If that is the case bamdst tool will stop and ask if you want to replace the original bamdst outputs, you can digit 2x y (yes) in the terminal that was running the ViralFlow or delete the previous results and re-run pipeline.

* My ViralFlow run reports an error in the pangolin update step, why?

Answer: It can happen for 2 reasons: The first one is related to problems in your local internet, the pangolin update should be performed in an environment with internet access. The second one is owing to possible new versions of the pangolin tool, which depends on new dependencies. Our team normally fixes it in a meantime of 2 days after the new pangolin versions.

* Can I use the ViralFlow for other viruses?

Answer: Yes, but you will need to change the genomic regions on the intrahost_script.py because it is hardcoded with Wuhan SARS-CoV-2 reference genome NC_045512.2 positions. You will need to replace these genomic positions with the annotation of the virus you are interested in. Moreover, you should ignore the pangolin and nextclade outputs since they are tailored for SARS-CoV-2 for now. We are working on a more flexible version of the workflow where the users will be able to more easily work with amplicon based sequencing of other viruses.

=====
Disclaimer
=====
* The adapters and reference file should be in the same directory of fastq files.
* The minor consensus version is based only on replacing the nucleotide from the consensus (majority consensus) with the minor allele (supported by 5 to 49% of the reads), ViralFlow is not based on statistical methods to reconstruct quasispecies genomic populations. For minor variants with a percentage near 50%, the results of this step should be curated manually owing to the possibility of different frequencies from ivar and bamreadcount analysis.
* More information `Here <https://dezordi.github.io/>`_;
