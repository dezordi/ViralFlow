IAM_SARS-CoV-2
=========


This repository contains general scripts to automatize the routine tasks of
my work at WallauLab. For some scripts, I added test files in the
directory /test\_files

=====
Usage
=====


**To format LTR\_finder output:**

.. code:: bash
    
    python LTR\_finder\_to\_csv.py -in

**To assembly mitogenomes using a reference-guided approach with MITObin:**

.. code:: bash

    bash mito.sh <reference_genome> <adapter_path/adapter_file> <trimmomatic_path/trimmomatic_jar> </home/user/mira_path/> <num_threads>
    #example line: 
    bash mito.sh ref_genome1.fasta sracode_1.fastq sracode_2.fastq Trimmomatic-0.36/adapters/adapter_file.fa ~/softwares/Trimmomatic-0.36/trimmomatic-0.36.jar /home/user/softwares/mira_4.0.2/ 8


**To automatize auto_mito.sh**

.. code:: bash

    bash auto_mito.sh <genomes_fastq_adapters.lst> <path/trimmomatic> <path/mira> <num_threads>
    #example line
    bash auto_mito.sh genomes_fastq_adapters.lst ~/softwares/Trimmomatic-0.36/trimmomatic-0.36.jar /home/user/softwares/mira_4.0.2/ 8

**To download fastq files from a file with a number or SRA codes**

.. code:: bash

    bash sra.sh <sra_file_list.txt>

**To execute cd-hit getting a tsv formatted file of clusters**

.. code:: bash

    python cd_hit_formated.py -in <fasta_file> -p -m -cd "cd-hit line" -md
    #example line for nucleotides: 
    python cd_hit_formated.py -in sequences.fasta -p 16 -m 60000 -cd "-c 0.8 -aL 0.8 -g 1 -d 200"
    #example line for proteins:
    python cd_hit_formated.py -in sequences.fasta -p 16 -m 60000 -cd "-c 0.8 -aL 0.8 -g 1 -d 200" -md prot

==========
Disclaimer
==========
- I'm not a computer engineer or some related profession, I'm just writing this script to study python and to automatize some bioinformatics tasks. So feel free to commit changes that make the code more efficient or cleaner.

- This repository will be constantly updated.