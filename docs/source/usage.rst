Usage
=====

.. _installation_mac:

Installation (MacOS)
--------------------

To use ViralFlow, first install it. Due to the limitation of using singularity on MacOS, to run ViralFlow on this type of system, we suggest using a Linux virtualization software called Lima. In this way, the user must follow three steps for installing Lima, and then satisfied with the Ubuntu installation guide.

.. code-block:: console

   # install lima
   brew install lima
   # create an ubuntu instance
   limactl start
   #start ubuntu
   lima

After that, follow the steps of ubuntu install :ref:`installation_ubuntu`

.. _installation_ubuntu:

Installation (Ubuntu)
---------------------

To install ViralFlow, four steps are necessary: Install system dependencies, in case you haven't installed them; Install Conda; install ViralFlow and, finally, assemble the containers for the analyses. This process is performed only once.

ViralFlow was developed and tested for the following operational systems and versions:

* Ubuntu 20.04 LTS;
* Ubuntu 22.04 LTS.

Installing system dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you don't have the pip dependency installer, the git version control system, and the uidmap package, you must install it with the following lines:

.. code-block:: console

   sudo apt update
   sudo apt upgrade
   sudo apt install git
   sudo apt install python3-pip
   sudo apt-get install uidmap

Installing and configuring the conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you don't have the conda environment manager installed, you can install and configure it with the following lines:

.. code-block:: console

   cd $HOME
   wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
   bash Miniconda3-latest-Linux-x86_64.sh -b
   $HOME/miniconda3/bin/conda init
   source $HOME/.bashrc
   conda update conda -y


Installing ViralFlow
~~~~~~~~~~~~~~~~~~~~

If you already have the aforementioned dependencies and conda installed, you can install ViralFlow with 5 lines of code:

.. code-block:: console

   git clone https://github.com/dezordi/ViralFlow.git
   cd ViralFlow/
   conda env create -f envs/env.yml
   conda activate viralflow
   pip install -e .

Building the containers
~~~~~~~~~~~~~~~~~~~~~~~

All steps of ViralFlow are performed in controlled environments, where each tool will have the same version regardless of the research group that is using the tool. ViralFlow has its own method to carry out all this construction of environments, running just one line of code.

.. code-block:: console

   viralflow -build_containers

.. _running:

Running
-------

ViralFlow provides 2 usage modes: sars-cov2 and custom. Regardless of the mode, the user must provide the absolute paths (the entire path to the directory or file to be indicated for the pipeline eg /home/user/test/) for each input file in the command line, otherwise the pipeline will be interrupted during execution

Customizing the snpEff database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the snpEff tool in ViralFlow is configured with the NC_045512.2 genome of the SARS-CoV-2 virus only. If you want to include the snpEff analysis for other viruses, you must update the snpEff database with the following line, example with Dengue:

.. code-block:: console

   viralflow -add_entry_to_snpeff --org_name Dengue --genome_code NC_001474.2

SARS-CoV-2
~~~~~~~~~~

In this model, the analysis is performed based on the reference genome NC_045512.2, and has, as additional analysis, the signature of strains with the Pangolin tool, and signing clades and mutations with the Nextclade tool. For this, the user needs to build a file with the analysis parameters,  `an example can be seen here <https://viralflow.github.io/index-en.html#:~:text=In%20this%20model,seen%20here.>`_.

.. code-block:: console

   viralflow -run --params_file test_files/sars-cov-2.params


Custom
~~~~~~

In this model, the analysis is performed based on the files for the virus that the user wants to analyze. In this mode, the user is responsible for providing each of the files necessary for the analysis. If the user wants to perform the snpEff analysis, he must pass the refseq code viral genome.

.. code-block:: console

   viralflow -run --params_file test_files/denv.params

Pangolin update
~~~~~~~~~~~~~~~

Periodically the pangolin tool updates the lineage database, as well as the usher classification phylogeny, the scorpion mutation constellations, and the pangoLearn trained model. To update the tool and/or it's databases, just run ViralFlow with one of the commands:

.. code-block:: console

   #update the tool and databases
   viralflow -update_pangolin

   #update only the tool
   viralflow -update_pangolin_data
