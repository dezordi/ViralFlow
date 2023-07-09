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

   git clone https://github.com/dezordi/ViralFlow.git -b vfnext
   cd ViralFlow/
   conda env create -f envs/env.yml
   conda activate viralflow
   pip install -e .

Building the containers
~~~~~~~~~~~~~~~~~~~~~~~

All steps of ViralFlow are performed in controlled environments, where each tool will have the same version regardless of the research group that is using the tool. ViralFlow has its own method to carry out all this construction of environments, running just one line of code.

viralflow_dev -build_containers