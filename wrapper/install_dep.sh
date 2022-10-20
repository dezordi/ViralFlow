# install miniconda
echo "@ Install Miniconda..."
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh -b

echo "@ Config Miniconda..."
$HOME/miniconda3/bin/conda init
exec bash
$HOME/miniconda3/bin/conda update conda -y
$HOME/miniconda3/bin/conda create -n viralflow python=3.9 -y
$HOME/miniconda3/bin/conda activate viralflow

# install singularity 
echo "@ install singularity"
conda install -c "conda-forge/label/main" singularity -y

# install nextflow 
echo "@ install nextflow"
conda install -c bioconda nextflow -y
nextflow self-update
