# install miniconda
echo "@ Install Miniconda..."
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh

# install singularity 
echo "@ install singularity"
conda install -c "conda-forge/label/main" singularity

# install nextflow 
echo "@ install nextflow"
conda install -c bioconda nextflow
nextflow self-update
