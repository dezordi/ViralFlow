import os
import subprocess
__author__ = "Antonio Marinho da Silva Neto"
__copyright__ = "Copyright 2021, Rede Genomica Fiocruz"
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Antonio Marinho da Silva Neto"
__email__ = "antonio.marinho@fiocruz.br"
__status__ = "Prototype"
__username__ = 'AMarinhoSN'

'''
Functions to build and running containers are defined and provided
'''


def buildSing(output_dir, singfl_path, container_name='viralflow_container',
              sing_path='/usr/local/bin/singularity',
              sing_opt="--fakeroot --sandbox"):
    '''
    build containers singularity container

    Parameters
    ----------
    output_dir:<path>
        directory path to write the container
    singfl_path:<path>
        path to singularity file
    container_name:<str>
        name of the container (default = 'viralflow_container')
    sing_path:<path>
        singularity binary path (default = '/usr/local/bin/singularity')
    sing_opt:<str>
        singularity build options (default = '--fakeroot --sandbox')
    '''
    # ---- sanity check -------------------------------------------------------
    # check if singularity is installed
    try:
        assert(os.path.exists(sing_path))
    except(AssertionError):
        msg_1 = sing_path + ' does not exist, be sure singularity is available'
        msg_2 = ' and the path provided is correct'
        raise Exception(msg_1+msg_2)

    # check if singularity file path is correct
    try:
        assert(os.path.exists(singfl_path))
    except(AssertionError):
        msg_1 = singfl_path+' does not exist. Be sure a valid Singularityfile '
        msg_2 = ' is provided.'
        raise Exception(msg_1+msg_2)

    # check output_dir
    try:
        assert(os.path.isdir(output_dir))
    except(AssertionError):
        raise Exception(output_dir+' is not a valid directory.')

    # -------------------------------------------------------------------------
    # run command
    cmd_str = sing_path+' build '+sing_opt+' ' + output_dir+container_name
    cmd_str += ' '+singfl_path
    os.system(cmd_str)


def run_sing_container(container_img, inArgsFile, input_dir, ref_gnm,
                       adapters_file, threads=1, depth=5, min_len=75,
                       min_dp_intrahost=100, trim_len=0, cpus_pprc=0,
                       sing_call='singularity', dry=False):
    '''
    run singularity viralflow singularity container on a pair of fastq files

    Parameters
    ----------
    container_img:<path>
        Path to viralflow singularity container

    input_dir:<path>
        Path for directory containing input files

    ref_gnm:<path>
        Path for reference genome (must be relative to input_dir)

    fastq_R1:<path>
        Path for fastq R1 filename (must be relative to input_dir)

    fastq_R2:<path>
        Path for fastq R2 filename (must be relative to input_dir)

    prefix_out:<str>
        Prefix for output filenames (should not have '/')

    adapters_file:<path>
        Path for adapters (primers) fasta file (must be relative to input dir)

    threads:<int>
        Number of threads to use on the sample process (default=1)

    depth:<int>
        Minimum depth to mask unanssembled regions

    min_len:<int>
        Minimum length to trimm sequences

    min_dp_intrahost:<int>
        Minimum depth value to consider intrahost minor allele (default = 100)

    trim_len:<int>
        Length to trimm front and tail of reads on fastp analysis (default = 0)

    sing_call:<str>
        string to call singularity, must be absolute path or 'singularity'
        (default = 'singularity')
    '''
    # assumes all files are at input dir
    argfl = inArgsFile.split('/')[-1]
    # get command line
    cmd = sing_call+' run --bind '
    cmd += input_dir+':/data/ '
    #cmd += '--env ARGS_FILE='+argfl+' '
    cmd += '--env FASTA='+argfl+' '
    #cmd += '--env FASTQ1='+fastq_R1+' '
    #cmd += '--env FASTQ2='+fastq_R2+' '
    #cmd += '--env PREFIXOUT='+prefix_out+' '
    #cmd += '--env THREADS='+str(threads)+' '
    #cmd += '--env DEPTH='+str(depth)+' '
    #cmd += '--env MIN_LEN='+str(min_len)+' '
    #cmd += '--env ADAPTERS='+adapters_file+' '
    #cmd += '--env CPUS_PSMPL='+str(cpus_pprc)+' '
    #cmd += '--env DP_INTRAHOST='+str(min_dp_intrahost)+' '
    #cmd += '--env TRIMM_LEN='+str(trim_len)+' '
    cmd += '--writable-tmpfs '+container_img
    # get write output
    #cmd += ' > '+input_dir+prefix_out+'.log'
    # run command
    if dry is True:
        print(cmd)
        return None
    #os.system(cmd)
    subprocess.run(cmd, shell=True, check=True)
