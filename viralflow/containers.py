import os

'''
Here methods to handles build and running containers are defined and provided
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
