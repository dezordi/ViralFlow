import subprocess
import shlex
import os
from pathlib import Path


def __run_command(cmd_str):
    '''
    run command line string

    Parameters
    ----------
    cmd_str: <str>
        a string of a given command line

    Returns
    -------
    subprocess.CompletedProcess instance
    '''
    # run command
    cmd_str = shlex.split(cmd_str)
    cmd = subprocess.Popen(cmd_str, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    cmd.wait()

    return cmd


def __write_popen_logs(proc, prefixout):
    '''
    write subprocess output to log files

    Parameters
    ---------
    proc:<subprocess.CompletedProcess instance>
        subprocess object

    prefixout:<str>
        prefix to be used on output files

    '''
    stdout, stderr = proc.communicate()
    out_f = open(prefixout+'_out.log', 'wb')
    out_f.write(stdout)
    out_f.close()

    err_f = open(prefixout+'_err.log', 'wb')
    err_f.write(stderr)
    out_f.close()


def __check_status(proc, log_file, proc_label):
    '''
    raise error if process return code is not 0

    Parameters
    ----------
    proc:<subprocess.CompletedProcess instance>
        subprocess object

    log_file:<str>
        path for log file to show at error message

    proc_label:<str>
        label for the command to show at error message

    '''

    if proc.returncode != 0:
        l_nm = log_file
        m = "Error in "+proc_label+", check log file (" + l_nm + ')'
        print(m)
        exit(1)
        #raise Exception(f'Invalid result: { proc.returncode }')


def run_fastp(fastq1, fastq2, prefixout, threads, adapters, min_len, trimm,
              out_dir):
    '''
    Perform bwa index analysis using fastp.

    Parameters
    ----------
    fastq1:<path>
        path for fastq R1 file

    fastq2:<path>
        path for fastq R2 file

    prefix_out:<str>
        prefix for ouput files

    threads:<int>
        number of cpu threads to run fastp

    adatpers:<path>
        valid path for a fasta file containing adapters to process

    trim:<int>
        Length to trimm front and tail of reads on fastp analysis

    '''
    # --- Sanity check --------------------------------------------------------
    if out_dir.endswith('/') is False:
        out_dir += '/'

    assert(Path(out_dir).is_dir()), out_dir+" does not exist"

    # -------------------------------------------------------------------------
    # mount command line
    prfx_wdir = out_dir+prefixout
    input_opts = f'-i {fastq1} -I {fastq2} --adapter_fasta {adapters} '
    threads_opt = f'--thread {threads} '
    prefix_opt = f'-o {prfx_wdir}.R1.fq.gz -O {prfx_wdir}.R2.fq.gz -h {prfx_wdir}.quality.html '
    var_opt = f'-l {min_len} -f {trimm} -t {trimm} -F {trimm} -T {trimm} '
    cte_opt = '--cut_front --cut_tail --qualified_quality_phred 20 '
    fastp = ("fastp  "+input_opts+prefix_opt+cte_opt+var_opt+threads_opt)
    #fastp = shlex.split(fastp)

    # run command
    cmd_fastp = __run_command(fastp)
    # write log file
    __write_popen_logs(cmd_fastp, prfx_wdir+'_fastp')
    # check if errors
    __check_status(cmd_fastp, prfx_wdir+'_fastp_err.log', 'FASTP')


def align_reads2ref(reference_fasta, R1_fq, R2_fq, prefixout, threads, out_dir):
    '''
    Align reads to reference sequence using Burrows-Wheeler (BWA-MEM) Aligner
    algorithm. In addition, sort and index BAM files generated, via SAMTOOLS.

    Parameters
    ----------
    reference_fasta:<path>
        path for fasta file containing the reference genome

    R1_fq:<path>
        path for R1 sense index file (usually '<something>.R1.fq.gz')

    R2_fq:<path>
        path for R2 sense index file (usually '<something>.R2.fq.gz')

    prefix_out:<str>
        prefix for output files

    threads:<int>
        number of cpus threads to use
    '''
    # --- Sanity check --------------------------------------------------------
    if out_dir.endswith('/') is False:
        out_dir += '/'

    assert(Path(out_dir).is_dir()), out_dir+" does not exist"
    # -------------------------------------------------------------------------
    # get output naming pattern
    prfx_wdir = out_dir+prefixout

    # run aligner -------------------------------------------------------------
    t_opt = f'-t {threads} '
    files = f"{reference_fasta} {R1_fq} {R2_fq} "
    output = f'-o {out_dir+prefixout}.bam'
    bwa_mem = ("bwa mem "+t_opt+files+output)

    p_bwa_mem = __run_command(bwa_mem)
    __write_popen_logs(p_bwa_mem, prfx_wdir+'_bwa_mem')
    # check if there was some error
    __check_status(p_bwa_mem, prfx_wdir+'_bwa_mem_err.log', 'BWA MEM')

    # -------------------------------------------------------------------------

    # Sort alignments by leftmost coordinates ---------------------------------
    sam_sort = (f"samtools sort -o {prfx_wdir}.sorted.bam {prfx_wdir}.bam")
    p_sam_sort = __run_command(sam_sort)
    __write_popen_logs(p_sam_sort, prfx_wdir+'_sam_sort')
    __check_status(p_sam_sort, prfx_wdir+'_sam_sort_err.log', 'SAM SORT')

    # -------------------------------------------------------------------------

    # Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast
    # random access.
    sam_index = (f"samtools index {prfx_wdir}.sorted.bam")
    p_sam_index = __run_command(sam_index)
    __write_popen_logs(p_sam_index, prfx_wdir+'_sam_index')
    __check_status(p_sam_index, prfx_wdir+'_sam_index_err.log', 'SAM INDEX')
    # remove bam files
    os.remove(prfx_wdir+'.bam')
