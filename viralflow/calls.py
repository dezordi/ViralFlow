import subprocess
import shlex


def run_fastp(fastq1, fastq2, prefixout, threads, adapters, min_len, trimm):
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
    # mount command line
    input_opts = f'-i {fastq1} -I {fastq2} --adapter_fasta {adapters} '
    threads_opt = f'--thread {threads} '
    prefix_opt = f'-o {prefixout}.R1.fq.gz -O {prefixout}.R2.fq.gz -h {prefixout}.quality.html '
    var_opt = f'-l {min_len} -f {trimm} -t {trimm} -F {trimm} -T {trimm} '
    cte_opt = '--cut_front --cut_tail --qualified_quality_phred 20 '
    fastp = ("fastp  "+input_opts+prefix_opt+cte_opt+var_opt+threads_opt)
    fastp = shlex.split(fastp)

    # run command
    cmd_fastp = subprocess.Popen(fastp, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
    cmd_fastp.wait()

    # write log file
    stdout, stderr = cmd_fastp.communicate()
    out_f = open(prefixout+'_fastp_out.log', 'wb')
    out_f.write(stdout)
    out_f.close()

    err_f = open(prefixout+'_fastp_err.log', 'wb')
    err_f.write(stderr)
    out_f.close()

    if cmd_fastp.returncode != 0:
        raise Exception(f'Invalid result: { cmd_fastp.returncode }')
        m = "Error in fastp, check log file (" + prefixout+'_fastp_err.log)'
        raise Exception(m)
