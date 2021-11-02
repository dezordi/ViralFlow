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


def run_bwa_idex(input_file, out_dir, prefixout):
    '''Perform bwa index analysis'''
    prfx_wdir = out_dir + prefixout
    bwa_index = (f"bwa index {input_file}")
    p_bwa_idx = __run_command(bwa_index)
    # write log file
    __write_popen_logs(p_bwa_idx, prfx_wdir+'_1_bwa_idx')
    # check if errors
    __check_status(p_bwa_idx, prfx_wdir+'_1_bwa_idx_err.log', 'BWA INDEX')


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

    Returns
    -------
    <dct>
        a dictionary with output files names expected to be obtained
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
    __write_popen_logs(cmd_fastp, prfx_wdir+'_1_fastp')
    # check if errors
    __check_status(cmd_fastp, prfx_wdir+'_1_fastp_err.log', 'FASTP')

    # check if expected output files were obtained
    output_flnm = {'R1': prfx_wdir+'.R1.fq.gz',
                   'R2': prfx_wdir+'.R2.fq.gz',
                   'quality_html': prfx_wdir+'.quality.html'}
    return output_flnm


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
    __check_status(p_bwa_mem, prfx_wdir+'_2_bwa_mem_err.log', 'BWA MEM')

    # -------------------------------------------------------------------------

    # Sort alignments by leftmost coordinates ---------------------------------
    sam_sort = (f"samtools sort -o {prfx_wdir}.sorted.bam {prfx_wdir}.bam")
    p_sam_sort = __run_command(sam_sort)
    __write_popen_logs(p_sam_sort, prfx_wdir+'_sam_sort')
    __check_status(p_sam_sort, prfx_wdir+'_2_sam_sort_err.log', 'SAM SORT')

    # -------------------------------------------------------------------------

    # Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast
    # random access.
    sam_index = (f"samtools index {prfx_wdir}.sorted.bam")
    p_sam_index = __run_command(sam_index)
    __write_popen_logs(p_sam_index, prfx_wdir+'_sam_index')
    __check_status(p_sam_index, prfx_wdir+'_2_sam_index_err.log', 'SAM INDEX')
    # remove bam files
    os.remove(prfx_wdir+'.bam')


def run_ivar(reference_fasta, prefixout, out_dir, depth):
    '''
    run ivar steps
    '''
    # --- Sanity check --------------------------------------------------------
    if out_dir.endswith('/') is False:
        out_dir += '/'

    assert(Path(out_dir).is_dir()), out_dir+" does not exist"
    # -------------------------------------------------------------------------
    # mpileup
    # ivar gets output from samtools mpileup, so mount samtools cmd string
    prfx_wdir = out_dir+prefixout
    mpi_p = f"-d 50000 --reference {reference_fasta} -a "
    out_flnm = f"-B {prfx_wdir}.sorted.bam "
    mpileup_str = "samtools mpileup "+mpi_p+out_flnm

    # IVAR STEP 1 -------------------------------------------------------------
    # Call variants
    ivar_1 = f"| ivar variants -p {prfx_wdir} -q 30 -t 0.05 "
    cmd_ivar_1 = mpileup_str+ivar_1
    #print("COMMAND:")
    #print(cmd_ivar_1)

    # WARNING #################################################################
    # for some reason, using shellx using "|" doesn't work very well. For now
    # os.system() works just fine, but THIS SHOULD BE FIX in the future
    # so we get more log output control
    os.system(cmd_ivar_1)
    #p_ivar_1 = __run_command(cmd_ivar_1)
    #__write_popen_logs(p_ivar_1, prfx_wdir+'_ivar_1')
    #__check_status(p_ivar_1, prfx_wdir+'_ivar_1_err.log', 'IVAR_variants')

    # IVAR STEP 2 -------------------------------------------------------------
    # Mount consensus for a given minimum depth
    ivar_2 = f"| ivar consensus -p {prfx_wdir} -q 30 -t 0 -m {depth} -n N"
    cmd_ivar_2 = mpileup_str+ivar_2
    #print("COMMAND:")
    #print(cmd_ivar_2)
    os.system(cmd_ivar_2)
    #p_ivar_2 = __run_command(cmd_ivar_2)
    #__write_popen_logs(p_ivar_2, prfx_wdir+'_ivar_2')
    #__check_status(p_ivar_2, prfx_wdir+'_ivar_2_err.log', "IVAR_consensus")

    # IVAR STEP 3 -------------------------------------------------------------
    d = f"-m {depth}"
    ivar_3 = f"| ivar consensus -p {prfx_wdir}.ivar060 -q 30 -t 0.60 -n N "+d
    cmd_ivar_3 = mpileup_str+ivar_3
    #print('COMMAND:')
    #print(cmd_ivar_3)
    os.system(cmd_ivar_3)

    #p_ivar_3 = __run_command(cmd_ivar_3)
    #__write_popen_logs(p_ivar_3, prfx_wdir+'_ivar_3')
    #__check_status(p_ivar_3, prfx_wdir+'_ivar_3_err.log', 'IVAR_consencus_60')

    # EDIT FILES --------------------------------------------------------------
    #
    mv_str = f"mv {prfx_wdir}.fa {prfx_wdir}.depth{depth}.fa"
    os.system(mv_str)

    mv_str = f"mv {prfx_wdir}.ivar060.fa {prfx_wdir}.depth{depth}.amb.fa"
    os.system(mv_str)

    '''

    s1_str = f"sed -i -e 's/>.*/>|'{prfx_wdir}'/g' {prfx_wdir}.depth{depth}.fa"
    os.system(s1_str)

    s2_str = f"sed -i -e 's/__/\//g' -e 's/--/|/g' {prfx_wdir}.depth{depth}.fa"
    os.system(s2_str)

    mv_str = f"mv {prfx_wdir}.ivar060.fa {prfx_wdir}.depth{depth}.amb.fa"
    os.system(mv_str)

    s3_str = f"sed -i -e 's/>.*/>|'{prfx_wdir}'/g' {prfx_wdir}.depth{depth}.amb.fa"
    os.system(s3_str)

    s4_str = f"sed -i -e 's/__/\//g' -e 's/--/|/g' {prfx_wdir}.depth{depth}.amb.fa"

    os.system(s4_str)
    '''
    #p_ivar_4 = __run_command(edit_str)
    #__write_popen_logs(p_ivar_4, prfx_wdir+'_ivar_4')
    #__check_status(p_ivar_4, prfx_wdir+'_ivar_4_err.log', 'IVAR_edit_file')


def get_minor_variants(reference_genome, prefixout, depth, threads, outdir):
    '''

    '''
    prfx_wdir = outdir+prefixout
    # BAM READCOUNT -----------------------------------------------------------
    # abm = annotation bed merged
    bamread_output = open(f"{prfx_wdir}.depth{depth}.fa.bc", "wb")
    input_ref = f"-f {reference_genome} "
    bam_out = f"{prfx_wdir}.sorted.bam "
    bwa_read = f"bam-readcount -d 50000 -q 30 -w 0 "+input_ref+bam_out
    # run command and store sdtout on bam read
    #print("COMMAND:")
    #print(bwa_read)
    # WARNING #################################################################
    # For some reason, using the __run_command stuck the pipeline at cmd.wait()
    # for now, os.system get the job done, but this needs to be fixed.
    # ------------------------------------------------------------------------
    #p_bwa_read = __run_command(bwa_read)
    # write bwa readcount stidout on fa.bc
    #stdout, stderr = p_bwa_read.communicate()
    #bamread_output.write(stdout)
    #bamread_output.close()
    # CURRENT
    os.system(bwa_read)

    # MAFFT -------------------------------------------------------------------
    #mafft_output = open(f"{prfx_wdir}.depth{depth}.fa.algn", "w")
    input_fa = f"--add {prfx_wdir}.depth{depth}.fa "
    t_str = f"--thread {threads}"
    out = f" --clustalout {prfx_wdir}.depth{depth}.fa.algn"
    mafft = "mafft --keeplength "+input_fa+t_str + \
        f" {reference_genome}"  # > {prfx_wdir}.depth{depth}.fa.algn"
    #print("COMMAND:")
    #print(mafft)
    cmd_mafft = __run_command(mafft)
    log_prfx = prfx_wdir+'_2_mafft'
    __write_popen_logs(cmd_mafft, log_prfx)
    __check_status(cmd_mafft, log_prfx+"_err.log", "MAFFT")

    # mafft print the output on the terminal
    os.system(f'cp {log_prfx}_out.log {prfx_wdir}.depth{depth}.fa.algn')


def get_variant_naming(prefixout, depth, threads, out_dir, ref_gnm,
                       nxt_dataset, nxt_bin='nextclade'):
    '''
    Define clade and variant for a given SARS-CoV2 sequence.
    The variant and clade definition is based on Nextclade and Pangolin.
    '''
    # standard input names
    prfx_wdir = out_dir+prefixout
    intrahost_tsv = f"{prfx_wdir}.depth{depth}.fa.bc.intrahost.short.tsv"
    consensus_fa = f"{prfx_wdir}.depth{depth}.fa "
    # get standard variant files output prefix
    output_prfx = f"{prfx_wdir}.depth{depth}"

    # check for single or multiple sequences at intrahost ---------------------
    num_lines = sum(1 for line in open(intrahost_tsv, 'r'))

    # if more than a single variant, then compile variants on a single file
    if num_lines > 1:
        minor_var_fa = f"{prfx_wdir}.depth{depth}.fa.algn.minor.fa"
        all_fa = f"{prfx_wdir}.depth{depth}.all.fa "
        output_nxt = f"{prfx_wdir}.depth{depth}.all.nextclade.csv"
        output_pgl = output_prfx+".all.fa.pangolin.csv"
        # compile
        cmd_str = "cat "+consensus_fa + minor_var_fa + f" > "+all_fa
        os.system(cmd_str)
        #p_cat = __run_command(cmd_str)
        #__write_popen_logs(p_cat, prfx_wdir+'_4_cat')
        #__check_status(p_cat, prfx_wdir+'_4_cat', 'COMPILE VARIANTS')

    # single sequence, use consensus fasta as input
    if num_lines == 1:
        input_seqs = consensus_fa
        output_nxt = f"{prfx_wdir}.depth{depth}.nextclade.csv"
        #output_prfx+".fa.nextclade.csv"
        output_pgl = output_prfx+".fa.pangolin.csv"
    # if no sequence, something went wrong on previous step
    if num_lines == 0:
        raise Exception("No sequence at "+intrahost_tsv)
    # -------------------------------------------------------------------------

    # NEXTCLADE ---------------------------------------------------------------
    # get clades for each sequence, via nextclade
    in_str = f" -i {input_seqs}"
    job_str = f" --jobs {threads} "
    ref_str = f" --input-root-seq={ref_gnm}"
    nxtdst_str = f" --input-dataset={nxt_dataset}"
    out_csv = f" --output-csv={prfx_wdir}.depth{depth}.all.nextclade.csv"
    out_str = f" --output-dir={out_dir}"
    # nxt_str = "nextclade -i "+input_seqs + job_str + " -c  "+output_nxt
    nxt_str = nxt_bin + in_str + job_str + ref_str + nxtdst_str + out_csv \
        + out_str
    #print("COMMAND NEXTCLADE:")
    #print(nxt_str)
    #os.system(nxt_str)
    # print(nxt_str)
    p_nxt = __run_command(nxt_str)
    __write_popen_logs(p_nxt, prfx_wdir+'_4_nextclade')
    __check_status(p_nxt, prfx_wdir+'_4_nextclade_err.log', "NEXTCLADE")
    # -------------------------------------------------------------------------

    # PANGOLIN ----------------------------------------------------------------
    # get variant namings, via pangolin
    out_str = " --outfile " + output_pgl
    pgl_str = "pangolin  " + input_seqs + f" -t {threads}" + out_str

    p_pgl = __run_command(pgl_str)
    __write_popen_logs(p_pgl, prfx_wdir+'_4_pangolin')
    __check_status(p_pgl, prfx_wdir+'_4_pangolin', "PANGOLIN")


def do_assembly_metrics(prefixout, outdir):
    '''
    do assembly metrics
    '''
    prfx_wdir = outdir+prefixout
    bam_in = f"{prfx_wdir}.sorted.bam"
    bed_out = f"{prfx_wdir}.sorted.bed"

    # BEDTOOLS BAMTOBED

    cmd_bedtools = f"bedtools bamtobed -i {bam_in} > {bed_out}"
    print('COMMAND:')
    print(cmd_bedtools)
    os.system(cmd_bedtools)

    # SAMTOOLS VIEW + BAMDST
    cmd_samtools = f"samtools view {bam_in} -u | bamdst -p {bed_out} -o {outdir}"
    print(cmd_samtools)
    os.system(cmd_samtools)

    # GUNZIP
    gzp_1 = f"gunzip {outdir}region.tsv.gz"
    __run_command(gzp_1)
    gzp_2 = f"gunzip {outdir}depth.tsv.gz"
    __run_command(gzp_2)
