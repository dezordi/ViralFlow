from pathlib import Path
import time
import os
import viralflow.calls

'''
Here we define some handy functions to run individual steps or the whole
pipeline. The functions defined here are essentially wrappers for functions
defined at 'viralflow.calls' added of some sanity check for input files.
'''


def run_step_1(fastq1, fastq2, adapters, prefixout, threads, min_len,
               trim, outdir):
    '''
    Run step 1 of ViralFlow, check github README and/or paper for further
    details
    '''
    print('>--- STEP 1 ---<')
    # --- sanity check of arguments ----------------------------------------
    try:
        assert(Path(outdir).is_dir())
    except(AssertionError):
        print("WARNING: "+outdir
              + " does not exist. Don't worry it will be created then. ")
        os.system('mkdir ' + outdir)
    # ----------------------------------------------------------------------
    print('@ running fastp...')
    s = time.time()
    viralflow.calls.run_fastp(fastq1, fastq2, prefixout, threads, adapters,
                              min_len, trim, outdir)
    f = time.time()
    print(' > total execution time : ', f - s)
    # TODO check output for possible errors


def run_step_2(ref_gnm, prefixout, threads, outdir, depth):
    '''
    Run step 2 of ViralFlow, check github README and/or paper for further
    details
    '''
    print('>--- STEP 2 ---<')

    print('@ aligning reads to reference genome...')
    fq1 = outdir+prefixout+'.R1.fq.gz'
    fq2 = outdir+prefixout+'.R2.fq.gz'

    s = time.time()
    #viralflow.calls.align_reads2ref(
    #    ref_gnm, fq1, fq2, prefixout, threads, outdir)
    f = time.time()
    print(' > total execution time : ', f - s)

    print('@  Get consensus and minor variant analysis...')

    s = time.time()
    viralflow.calls.run_ivar(ref_gnm, prefixout, outdir, depth)
    f = time.time()
    print(' > total execution time : ', f - s)

    print('@ get minor variants data...')
    viralflow.calls.get_minor_variants(ref_gnm, prefixout, depth, threads,
                                       outdir)
    # TODO check output for possible errors


def run_step_3(prefixout, depth, intrahost_depth, out_dir,
               depth_number=100, per_limit=0.05):
    '''
    Run step 3 of ViralFlow, check github README and/or paper for further
    details
    '''
    # --- sanity check of arguments ----------------------------------------
    try:
        assert(Path(out_dir).is_dir())
    except(AssertionError):
        print("WARNING: "+out_dir
              + " does not exist. Don't worry it will be created then. ")
        os.system('mkdir ' + out_dir)
    # ----------------------------------------------------------------------

    print('>--- STEP 3 ---<')

    print("@ minor variants identification...")
    print('   > minor variants only if depth >= ', depth_number)
    print('   > putative minor if Depth Base / Total Base >=', per_limit)
    print('   > INPUT FILES:')
    prfx_wdir = out_dir+prefixout
    bam_rc_file = f"{prfx_wdir}.depth{depth}.fa.bc"
    alignment_file = f"{prfx_wdir}.depth{depth}.fa.algn"
    # sanity check ----------------------------------------------------------
    assert(Path(bam_rc_file).is_file()), bam_rc_file+' does not exist.'
    assert(Path(alignment_file).is_file()), alignment_file+' does not exist.'
    # -----------------------------------------------------------------------
    print('     - ', bam_rc_file)
    print('     - ', alignment_file)
    s = time.time()
    viralflow.intrahost.get_instrahost_tsv(bam_rc_file, alignment_file,
                                           per_limit=per_limit,
                                           depth_value=depth_number)
    f = time.time()
    print(' > total execution time : ', f - s)
    # TODO check output for possible errors


def run_step_4(prefixout, depth, threads, out_dir, ref_gnm, nxtdst_str):
    '''
    Run step 4 of ViralFlow, check github README and/or paper for further
    details

    '''
    print('>--- STEP 4 ---<')

    print("@ Naming variants identified...")
    viralflow.calls.get_variant_naming(prefixout, depth, threads, out_dir,
                                       ref_gnm, nxtdst_str)
    # TODO check output for possible errors


def run_step_5(prefixout, outdir):
    '''
    '''
    print('>--- STEP 5 ---<')
    print('@ computings assembly metrics...')
    viralflow.calls.do_assembly_metrics(prefixout, outdir)
    # TODO check output for possible errors
