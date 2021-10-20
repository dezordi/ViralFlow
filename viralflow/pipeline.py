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


def run_step_2(ref_gnm, fq1, fq2, prefixout, threads, outdir, ):
    '''
    Run step 2 of ViralFlow, check github README and/or paper for further
    details
    '''
    print('>--- STEP 2 ---<')

    print('@ aligning reads to reference genome...')

    s = time.time()
    viralflow.calls.align_reads2ref(
        ref_gnm, fq1, fq2, prefixout, threads, outdir)
    f = time.time()

    print(' > total execution time : ', f - s)
    # TODO run ivar
