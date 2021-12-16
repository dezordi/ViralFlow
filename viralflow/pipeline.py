from pathlib import Path
import time
import os
import viralflow.calls

"""
Here we define some handy functions to run individual steps or the whole
pipeline. The functions defined here are essentially wrappers for functions
defined at 'viralflow.calls' added of some sanity check for input files.
"""
script_file = _path = os.path.realpath(__file__)
VIRALFLOW_PATH = "/".join(script_file.split("/")[0:-2]) + "/"


def run_step_0(reference_genome, outdir, verbose=False):
    """
    reference genome index mapping
    """
    if verbose is True:
        print("@ running pangolin update...")
    viralflow.calls.pango_update()
    if verbose is True:
        print("@ running bwa index...")
    s = time.time()
    viralflow.calls.run_bwa_idex(reference_genome, outdir)
    f = time.time()
    if verbose is True:
        print(" > total execution time : ", f - s)


def run_step_1(
    fastq1,
    fastq2,
    adapters,
    prefixout,
    threads,
    min_len,
    trim,
    outdir,
    verbose=True,
):
    """
    Run step 1 of ViralFlow, check github README and/or paper for further
    details
    """
    if verbose is True:
        print(">--- STEP 1 ---<")
    # --- sanity check of arguments ----------------------------------------
    try:
        assert Path(outdir).is_dir()
    except (AssertionError):
        print(
            "WARNING: "
            + outdir
            + " does not exist. Don't worry it will be created then. "
        )
        os.system("mkdir " + outdir)
    # ----------------------------------------------------------------------

    if verbose is True:
        print("@ running fastp...")
    s = time.time()
    viralflow.calls.run_fastp(
        fastq1, fastq2, prefixout, threads, adapters, min_len, trim, outdir
    )
    f = time.time()
    if verbose is True:
        print(" > total execution time : ", f - s)
    # TODO : check output for possible errors ################################


def run_step_2(ref_gnm, prefixout, threads, outdir, depth, verbose=True):
    """
    Run step 2 of ViralFlow, check github README and/or paper for further
    details
    """
    if verbose is True:
        print(">--- STEP 2 ---<")

        print("@ aligning reads to reference genome...")
    fq1 = outdir + prefixout + ".R1.fq.gz"
    fq2 = outdir + prefixout + ".R2.fq.gz"

    s = time.time()
    viralflow.calls.align_reads2ref(ref_gnm, fq1, fq2, prefixout, threads, outdir)
    f = time.time()
    if verbose is True:
        print(" > total execution time : ", f - s)
        print("@  Get consensus and minor variant analysis...")

    s = time.time()
    viralflow.calls.run_ivar(ref_gnm, prefixout, outdir, depth)
    f = time.time()

    if verbose is True:
        print(" > total execution time : ", f - s)
        print("@ get minor variants data...")

    viralflow.calls.get_minor_variants(ref_gnm, prefixout, depth, threads, outdir)
    # TODO check output for possible errors ##################################


def run_step_3(
    prefixout,
    depth,
    intrahost_depth,
    out_dir,
    depth_number=100,
    per_limit=0.05,
    verbose=True,
):
    """
    Run step 3 of ViralFlow, check github README and/or paper for further
    details
    """
    # --- sanity check of arguments ----------------------------------------
    try:
        assert Path(out_dir).is_dir()
    except (AssertionError):
        print(
            "WARNING: "
            + out_dir
            + " does not exist. Don't worry it will be created then. "
        )
        os.system("mkdir " + out_dir)
    # ----------------------------------------------------------------------
    if verbose is True:
        print(">--- STEP 3 ---<")

        print("@ minor variants identification...")
        print("   > minor variants only if depth >= ", depth_number)
        print("   > putative minor if Depth Base / Total Base >=", per_limit)
        print("   > INPUT FILES:")
    prfx_wdir = out_dir + prefixout
    bam_rc_file = f"{prfx_wdir}.depth{depth}.fa.bc"
    alignment_file = f"{prfx_wdir}.depth{depth}.fa.algn"
    # sanity check ----------------------------------------------------------
    assert Path(bam_rc_file).is_file(), bam_rc_file + " does not exist."
    assert Path(alignment_file).is_file(), alignment_file + " does not exist."
    # -----------------------------------------------------------------------
    if verbose is True:
        print("     - ", bam_rc_file)
        print("     - ", alignment_file)

    s = time.time()
    # WARNING -----------------------------------------------------------------
    # For now we gonna keep using the intrahost_script, but this should be
    # rewritten as a callable function
    os.system(
        f"python {VIRALFLOW_PATH}/viralflow/intrahost_script.py \
    -in {prfx_wdir}.depth{depth}.fa.bc -al {prfx_wdir}.depth{depth}.fa.algn \
    -dp {intrahost_depth}"
    )
    # TODO : Filipe, this one is yours. #######################################
    # viralflow.intrahost.get_instrahost_tsv(bam_rc_file, alignment_file,
    #                                       per_limit=per_limit,
    #                                       depth_value=depth_number)
    f = time.time()
    if verbose is True:
        print(" > total execution time : ", f - s)
    # TODO: check output for possible errors ##################################


def run_step_4(
    prefixout,
    depth,
    threads,
    out_dir,
    ref_gnm,
    nxtdst_str,
    nxt_bin="nextclade",
    verbose=True,
):
    """
    Run step 4 of ViralFlow, check github README and/or paper for further
    details

    """
    if verbose is True:
        print(">--- STEP 4 ---<")

        print("@ Naming variants identified...")
    viralflow.calls.get_variant_naming(
        prefixout, depth, threads, out_dir, ref_gnm, nxtdst_str, nxt_bin=nxt_bin
    )
    # TODO: check output for possible errors ##################################


def run_step_5(prefixout, outdir, verbose=True):
    """ """
    if verbose is True:
        print(">--- STEP 5 ---<")
        print("@ computings assembly metrics...")
    viralflow.calls.do_assembly_metrics(prefixout, outdir)
    # TODO: check output for possible errors #################################


def run_pipeline(
    outdir,
    ref_gnm,
    fastq1,
    fastq2,
    adapters,
    prefixout,
    threads,
    depth,
    min_len,
    trim,
    intrahost_depth,
    nxt_dataset,
    nxt_bin="nextclade",
    verbose=True,
    skip_0=False,
):
    """
    Run the full pipeline on a pair of fastq
    """
    try:
        assert Path(outdir).is_dir()
    except (AssertionError):
        print(
            "WARNING: "
            + outdir
            + " does not exist. Don't worry it will be created then. "
        )
        os.system("mkdir " + outdir)
    if skip_0 is False:
        run_step_0(ref_gnm, outdir, prefixout, verbose=verbose)

    run_step_1(
        fastq1,
        fastq2,
        adapters,
        prefixout,
        threads,
        min_len,
        trim,
        outdir,
        verbose=verbose,
    )

    run_step_2(ref_gnm, prefixout, threads, outdir, depth, verbose=verbose)

    run_step_3(prefixout, depth, intrahost_depth, outdir, verbose=verbose)

    run_step_4(
        prefixout,
        depth,
        threads,
        outdir,
        ref_gnm,
        nxt_dataset,
        nxt_bin=nxt_bin,
        verbose=verbose,
    )

    run_step_5(prefixout, outdir, verbose=verbose)