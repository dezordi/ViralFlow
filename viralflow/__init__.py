# path handling stuff
import os
from pathlib import Path

# multiprocessing stuff
import multiprocessing as mp
import tqdm
import viralflow.containers
import viralflow.calls
import viralflow.pipeline
import viralflow.intrahost
import viralflow.output

# --- FUNCTIONS ---------------------------------------------------------------


def get_fastq_pairs(input_dir, format=tuple(["fastq.gz", "fq.gz"])):
    """
    Get fastq pairs files (R1 and R2) on a given dir.
    This function assumes pairs have the same name differing only that instead
    by the 'R1' and 'R2'.

    Parameters
    ----------
    input_dir:<str>
        path to directory containing input files

    format:<str>
        fastq format (default = 'fastq.gz')

    Return
    ------
    <lst>
        list of fastq pairs detected
    """
    # --- LOCAL FUNCTION ------------------------------------------------------
    def get_R2_file(fastq_R1):
        """
        get R2 pair name and check if file actually exist

        Parameters
        ----------
        fastq_R1:<str>
            R1 fastq filename

        Returns
        -------
        <str>
            input filename, but 'R1' is changed to 'R2'
        """
        # get index of 'R1' at fastq file name
        R_idx = fastq_R1.find("R1")

        # get 'R2' file
        new_name = ""
        for i, char in enumerate(fastq_R1):
            if i is R_idx:
                new_name += "R"
                continue
            if i is R_idx + 1:
                new_name += "2"
                continue
            new_name += char
        return new_name

    # -------------------------------------------------------------------------

    # get fastq.gz files
    files_at_dir = os.listdir(input_dir)
    fastq_lst = [x for x in files_at_dir if x.endswith(format)]

    # get R1 and R2 pairs
    pairs = []
    for fastq in fastq_lst:
        # check if R1, if true get R2 filename
        if "R1" in fastq:
            fastq_R2_fl = get_R2_file(fastq)
            try:
                assert Path(input_dir + fastq_R2_fl).is_file()
            except (AssertionError):
                msg = (
                    fastq_R2_fl + " does not exist, but is the expected R2 for " + fastq
                )
                raise Exception(msg)
            # store pairs
            pairs.append([fastq, fastq_R2_fl])

    return pairs


def run_sing_cont_pp(kwargs):
    """
    wrapper of containers run sing container function
    """
    container_img = kwargs["container_img"]
    input_dir = kwargs["input_dir"]
    ref_gnm = kwargs["ref_gnm"]
    prefix_out = kwargs["prefix_out"]
    adapters_file = kwargs["adapters_file"]
    threads = kwargs["threads"]
    depth = kwargs["depth"]
    min_len = kwargs["min_len"]
    min_dp_intrahost = kwargs["min_dp_intrahost"]
    trim_len = kwargs["trim_len"]
    sing_call = kwargs["sing_call"]
    fastq_R1 = kwargs["fastq_R1"]
    fastq_R2 = kwargs["fastq_R2"]

    viralflow.containers.run_sing_container(
        container_img,
        input_dir,
        ref_gnm,
        fastq_R1,
        fastq_R2,
        prefix_out,
        adapters_file,
        threads=threads,
        depth=depth,
        min_len=min_len,
        min_dp_intrahost=min_dp_intrahost,
        trim_len=trim_len,
        sing_call=sing_call,
        dry=False,
    )


def run_pipeline_pp(kwargs):
    """
    wrapper for run pipeline function
    """
    outdir = kwargs["outdir"]
    ref_gnm = kwargs["ref_gnm"]
    fastq1 = kwargs["fastq_R1"]
    fastq2 = kwargs["fastq_R2"]
    adapters = kwargs["adapters_file"]
    prefixout = kwargs["prefixout"]
    threads = kwargs["threads"]
    depth = kwargs["depth"]
    min_len = kwargs["min_len"]
    trim = kwargs["trim_len"]
    intrahost_depth = kwargs["min_dp_intrahost"]
    nxt_dataset = kwargs["nxt_dataset"]
    verbose = kwargs["verbose"]
    nxt_bin = kwargs["nxt_bin"]

    viralflow.pipeline.run_pipeline(
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
        nxt_bin=nxt_bin,
        verbose=verbose,
        skip_0=True,
    )


def run_viralflow_pp(
    input_dir,
    ref_gnm,
    adapters_file,
    depth,
    min_len,
    min_dp_intrahost,
    trim_len,
    nxt_dataset,
    nxt_bin="nextclade",
    cpus_total=mp.cpu_count(),
    cpus_pprc=None,
    verbose=False,
):
    """
    Run multiple viralflow call to all fastq pais (R1 and R2) on a given
    directory.
    """
    # santy check --------------------------------------------------------------
    # be sure the input is a valid directory
    if input_dir.endswith("/") is False:
        input_dir += "/"
    assert Path(input_dir).is_dir(), input_dir + " is not a dir"

    # check if files exists
    files = [adapters_file, ref_gnm]
    for f in files:
        assert Path(input_dir + f).exists(), f + " does not exist at " + input_dir

    # --------------------------------------------------------------------------
    # get pairs of fastq available
    print("@ mapping fastq pairs to process...")
    pairs_lst = viralflow.get_fastq_pairs(input_dir)
    assert len(pairs_lst) > 0, "no fastq pairs detected "
    n_runs = len(pairs_lst)
    print("  >> ", n_runs, " total samples detected")

    # get number of parallel run calls by considering the total number of cpus
    # available and number of cpus per run
    print("@ multiprocessing prep")
    # total number of cpus to use
    cpus_avlbl = mp.cpu_count()
    if cpus_total > cpus_avlbl:
        print(
            "WARNING: ",
            cpus_total,
            " requested but only ",
            cpus_avlbl,
            "are available. Setting total cpus as ",
            cpus_avlbl,
        )
        cpus_total = cpus_avlbl

    print("  >> total cpus = ", cpus_total)
    # total number of cpus per simultaneos process

    if cpus_pprc not in [None, 0]:
        # check if multiple call number is higher then number of samples
        assert type(cpus_pprc) == int, "cpus_pprc must be a integer"
        n_parallel_calls = round(cpus_total / cpus_pprc)
        if n_parallel_calls > n_runs:
            print(
                "WARNING: the ratio of total cpus and cpus per sample",
                f"requested({cpus_total}/{cpus_pprc}={n_parallel_calls})",
                f" is bigger then the total number of samples ({n_runs}).",
                "Setting cpus per samples to AUTO",
            )
            cpus_pprc = 0

    if cpus_pprc in [None, 0]:
        cpus_pprc = round(cpus_total / n_runs)

    if cpus_pprc > cpus_total:
        print(
            "WARNING: ",
            cpus_pprc,
            " requested but only ",
            cpus_total,
            " will be used. Setting total cpus per process automatically",
        )
        cpus_pprc = round(cpus_total / n_runs)

    if cpus_pprc < 1:
        cpus_pprc = 1
    print("  >> cpus per call = ", cpus_pprc)

    n_parallel_calls = round(cpus_total / cpus_pprc)
    print("  > number of multiple call = ", n_parallel_calls)

    # get keyword arguments per pairs of fastq
    dct_fix = {
        "input_dir": input_dir,
        "ref_gnm": input_dir + ref_gnm,
        "adapters_file": input_dir + adapters_file,
        "threads": cpus_pprc,
        "depth": depth,
        "min_len": min_len,
        "min_dp_intrahost": min_dp_intrahost,
        "trim_len": trim_len,
        "nxt_dataset": nxt_dataset,
        "nxt_bin": nxt_bin,
        "verbose": verbose,
    }

    # do the ref genome indexing
    rgm = input_dir + ref_gnm
    viralflow.pipeline.run_step_0(rgm, input_dir, verbose=verbose)
    # get kwargs list
    kwargs_lst = []

    for pair in pairs_lst:
        R1 = pair[0]
        R2 = pair[1]
        prefix_out = R1.split(".")[0]
        kwarg = {
            "fastq_R1": input_dir + R1,
            "fastq_R2": input_dir + R2,
            "prefixout": prefix_out,
            "outdir": input_dir + prefix_out + ".results/",
        }
        kwarg_f = {**kwarg, **dct_fix}
        kwargs_lst.append(kwarg_f)

    # run multiple calls in parallels
    workers = mp.Pool(n_parallel_calls)
    max_ = len(kwargs_lst)

    with tqdm.tqdm(total=max_) as pbar:
        for i, chunk_results in enumerate(
            workers.imap(run_pipeline_pp, kwargs_lst)  # run_sing_cont_pp,
        ):
            pbar.update()
        workers.close()
