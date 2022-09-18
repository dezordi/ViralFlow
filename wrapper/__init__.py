from logging import root
import os

def build_containers(root_path):
    """
    run script to build container for vfnext
    """
    cd_to_dir= f"cd {root_path}/vfnext/containers/" 
    rm_old_cntnrs = "rm *.sif"
    run_bash = f"bash ./setupContainers.sh"
    print(cd_to_dir+';'+rm_old_cntnrs+';'+run_bash)
    os.system(cd_to_dir+';'+rm_old_cntnrs+';'+run_bash)

def install_dependencies(root_path):
    """
    install nextflow and singularity via conda
    """
    run_bash = f"bash {root_path}/wrapper/install_dep.sh"
    print(run_bash)
    os.system(run_bash)

# input args file load
def parse_params(in_flpath):
    """
    load text file containing viralflow arguments
    """
    valid_args = [
        "virus",
        "adaptersFile",
        "outDir",
        "inDir",
        "runSnpEff",
        "writeMappedReads",
        "minLen",
        "depth",
        "minDpIntrahost",
        "trimLen",
        "runSnpEff",
        "refGenomeCode",
        "referenceGFF",
        "referenceGenome",
        "nextflowSimCalls",
        "fastp_threads",
        "bwa_threads",
        "mafft_threads",
        "nxtclade_jobs"
    ]
    path_params = ["inDir", "outDir", "referenceGFF", "referenceGenome", "adaptersFile"]
    in_file = open(in_flpath, "r")
    dct = {}
    for l in in_file:
        # skip lines
        if (l in ["", " ", "\n"]) or l.startswith("#"):
            continue

        # get line data
        l_dt = l.replace("\n", "").split(" ")
        
        # get content
        key = l_dt[0]
        if (key not in valid_args):
            raise Exception(f"ERROR: {key} not a valid argument")
        # fill dict
        if key in valid_args:
        
            vls_1 = l_dt[1 : len(l_dt)]
            vls = []
        
            for v in vls_1:
                if v in [""]:
                    continue
                vls.append(v)
            # if single value
            if len(vls) == 1:
                # skip null values
                if vls[0] == "null":
                    continue
                # be sure paths are absolute
                if key in path_params:
                    dct[key] = os.path.abspath(vls[0])
                    continue
                dct[key] = vls[0]
            # if a list of values
            if len(vls) > 1:
                dct[key] = vls
            continue
    # get arguments for nextflow
    
    args_str = ""
    for key in dct:
        args_str += f"--{key} {dct[key]} "
    args_str += "-resume"
    return args_str


def run_vfnext(root_path, params_fl):
    # get nextflow arguments
    args_str = parse_params(params_fl)
    nxtflw_ver="22.04.0"
    run_nxtfl_cmd = f"NXF_VER={nxtflw_ver} nextflow run {root_path}/vfnext/main.nf {args_str}"
    print(run_nxtfl_cmd)
    os.system(run_nxtfl_cmd)
