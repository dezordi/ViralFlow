from distutils.command.build_scripts import first_line_re
from logging import root
import os


def add_entries_to_DB(root_path, org_name, refseq_code):
    """
    add entries provided to snpeff database
    """
    run_bash = f"bash {root_path}/vfnext/containers/add_entries_SnpeffDB.sh"
    print(f"{run_bash} {org_name} {refseq_code}")
    os.system(f"{run_bash} {org_name} {refseq_code}")

def parse_csv(csv_flpath):
    with open(csv_flpath, "r") as csv_fl:
        first_line = True
        entries_lst = []
        for line in csv_fl:
            # skip header
            if first_line == True:
                first_line = False
                continue
            ln_data = line.split(",")
            entry = [ln_data[0], ln_data[1].replace("\n","")]
            entries_lst.append(entry)
    return entries_lst

def build_containers(root_path):
    """
    run script to build container for vfnext
    """
    # build containers
    cd_to_dir= f"cd {root_path}/vfnext/containers/" 
    rm_old_cntnrs = "rm -rf *.sif"
    run_bash = f"bash ./setupContainers.sh"
    print(cd_to_dir+';'+rm_old_cntnrs+';'+run_bash)
    os.system(cd_to_dir+';'+rm_old_cntnrs+';'+run_bash)
    # add new entries to snpeff
    entries = parse_csv(f"{root_path}/vfnext/databases/add_to_snpeff_db.csv")
    print("@ add entries to SnpEff DB")
    for entry in entries:
        print(f"..> {entry[0]}:{entry[1]}")
        add_entries_to_DB(root_path, entry[0], entry[1])
    
def install_dependencies(root_path):
    """
    install nextflow and singularity via conda
    """
    # install software dependencies
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

def update_pangolin(root_path):
    cd_to_dir= f"cd {root_path}/vfnext/containers/" 
    run_update = "singularity exec ./pangolin_latest.sif pangolin --update"
    os.system(cd_to_dir+';'+run_update)

def run_vfnext(root_path, params_fl):
    # get nextflow arguments
    args_str = parse_params(params_fl)
    nxtflw_ver="22.04.0"
    run_nxtfl_cmd = f"NXF_VER={nxtflw_ver} nextflow run {root_path}/vfnext/main.nf {args_str}"
    print(run_nxtfl_cmd)
    os.system(run_nxtfl_cmd)
