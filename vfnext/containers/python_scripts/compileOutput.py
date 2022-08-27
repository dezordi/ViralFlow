from ast import Assert
import os
import pandas as pd
import sys
from Bio import SeqIO
import argparse

__author__ = "Antonio Marinho"
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Antonio Marinho"
__email__ = "amarinhosn@pm.me"
__date__ = "2022/06/07"
__username__ = "amarinhosn"


# --- FUNCTIONS ---------------------------------------------------------------
def get_consensus_seq(fasta_path, cod):
    """
    get sequence line from a single sequence fasta file
    """
    with open(fasta_path, "r") as fasta_fl:
        for line in fasta_fl:
            if line.startswith(">"):
                # sanity check
                assert cod in line
            else:
                return line

def parse_fabc(fbc_path):
    """
    parse depth data on 'fa.bc' file
    """
    # load file
    fl = open(fbc_path, "r")
    # get only pos and depth per line and return data as dataframe
    pos_lst = []
    depth_lst = []
    for l in fl:
        l_dt = l.split("\t")
        pos_lst.append(int(l_dt[1]))
        depth_lst.append(int(l_dt[3]))

    df = pd.DataFrame.from_dict({"pos": [pos_lst], "depth": [depth_lst]})
    return df

def __check_if_outfls_exist(cod, out_fls_lst, mut_fl):
    """
    Check if expected files exist
    """
    # for a list of expected files, check if they exist at specified path
    # if not, then print a warning, store missing error data on a dictionary
    # and return it
    err_dct_lst = []
    for fl in out_fls_lst:
        try:
            assert os.path.isfile(fl)
        except (AssertionError):
            srvc = fl.split(".")[-2]
            if fl == mut_fl:
                srvc = "Mutations"
            print(f"  :: WARNING: missing {srvc} output for {cod} [expects {fl}]")
            type = "ERROR"
            kind = f"No {srvc} output [expected={fl}]"
            err_dct = {"cod": cod, "type": type, "kind": kind}
            err_dct_lst.append(err_dct)
            continue
    return err_dct_lst

def __parse_xsv(cod, xsv_path, sep=","):
    """
    get dataframe from tabular data (ex: csv and tsv)
    """
    df = pd.read_csv(xsv_path, sep=sep)
    df["cod"] = [cod for i in range(0, len(df))]
    return df

def get_mut(row):
    """
    get mutation data on a single string
    """
    return row["REF"] + str(row["POS"]) + row["ALT"]


def compile_output_fls(data_dir, out_dir, depth):
    """
    Compile Viralflow output data files on subdirs
    Parameters
    ==========
    data_dir:<path>
        path for dir containing results output subdirs
    out_dir:<path>
        output dir to write compiled data
    depth:<depth>
        depth prefix used
    """
    # create outdir if does not exist
    try:
        assert os.path.exists(out_dir)
    except (AssertionError):
        print("@ creating output dir")
        os.system("mkdir " + out_dir)
    if out_dir.endswith("/") is False:
        out_dir += "/"
    print("@ compiling output files")
    # create multifasta file
    multifas_fl = open(out_dir + "seqbatch.fa", "w")

    # create list to store data
    pango_df_lst = []
    chrms_df_lst = []
    nxtcd_df_lst = []
    ithst_df_lst = []
    mut_df_lst = []
    err_dct_lst = []
    skip_lst = []
    dpth_lst = []
    # iterate over subdirs on data dir
    c = 0
    for path, subdirs, files in os.walk(data_dir):
        for name in files:
            # check only files at *.results subdirs
            sub_dir = path.split("/")[-1]

            if sub_dir.endswith("_results") is False:
                continue

            # get consensus
            cod = sub_dir.split("_")[0]
            prfx = path + f"/{cod}."
            if cod in skip_lst:
                continue

            # look for consensus sequence fasta, if yes, then check for
            # 1) expected files
            # 2) parse and compiled results in single files
            if name.endswith(f"depth{depth}.fa"):
                c += 1
                print(f"  > {c} samples processed", end="\r")
                fasta_path = os.path.join(path, name)
                # sanity
                assert(cod in name)
                # write sequence to new
                seq = get_consensus_seq(fasta_path, cod)
                if len(seq) == 1:
                    type = "ERROR"
                    kind = "No consensus sequence obtained"
                    print(f"  :: WARNING: No consensus sequence obtained for {cod}")
                    err_dct_lst.append({"cod": cod, "type": type, "kind": kind})
                    skip_lst.append(cod)
                    continue
                # write to multifasta
                multifas_fl.write(">" + cod + "\n")
                multifas_fl.write(seq + "\n")
                # continue
                # -- SET EXPECTED OUTPUT FILES --------------------------------
                # depth data
                depth_fl = f"{prfx}depth{depth}.fa.bc"
                # pangolin
                try:
                    pango_fl = f"{prfx}fa.pango.out.csv"
                    assert os.path.isfile(pango_fl)
                except (AssertionError):
                    pango_fl = f"{prfx}all.fa.pango.out.csv"
                # chromosomes
                #chrms_fl = f"{path}/chromosomes.report"
                # intrahost
                ithst_fl = f"{prfx}depth{depth}.fa.bc.intrahost.tsv"
                # mutation
                mut_fl = f"{prfx}tsv"
                # nextclade
                try:
                    nxtcl_fl = f"{prfx}depth{depth}.fa.nextclade.csv"
                    assert os.path.isfile(nxtcl_fl)
                except (AssertionError):
                    nxtcl_fl = f"{prfx}depth{depth}.all.fa.nextclade.csv"

                # --- CHECK IF FILES EXIST ------------------------------------
                out_fls_lst = [pango_fl, ithst_fl, mut_fl]#,depth_fl, chrms_fl]
                err_dct_ = __check_if_outfls_exist(cod, out_fls_lst, mut_fl)
                # if any file is missing, skip current sample
                if len(err_dct_) > 0:
                    for err in err_dct_:
                        err_dct_lst.append(err)
                    break

                # --- PARSE DATA ----------------------------------------------
                # pangolin
                pango_df_lst.append(__parse_xsv(cod, pango_fl, sep=","))
                # chromossomes
                #chrms_df_lst.append(__parse_xsv(cod, chrms_fl, sep="\t"))
                # nextclade
                nxtcd_df_lst.append(__parse_xsv(cod, nxtcl_fl, sep=";"))
                # intrahost
                ithst_df_lst.append(__parse_xsv(cod, ithst_fl, sep="\t"))

                # Mutations
                mut_df = pd.read_csv(mut_fl, sep="\t")
                val_mut_df = mut_df.loc[mut_df["PASS"] == True]
                # check if no mutation data
                if len(val_mut_df) == 0:
                    print(f"WARNING: no mutation data for {cod}")
                    mut_lst = []
                if len(val_mut_df) > 0:
                    mut_lst = val_mut_df.apply(get_mut, axis=1).values

                dct = {"cod": cod, "mut": mut_lst}
                mut_df_lst.append(pd.DataFrame(dct))
                # depth data
                dpth_df = parse_fabc(depth_fl)
                dpth_df["cod"] = [cod]
                dpth_lst.append(dpth_df)

    print(f"  > Total {c} samples processed")

    # mount compiled dataframes
    if c == 0:
        print("ERROR: No ViralFlow output files found")
        sys.exit(1)
    if c > 0:
        print("@ writing compiled data")

        try:
            assert(len(pango_df_lst)>0)
            all_pango_df = pd.concat(pango_df_lst, ignore_index=True)
            all_pango_df.to_csv(out_dir + "/pango.csv", index=False)
            print("  > pango.csv")
        except(AssertionError):
            print("WARN: No data from pangolin")
        
        try:
            assert(len(nxtcd_df_lst)>0)
            all_nxtcd_df = pd.concat(nxtcd_df_lst, ignore_index=True)
            all_nxtcd_df.to_csv(out_dir + "/nextclade.csv", index=False)
            print("  > nextclade.csv")
        except(AssertionError):
            print("WARN: No data from Nextclade")
        
        try:
            assert(len(mut_df_lst)>0)
            all_mut_df = pd.concat(mut_df_lst, ignore_index=True)
            all_mut_df.to_csv(out_dir + "/mutations.csv", index=False)
            print("  > mutations.csv")
        except(AssertionError):
            print("WARN: No mutation data")

        errors_df = pd.DataFrame(err_dct_lst)
        errors_df.to_csv(out_dir + "/errors_detected.csv", index=False)
        print("  > errors_detected.csv")

        # write csvs
        #all_chrms_df = pd.concat(chrms_df_lst, ignore_index=True)
        #all_dpth_df = pd.concat(dpth_lst, ignore_index=True)
        #all_chrms_df.to_csv(out_dir + "/chromossomes.csv", index=False)
        #print("  > chromossomes.csv")
        #all_dpth_df.to_csv(out_dir + "/depth.csv", index=False)
        #print("  > depth.csv")
        
    print(" :: DONE ::")
# -----------------------------------------------------------------------------

if __name__=="__main__":
  parser = argparse.ArgumentParser(
   description = 'This scripts compile outputs of samples processed on viralow on single files')
  parser.add_argument("-dD","--dataDir", required=True,
    help="directory containing viralflow output files")
  parser.add_argument("-oD","--outputDir", required=True,
    help="directory to output compiled files")
  parser.add_argument("-dp","--depth", type=int, default= 5,
    help="Minimum depth value to consider set as input to ViralFlow as minor variant, default = 5")
  args = parser.parse_args()

  compile_output_fls(args.dataDir, args.outputDir, args.depth)
