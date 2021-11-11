import os
import pandas as pd
import sys
# get list of results dir
# === FUNCTIONS ===============================================================
# --- Compile results ---------------------------------------------------------


def load_short_summary_df(chrm_df, pang_df):
    '''
    '''
    # get slices
    chrms_slice = chrm_df[['cod', ' Avg depth', ' Coverage%']]
    pang_slice = pang_df[['cod', 'taxon', 'lineage']]
    # merge dataframes
    short_summary = pd.merge(chrms_slice, pang_slice,
                             right_on='cod', left_on='cod')
    # return df
    return short_summary


def get_consensus_seq(fasta_path, cod):
    '''
    get sequence line from a single sequence fasta file
    '''
    with open(fasta_path, 'r') as fasta_fl:
        for line in fasta_fl:
            if line.startswith('>'):
                # sanity check
                assert(cod in line)
            else:
                return line


def get_mut(row):
    return row['REF']+str(row['POS'])+row['ALT']


def parse_fabc(fbc_path):
    '''
    parse depth data on 'fa.bc' file
    '''
    fl = open(fbc_path, 'r')

    pos_lst = []
    depth_lst = []
    for l in fl:
        l_dt = l.split('\t')
        pos_lst.append(int(l_dt[1]))
        depth_lst.append(int(l_dt[3]))

    df = pd.DataFrame.from_dict({'pos': [pos_lst],
                                 'depth': [depth_lst]})
    return df


def __parse_xsv(cod, xsv_path, sep=','):
    '''
    get dataframe from tabular data (ex: csv and tsv)
    '''
    df = pd.read_csv(xsv_path, sep=sep)
    df['cod'] = [cod for i in range(0, len(df))]
    return df


def __check_if_outfls_exist(cod, out_fls_lst, mut_fl):
    '''
    Check if expected files exist
    '''
    err_dct_lst = []
    for fl in out_fls_lst:
        try:
            assert(os.path.isfile(fl))
        except(AssertionError):
            srvc = fl.split('.')[-2]
            if fl == mut_fl:
                srvc = 'Mutations'
            print(
                f'  :: WARNING: missing {srvc} output for {cod} [expects {fl}]'
            )
            type = 'ERROR'
            kind = f'No {srvc} output [expected={fl}]'
            err_dct = {'cod': cod, 'type': type, 'kind': kind}
            err_dct_lst.append(err_dct)
            continue
    return err_dct_lst


def compile_output_fls(data_dir, out_dir, depth):
    '''
    Compile Viralflow output data files on subdirs
    Parameters
    ==========
    data_dir:<path>
        path for dir containing results output subdirs
    out_dir:<path>
        output dir to write compiled data
    depth:<depth>
        depth prefix used
    '''
    # create outdir if does not exist
    try:
        assert(os.path.exists(out_dir))
    except(AssertionError):
        print('@ creating output dir')
        os.system('mkdir '+out_dir)
    if out_dir.endswith('/'):
        out_dir += '/'
    print('@ compiling output files')
    # create multifasta file
    multifas_fl = open(data_dir+'seqbatch.fa', 'w')

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
            sub_dir = path.split('/')[-1]

            if sub_dir.endswith('.results') is False:
                continue
            # get consensus
            cod = sub_dir.split('.')[0]
            prfx = path+f'/{cod}.'
            if cod in skip_lst:
                continue

            # look for consensus sequence fasta, if yes, then check for
            # 1) expected files
            # 2) parse and compiled results in single files
            if name.endswith(f'depth{depth}.fa'):
                c += 1
                print(f'  > {c} samples processed', end='\r')
                fasta_path = os.path.join(path, name)
                # sanity
                assert(cod in name)
                # write sequence to new
                seq = get_consensus_seq(fasta_path, cod)
                if len(seq) == 1:
                    type = 'ERROR'
                    kind = 'No consensus sequence obtained'
                    print(
                      f'  :: WARNING: No consensus sequence obtained for {cod}'
                    )
                    err_dct_lst.append(
                        {'cod': cod, 'type': type, 'kind': kind})
                    skip_lst.append(cod)
                    continue
                # write to multifasta
                multifas_fl.write('> '+cod+'\n')
                multifas_fl.write(seq+'\n')
                # continue
                # -- SET EXPECTED OUTPUT FILES --------------------------------
                # depth data
                depth_fl = f'{prfx}depth{depth}.fa.bc'
                # pangolin
                try:
                    pango_fl = f'{prfx}depth{depth}.fa.pango.csv'
                    assert(os.path.isfile(pango_fl))
                except(AssertionError):
                    pango_fl = f'{prfx}depth{depth}.all.fa.pango.csv'
                # chromosomes
                chrms_fl = f'{path}/chromosomes.report'
                # intrahost
                ithst_fl = f'{prfx}depth{depth}.fa.bc.intrahost.tsv'
                # mutation
                mut_fl = f'{prfx}tsv'
                # nextclade
                try:
                    nxtcl_fl = f'{prfx}depth{depth}.fa.nextclade.csv'
                    assert(os.path.isfile(nxtcl_fl))
                except(AssertionError):
                    nxtcl_fl = f'{prfx}depth{depth}.all.fa.nextclade.csv'

                # --- CHECK IF FILES EXIST ------------------------------------
                out_fls_lst = [pango_fl, chrms_fl, ithst_fl, mut_fl, depth_fl]
                err_dct_ = __check_if_outfls_exist(cod, out_fls_lst, mut_fl)
                # if any file is missing, skip current sample
                if len(err_dct_) > 0:
                    for err in err_dct_:
                        err_dct_lst.append(err)
                    break

                # --- PARSE DATA ----------------------------------------------
                # pangolin
                pango_df_lst.append(__parse_xsv(cod, pango_fl, sep=','))
                # chromossomes
                chrms_df_lst.append(__parse_xsv(cod, chrms_fl, sep='\t'))
                # nextclade
                nxtcd_df_lst.append(__parse_xsv(cod, nxtcl_fl, sep=';'))
                # intrahost
                ithst_df_lst.append(__parse_xsv(cod, ithst_fl, sep='\t'))

                # Mutations
                mut_df = pd.read_csv(mut_fl, sep='\t')
                val_mut_df = mut_df.loc[mut_df['PASS'] == True]
                # check if no mutation data
                if len(val_mut_df) == 0:
                    print(f"WARNING: no mutation data for {cod}")
                    mut_lst = []
                if len(val_mut_df) >= 0:
                    mut_lst = val_mut_df.apply(get_mut, axis=1).values

                dct = {'cod': cod, 'mut': mut_lst}
                mut_df_lst.append(pd.DataFrame(dct))
                # depth data
                dpth_df = parse_fabc(depth_fl)
                dpth_df['cod'] = [cod]
                dpth_lst.append(dpth_df)
        # if name.endswith('intrahost.short.tsv')
    print(f'  > Total {c} samples processed')

    # mount compileds dataframes
    if c == 0:
        print("ERROR: No ViralFlow output files found")
        sys.exit(1)
    if c > 0:
        print('@ writing compiled data')
        all_pango_df = pd.concat(pango_df_lst, ignore_index=True)
        all_nxtcd_df = pd.concat(nxtcd_df_lst, ignore_index=True)
        all_chrms_df = pd.concat(chrms_df_lst, ignore_index=True)
        all_mut_df = pd.concat(mut_df_lst, ignore_index=True)
        all_dpth_df = pd.concat(dpth_lst, ignore_index=True)
        errors_df = pd.DataFrame(err_dct_lst)
        # write csvs
        all_pango_df.to_csv(out_dir+'/pango.csv')
        print('  > pango.csv')
        all_nxtcd_df.to_csv(out_dir+'/nextclade.csv')
        print('  > nextclade.csv')
        all_chrms_df.to_csv(out_dir+'/chromossomes.csv')
        print('  > chromossomes.csv')
        all_mut_df.to_csv(out_dir+'/mutations.csv')
        print('  > mutations.csv')
        all_dpth_df.to_csv(out_dir+'/depth.csv')
        print('  > depth.csv')
        errors_df.to_csv(out_dir+'/errors_detected.csv')
        print('  > errors_detected.csv')

    print(' :: DONE ::')

# --- Check controls ----------------------------------------------------------


def get_controls_lineages(control_lbl_lst, pango_df):
    '''
    get control unique lineages set
    '''
    def __get_ctrls(row):
        if row['cod'] in control_lbl_lst:
            return True
        else:
            return False
    ctrl_bool_tbl = pango_df.apply(__get_ctrls, axis=1)
    control_slice_df = pango_df.loc[ctrl_bool_tbl]
    control_lngs = control_slice_df.lineage.unique()
    return control_lngs


def sel_samples_ofLngs(pango_df, lng_lst):
    '''
    Get samples with a given set of lineages
    '''
    def get_samples_wCtrls_lng(row):
        if row['lineage'] in lng_lst:
            return True
        else:
            return False
    sample_wCtrl_lng_df = pango_df.loc[pango_df.apply(
        get_samples_wCtrls_lng, axis=1)]
    return sample_wCtrl_lng_df


def check_controls_lineages(control_lbl_lst, pango_csv):
    '''
    check if controls have lineages assigned.
    '''
    print('@ checking controls lineages...')
    # get controls lineages
    pango_df = pd.read_csv(pango_csv)
    control_lngs = get_controls_lineages(control_lbl_lst, pango_df)
    n_ctrl_lngs = len(control_lngs)
    if n_ctrl_lngs == 0:
        print('No lineages on negative controls')
        return None
    if n_ctrl_lngs > 0:
        # check if control lineages are different from samples (is it unique?)
        print('Lineages found on negative controls')
        print(f'  > {n_ctrl_lngs} found : {control_lngs}')
        print('@ looking for contamination evidence...')
        # check if lineages are present on samples
        suspected_df = sel_samples_ofLngs(pango_df, control_lngs)
        n_suspected = len(suspected_df)
        if n_suspected == 0:
            print(
                '  > No sample presenting same lineages of control. Everything is fine.')
            return None
        if n_suspected > 0:
            print(f'  > {n_suspected} suspected samples')
            print('    | cod | taxon | lineage ')
            for i in suspected_df[['cod', 'taxon', 'lineage']].values:
                print(f'    |{i[0]} | {i[1]} | {i[2]} ')
        return suspected_df

# --- Data summary ------------------------------------------------------------


def get_lineages_summary(pango_csv, chromosomes_csv, outdir):
    '''
    count lineages on a given pangolin dataframe
    '''
    print('@ compute lineage summary ')
    pango_df = pd.read_csv(pango_csv)
    print(f'  > {len(pango_df)} total samples')
    lineage_df = pango_df['lineage'].value_counts()
    lineage_df.to_csv(outdir+'/lineage_summary.csv')
    print(f'  > {outdir}lineage_summary.csv')
    print(lineage_df)

    print("@ generating short summary [sample, depth, coverage, lineage]...")
    chrm_df = pd.read_csv(chromosomes_csv, index_col=0)
    pang_df = pd.read_csv(pango_csv, index_col=0)
    # get short summary csv
    short_summary_df = load_short_summary_df(chrm_df, pang_df)
    short_summary_df.to_csv(outdir+'/short_summary.csv')
    print(f'  > {outdir}short_summary.csv')
    return lineage_df, short_summary_df

#def get_short_summary(pango_csv, chromosomes_csv, ):
    # load chromossomes and mix with pango strain
