""" imports """
# ------------------------------------------------------------------------ <editor-fold>

from Tool_File import *
from utiles import *
import math
starttime = time.time()
# ------------------------------------------------------------------------ </editor-fold>

""" paths """
# ------------------------------------------------------------------------ <editor-fold>

path_to_tool = '/home/sacharen/netMHCpan-4.1/'
path_with_data = '/home/sacharen/Documents/hlaeecovrege/'
path_to_save = path_with_data
if os.path.exists(path_to_tool) and os.path.exists(path_to_save) and os.path.exists(path_with_data):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')


# ------------------------------------------------------------------------ </editor-fold>

# ------------------------------functions--------------------------------- <editor-fold>


def scale_function(number):
    e = math.e
    return 1 / (e ** (number - e))


def proces_out_MHC_xls(pro_df):
    """
    this function takes in the output xlsx file from net MHCpan and
    cleans it up so ut can be merged into the main df
    :param pro_df:
    :return: clean df
    """
    new_col_dict = {}

    current_hla = ''
    for col in pro_df.columns:
        print(col)
        if col[0].startswith('HLA'):
            current_hla = col[0]
        if current_hla == '':
            new_col_dict[col] = col[1]
        if col[1] == 'EL_Rank':
            new_col_dict[col] = current_hla

    pro_df.columns = pro_df.columns.to_flat_index()
    pro_df.rename(columns=new_col_dict, inplace=True)
    for col in pro_df.columns:
        if type(col) is tuple:
            pro_df.drop(col, axis=1, inplace=True)

    allele_cols_scale_part = []
    for dcs in pro_df.columns:
        if 'HLA-' in dcs:
            allele_cols_scale_part.append(dcs)
    print(pro_df)
    pro_df[allele_cols_scale_part] = pro_df[allele_cols_scale_part].apply(scale_function)
    pro_df['mean'] = pro_df[allele_cols_scale_part].mean(axis=1)

    return pro_df


def is_not_og_aa(aa_id_f, ref_s):
    pos = int(aa_id_f[:-1])
    aa = aa_id_f[-1:]
    if ref_s[pos - 1] == aa:
        return False
    else:
        return True


# ------------------------------------------------------------------------ </editor-fold>

# ------------------------------input files----------------------------- <editor-fold>
prevalent_allel_file = path_with_data+'prevalent_alleles'
file_with_mut_peptides = path_with_data+'temp_mut_prp_file.txt'
temp_xlsx_file_and_path_f = path_with_data+'temp_file_netMHC.xlsx'
ref_df = pd.read_excel('/home/sacharen/Documents/find_ref/new_ref_df.xlsx')
ref_df.set_index(['GISAID name'], inplace=True)

# ------------------------------------------------------------------------ </editor-fold>

# ------------------------------params------------------------------------ <editor-fold>
for P in ['NSP10','NSP9','NSP8','NSP7','NSP6','NSP5','NSP3','NSP4','NSP2','NSP1','NS7b','NSP11','NSP16','NSP15','NSP14','NSP13','NSP12','N','NS8','NS7a','NS6','M','E','NS3','Spike']:
# for P in ['E']:
    print('prtotien: '+P+'   .................................')
    # P = 'M'
    prot_seq= ref_df.loc[P, 'sequence']


    # ------------------------------------------------------------------------ </editor-fold>

    # -------------------------------take all 9 mers from hla coverage file --- <editor-fold>

    df_f = pd.read_csv(path_with_data+P+'hla_coverage.csv',index_col=0)
    df_f['len pep'] = df_f['Peptide'].apply(lambda x:len(x))
    df_f = df_f[df_f['len pep']==9]
    df_f.drop('len pep',axis=1,inplace=True)
    unique_pep_list = list(df_f['Peptide'])

    # ------------------------------------------------------------------------ </editor-fold>


    # ----------------------- a mut list from the pfm no bias file------ <editor-fold>
    mut_pdf_df = pd.read_excel(path_with_data + P+'_entropy no bias PFM.xlsx')
    mut_pdf_df = pd.melt(mut_pdf_df, id_vars=['Unnamed: 0'],
                         value_vars=['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T',
                                     'W', 'Y', 'V', '-'])

    # mut_pdf_df = mut_pdf_df[mut_pdf_df['value'] > 0.0001]  ###### the mut cut off######
    mut_pdf_df['mus_id'] = mut_pdf_df['Unnamed: 0'].astype(str) + mut_pdf_df['variable']

    mut_pdf_df['for list'] = mut_pdf_df['mus_id'].apply(is_not_og_aa, ref_s=prot_seq)
    mut_pdf_df = mut_pdf_df[mut_pdf_df['for list']]

    mut_list = list(mut_pdf_df['mus_id'])

    # ----------------------------------------------------------------------------- </editor-fold>

    # creat the main df with all data  ----------------------------------------------------------------------

    main_df = pd.DataFrame(columns=['mutations position', 'mutations aa', 'Peptide', 'start pep', 'end pep'])

    main_idx = 0
    for mut in mut_list:
        mut_position = int(mut[:-1])
        mutation_aa = mut[-1:]
        for pep in unique_pep_list:
            if is_mutation_in_peptide(mut_pos_f=mut_position, pep_f=pep, ref_sequence_f=prot_seq):
                start_pos = find_peptide_position(pep, ref_seq=prot_seq)
                end_pos = start_pos + len(pep) - 1
                main_df.at[main_idx, 'mutations position'] = mut_position
                main_df.at[main_idx, 'mutations aa'] = mutation_aa
                main_df.at[main_idx, 'Peptide'] = pep
                main_df.at[main_idx, 'start pep'] = start_pos
                main_df.at[main_idx, 'end pep'] = end_pos
                main_idx += 1
    main_df = pd.merge(main_df, df_f, how="left", on=["Peptide"])

    # creat all the mutated peptides in the merged df ----------------------------------------------------------------------
    main_df['MUT peptides'] = main_df.apply(
        lambda x: creat_the_mutated_peptide(x['Peptide'], x['mutations position'], x['mutations aa'],
                                            ref_sequence_f=prot_seq), axis=1)

    up_dated_mut_list = list(main_df['MUT peptides'])
    creat_pep_file_from_list(pep_list_f=up_dated_mut_list, file_path_and_name_f=file_with_mut_peptides)

    HLA_str = read_in_allels(prevalent_allel_file)
    command = path_to_tool + "netMHCpan " + " -a " + HLA_str + " -p " + file_with_mut_peptides + " -xls " + " -xlsfile " + temp_xlsx_file_and_path_f

    subprocess.check_output('%s' % command, shell=True)

    mut_predict_df = pd.read_csv(temp_xlsx_file_and_path_f, sep='\t', header=[0, 1])

    mut_predict_df = proces_out_MHC_xls(pro_df=mut_predict_df)
    dict_mut_name = {}
    for col in mut_predict_df.columns:
        if 'HLA' in col:
            dict_mut_name[col] = 'MUT ' + col

    dict_mut_name['Peptide'] = "MUT peptides"
    dict_mut_name['mean'] = "MUT mean"

    mut_predict_df.rename(columns=dict_mut_name, inplace=True)

    mut_predict_df.to_excel(path_to_save + P +'_mut_df.xlsx')

    # mut_predict_df = pd.read_excel(path_with_data+'mut_df.xlsx',index_col=0)
    # df_for_mean = mut_predict_df.filter(regex='HLA-[A-Z]\d+:\d+', axis=1)
    # mut_predict_df['MUT mean'] = df_for_mean.mean(axis=1)

    final_merged_df = pd.merge(main_df, mut_predict_df, how="left", on=["MUT peptides"])

    final_merged_df.to_excel(path_to_save + P +'_final_df.xlsx')


################################## end time ##############################
endtime = time.time()
seconds = endtime-starttime
hours = seconds // (60*60)
seconds %= (60*60)
minutes = seconds // 60
seconds %= 60
print( 'end of script....')
print('time:')
print("%02i:%02i:%02i" % (hours, minutes, seconds))