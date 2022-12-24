""" imports """
# ------------------------------------------------------------------------ <editor-fold>
import subprocess
import os
import time
import pandas as pd
import math
from Tool_File import *

starttime = time.time()

# ------------------------------------------------------------------------ </editor-fold>

""" paths """
# ------------------------------------------------------------------------ <editor-fold>

path_to_tool = '/home/sacharen/netMHCpan-4.1/'
path_with_data = '/home/sacharen/Documents/hlaproject/'
path_to_save = '/home/sacharen/Documents/hlaproject/'
path_to_save = path_with_data
if os.path.exists(path_to_tool) and os.path.exists(path_to_save) and os.path.exists(path_with_data):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')
    exit()


# ------------------------------------------------------------------------ </editor-fold>

# ------------------------------FUNctIONs--------------------------------- <editor-fold>

def connect_patient_to_pep(peptide_f):
    """
    this function takes in a peptide and finds out were it is lokated in the protien
    and then it can find the ralevant mutation and petient that are connected to it
    :param peptide_f: from the df_f
    global: orf1ab
    :return: the patient and the mutation
    """
    global orf1ab_seq
    global patiant_dict
    return_list = []

    pep_start = find_peptide_position(peptide=peptide_f, ref_seq=orf1ab_seq)
    pep_end = pep_start + len(peptide_f)
    for key, val in patiant_dict.items():
        for mut in val:
            mut_pos = int(mut[1:-1])
            if mut_pos in range(pep_start, pep_end + 1):
                return_list.append({key: mut})

    return return_list


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


def read_in_allels(prealent_allel_file_and_path):
    """
    tis function will read in the file and return all the allelels in it in the following format:
    'HLA-B40:01,HLA-B58:01,HLA-B15:01'
    so the netHmcpan predicter can read it in
    :param prealent_allel_file_and_path:
    :return: string
    """
    with open(prealent_allel_file_and_path) as file_in:
        lines = ''
        for line in file_in:
            temp_line = line.rstrip() + ','
            if len(temp_line) < 2:
                continue
            else:
                lines += temp_line
    return lines[:-1]


# ------------------------------------------------------------------------ </editor-fold>


# ------------------------------------------------------------------------ </editor-fold>

""" global variables """
# ------------------------------------------------------------------------ <editor-fold>
prevalent_allel_file = path_with_data + 'prevalent_alleles'
# full_alele_file = pd.read_csv(path_with_data + 'all_netMHCpan_alleles.txt', sep='\s+', header=None)
path_with_all_protien_fa = '/home/sacharen/Documents/find_ref/sars_protien_fasta/'
ref_df = pd.read_excel('/home/sacharen/Documents/find_ref/new_ref_df.xlsx')



# ------------------------------------------------------------------------ </editor-fold>

# ------------------------------params------------------------------------- <editor-fold>
# ------------------------------------------------------------------------ </editor-fold>
#
# # ------------------------------MAIN for all mhc alelles------------------- <editor-fold>

#
# # devide the list into 90 alleles because the predictor cant take in more than 80 alleles at once
# list_of_all_the_alleles = list(full_alele_file[0])
#
#
# def chunks(lst, n):
#     """Yield successive n-sized chunks from lst."""
#     for i in range(0, len(lst), n):
#         yield lst[i:i + n]
#
#
# list_of_list = list(chunks(list_of_all_the_alleles,80))
#
# for p in ['NSP10', 'NSP9', 'NSP8', 'NSP7', 'NSP6', 'NSP5', 'NSP4', 'NSP3', 'NSP2', 'NSP1', 'NS7b', 'NSP11', 'NSP16',
#           'NSP15', 'NSP14'
#     , 'NSP13', 'NSP12', 'N', 'NS8', 'NS7a', 'NS6', 'M', 'E', 'NS3', 'Spike']:
#
#     print('processing protein '+p+'...........................')
#     file_and_path_with_fasta = path_with_all_protien_fa + p + "_fasta.fa"
#
#     # creat temp xlsx file
#     # with open(path_with_data + 'temp.xlsx', 'a') as fp:
#     #     pass
#
#     count = 0
#     for alleles in list_of_list:
#         print('loop number: '+str(count))
#         HLA_str = ","+','.join(alleles)+","
#         open(path_with_data + 'temp.xlsx', 'a').close()
#
#         temp_xlsx_file_and_path_f = path_with_data + 'temp.xlsx'
#
#         # HLA_str = read_in_allels(prevalent_allel_file)
#
#         # peptide_length = "8,9,10,11"  # -l
#         peptide_length = "9"  # -l
#
#         command = path_to_tool + "netMHCpan " + " -a " + HLA_str + " -f " + file_and_path_with_fasta + " -l " + "9" + " -xls " + " -xlsfile " + temp_xlsx_file_and_path_f
#
#         subprocess.check_output(command, shell=True)
#
#         df_f = pd.read_csv(temp_xlsx_file_and_path_f, sep='\t', header=[0, 1])
#
#         df_f = proces_out_MHC_xls(pro_df=df_f)
#         df_f.to_csv(path_to_save + p + 'hla_coverage_'+str(count)+'.csv')
#         # delete temp xlsx file
#         os.remove(temp_xlsx_file_and_path_f)
#         count += 1
# # ------------------------------------------------------------------------ </editor-fold>


# ------------------------------MAIN for all prevalent alleles------------- <editor-fold>
for p in ['NSP10', 'NSP9', 'NSP8', 'NSP7', 'NSP6', 'NSP5', 'NSP4', 'NSP3', 'NSP2', 'NSP1', 'NS7b', 'NSP11', 'NSP16',
          'NSP15', 'NSP14'
    , 'NSP13', 'NSP12', 'N', 'NS8', 'NS7a', 'NS6', 'M', 'E', 'NS3', 'Spike']:

    print('processing protein '+p+'...........................')
    file_and_path_with_fasta = path_with_all_protien_fa + p + "_fasta.fa"

    # creat temp xlsx file
    # with open(path_with_data + 'temp.xlsx', 'a') as fp:
    #     pass


    open(path_with_data + 'temp.xlsx', 'a').close()

    temp_xlsx_file_and_path_f = path_with_data + 'temp.xlsx'

    HLA_str = read_in_allels(prevalent_allel_file)

    # peptide_length = "8,9,10,11"  # -l
    peptide_length = "9"  # -l

    command = path_to_tool + "netMHCpan " + " -a " + HLA_str + " -f " + file_and_path_with_fasta + " -l " + "9" + " -xls " + " -xlsfile " + temp_xlsx_file_and_path_f

    subprocess.run(command, shell=True)

    df_f = pd.read_csv(temp_xlsx_file_and_path_f, sep='\t', header=[0, 1])

    df_f = proces_out_MHC_xls(pro_df=df_f)
    df_f.to_csv(path_to_save + p + 'hla_coverage.csv')
    # delete temp xlsx file
    os.remove(temp_xlsx_file_and_path_f)
# ------------------------------------------------------------------------ </editor-fold>


endtime = time.time()
seconds = endtime-starttime
hours = seconds // (60*60)
seconds %= (60*60)
minutes = seconds // 60
seconds %= 60
print( 'end of script....')
print('time:')
print("%02i:%02i:%02i" % (hours, minutes, seconds))