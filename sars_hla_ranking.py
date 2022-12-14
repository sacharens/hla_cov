"""
##############################################################################################
###################################### THE GRAND FINALLY T CELL SCORE ########################
##############################################################################################

this scrip will give us the top peptides that are suspected to have mutation that escape hla binding,
or TCR recognition.




by: Sinai Sacharen

"""

# ---------------------""" imports  & paths"""------------------------- <editor-fold>
import concurrent.futures

from Tool_File import *
import gc
# import statsmodels

starttime = time.time()

path_with_data = '/home/sacharen/Documents/hlaeecovrege/'
path_with_gisaid_data = '/home/sacharen/Documents/hla_ee/'
path_allele_iedb_delf = '/home/sacharen/Documents/sars_hla_ranking_data/'
path_to_save = '/home/sacharen/Documents/'

if os.path.exists(path_to_save) and os.path.exists(path_with_data):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')
    exit()

#
# ------------------------------------------------------------------------ </editor-fold>

# --------------------------------"""files to import"""---------------- <editor-fold>
# main_df = pd.read_csv(path_with_data + 'main_df.csv', index_col=0)
ref_df = pd.read_excel('/home/sacharen/Documents/find_ref/new_ref_df.xlsx')
ref_df.set_index(['GISAID name'], inplace=True)
monthly_df = pd.read_csv(path_with_gisaid_data + 'monthly_df_all_prot.csv', sep=',')

allele_data_omst = pd.read_csv(path_allele_iedb_delf + "Luo_SD_1_inferred_HLAallele_summary.csv", index_col=0,
                               sep='\s+', header=1)
allele_data_omst = allele_data_omst.filter(regex='^[A,B]/*', axis=0)
# allele_data_omst = allele_data_omst.filter(regex='G$',axis=0)
allele_data_omst = allele_data_omst[allele_data_omst['Freq'] > 0.01]
allele_data_omst.reset_index(inplace=True)
# sars t cell
iedb_t_cell = pd.read_csv(path_allele_iedb_delf + 'sars_t_cell.csv', skiprows=1)
t_cell_set = set(iedb_t_cell['Description'])
# sars hla
sars_hla = pd.read_csv(path_allele_iedb_delf + 'sars_hla.csv', skiprows=1)
sars_hla = set(sars_hla['Description'])

monthly_df_melt = pd.melt(monthly_df, id_vars=['Unnamed: 0', 'date', 'protein'],
                          value_vars=['A', 'R', 'N', 'D', 'C', 'Q',
                                      'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-'])

entropy_df_path = '/home/sacharen/Documents/hla_ee/'

Hydropathicity = {'Ala': 1.800, 'Arg': -4.500, 'Asn': -3.500, 'Asp': -3.500, 'Cys': 2.500, 'Gln': -3.500,
                  'Glu': -3.500, 'Gly': -0.400, 'His': -3.200, 'Ile': 4.500, 'Leu': 3.800, 'Lys': -3.900,
                  'Met': 1.900, 'Phe': 2.800, 'Pro': -1.600, 'Ser': -0.800, 'Thr': -0.700, 'Trp': -0.900,
                  'Tyr': -1.300, 'Val': 4.200}
Bulkiness = {'Ala': 11.500, 'Arg': 14.280, 'Asn': 12.820, 'Asp': 11.680, 'Cys': 13.460, 'Gln': 14.450,
             'Glu': 13.570, 'Gly': 3.400, 'His': 13.690, 'Ile': 21.400, 'Leu': 21.400, 'Lys': 15.710, 'Met': 16.250,
             'Phe': 19.800, 'Pro': 17.430, 'Ser': 9.470, 'Thr': 15.770, 'Trp': 21.670, 'Tyr': 18.030, 'Val': 21.570}
aa_three_to_one = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', 'G': 'Gly',
                   'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser',
                   'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'}


# ------------------------------------------------------------------------ </editor-fold>

# --------------------------------""" Functions """----------------------- <editor-fold>


def count_number_of_binders(score_list):
    score_list = list(score_list)
    count = 0
    # print(score_list)
    for score in score_list:
        if score >= 2:
            count += 1
    return count


def calc_pepe_mean_entropy(peptide_str, pf):
    """
    here we get a peptide and the spikes entropy per position as a global parameter
    so the mean entropy of a peptide can be calculated
    :param peptide_str: a peptide
    :return: the mean peptide Entropy score
    """
    global ref_df
    global entropy_df_path
    entropy_df = pd.read_excel(entropy_df_path+ pf + '_entropy_matrix.xlsx')

    prot_seq = ref_df.loc[pf, 'sequence']
    peptide_start = find_peptide_position(peptide_str, ref_seq=prot_seq)
    peptide_end = peptide_start + len(peptide_str)
    # print(peptide,peptide_start)
    pep_entropy = entropy_df.loc[peptide_start:peptide_end, 'entropy'].mean()
    return pep_entropy


def find_number_of_monthly_occurrences(main_df_row):
    """
    this function will take a row from the main df and will return the number of months
    that the peptide has been seen more than a set threshold
    and the mean frequency throughout the months
    and the month with the highest frequency
    :param main_df_row:
    :return: row with the new data
    """
    global monthly_df_melt

    position = main_df_row['mutations position']
    mut_aa = main_df_row['mutations aa']
    protein_f = main_df_row['protein']
    # slice the df
    sliced_df = monthly_df_melt[
        (monthly_df_melt['Unnamed: 0'] == position) & (monthly_df_melt['variable'] == mut_aa) & (
                monthly_df_melt['protein'] == protein_f) & (monthly_df_melt['value'] > 0)]
    number_of_months = list(sliced_df['value'])
    if len(number_of_months) == 0:
        main_df_row['max month for mutation'] = 'no max month'
        main_df_row['Mutation monthly AVG'] = 'no data'
        main_df_row['Number of months gt zero'] = 'no data'
        main_df_row['Number of months gt 0.01'] = 'no data'
        main_df_row['Number of months gt 0.1'] = 'no data'
    else:
        idx = sliced_df['value'].idxmax()
        max_mont = sliced_df.loc[idx, 'date']
        main_df_row['max month for mutation'] = max_mont

        month_mean = np.mean(number_of_months)
        main_df_row['Mutation monthly AVG'] = month_mean
        main_df_row['Number of months gt zero'] = len(number_of_months)
        main_df_row['Number of months gt 0.01'] = sliced_df[sliced_df['value'] > 0.01].shape[0]
        main_df_row['Number of months gt 0.1'] = sliced_df[sliced_df['value'] > 0.1].shape[0]

    return main_df_row


def phisical_change(row_data, phis='B', id=False):
    global ref_df
    global Hydropathicity
    global Bulkiness
    global aa_three_to_one
    pos = row_data[1]
    mut_aa = row_data[2]
    prot_name = row_data[0]
    # print(mut_aa)
    prot_seq = ref_df.loc[prot_name, 'sequence']
    wt_aa = prot_seq[pos - 1]
    if id == True:
        return wt_aa + str(pos) + mut_aa
    if mut_aa == '-':
        return np.nan
    else:
        if phis == "B":
            delta_phis = Bulkiness[aa_three_to_one[mut_aa]] - Bulkiness[aa_three_to_one[wt_aa]]
            return delta_phis
        if phis == "H":
            delta_phis = Hydropathicity[aa_three_to_one[mut_aa]] - Hydropathicity[aa_three_to_one[wt_aa]]
            return delta_phis


def get_allele_type(allel_str):
    # a = 'HLA-A02:02'
    return allel_str[4:5]


def row_anchor_paprams(row_data):
    mut_pos = row_data[0]
    statr_pos = row_data[2]
    end_pos = row_data[3]
    peptide = row_data[1]
    # if mut_pos == end_pos or mut_pos == statr_pos + 1 or mut_pos == statr_pos + 3 or mut_pos == statr_pos + 7:
    #     return 'P'+str(mut_pos-statr_pos+1)
    # if mut_pos == statr_pos:
    #     return 1
    # return peptide[2:-1]
    if mut_pos == statr_pos:
        return 'P1'
    else:
        return 'P' + str(mut_pos - statr_pos + 1)


def get_allele_g_type(allel_str):
    # a = 'HLA-A02:02'
    allel_str = allel_str[:7]
    allel_str = allel_str.replace('*', '')
    return 'HLA-' + allel_str


def is_upper_pep_in_genome(peptide):
    """
    this finction takes in  a peptide and aline for a file contanig a sequnce
    it returns true if positions 3-8 are in the file
    :param peptide:
    :return: bool
    """
    global path_to_save
    command = 'grep -c "[A-Z][A-Z]' + peptide[2:8] + '[A-Z]" ' + path_allele_iedb_delf + 'gencode_human.fa'
    try:
        subprocess.check_output(command, shell=True)
        # print('True')
        return True
    except subprocess.CalledProcessError:
        # print('Flse')
        return False


def parallelize_dataframe(df, func, n_cores=4):
    """
    This function takes a pandas DataFrame and applies a function to it in parallel using multiprocessing.
    using the thread pool executor from concurrent.futures library
    :param df: the df to apply the function on
    :param func: the function to apply
    :param n_cores: the number of cores to use
    :return: df with the function applied on it
    """
    df_split = np.array_split(df, n_cores * 2)
    with concurrent.futures.ThreadPoolExecutor() as executor:
        df_list = executor.map(lambda x: x.copy(deep=True).apply(func, axis=1), df_split)
    df_list = [df for df in df_list]
    print(len(df_list))

    df = pd.concat(df_list)
    return df


# ------------------------------------------------------------------------ </editor-fold>

# ---------------------MAIN----------------------------------------------- <editor-fold>



# --------------------------------------------------------------------------- </editor-fold>

# --------------------------------mia function------------------------------- <editor-fold>

def main_hla_ranking(main_df_f):
    """
    this function will addd all the information nedded for the ranking df
    :param main_df_f: a pandas df with all the data
    :return: main_df with additions
    """
    global t_cell_set
    global sars_hla
    global monthly_df_melt
    global allele_high_dict
    # sars t cell
    print('here we go....')

    main_df_f['IEDB T cell'] = main_df_f['Peptide'].apply(lambda x: True if x in t_cell_set else False)
    # sars hla

    main_df_f['IEDB SARS HLA'] = main_df_f['Peptide'].apply(lambda x: True if x in sars_hla else False)

    sup_list = ['HLA-A01:01', 'HLA-A02:01', 'HLA-A03:01', 'HLA-A24:02', 'HLA-A26:01', 'HLA-B07:02', 'HLA-B08:01',
                'HLA-B27:05', 'HLA-B40:01', 'HLA-B58:01', 'HLA-B15:01']

    wt_hla_list = []
    mut_hla_list = []
    for col in main_df_f.columns:
        if col.startswith('HLA-'):
            wt_hla_list.append(col)
        if col.startswith('MUT HLA-'):
            mut_hla_list.append(col)
    mut_sup_list = []
    for sup in sup_list:
        mut_sup_list.append('MUT ' + sup)

    main_df_f['Allele sup mean'] = main_df_f[sup_list].mean(axis=1)
    main_df_f['Allele sup MUT mean'] = main_df_f[mut_sup_list].mean(axis=1)
    main_df_f['Allele sup delta'] = main_df_f['Allele sup MUT mean'] - main_df_f['Allele sup mean']

    ### number of sup hlas that had a score grater than 2

    main_df_f['Number of SUP WT Binders'] = main_df_f[sup_list].apply(count_number_of_binders, axis=1)
    main_df_f['Number of SUP MUT Binders'] = main_df_f[mut_sup_list].apply(count_number_of_binders, axis=1)

    ### number of hlas that had a score grater than 2
    main_df_f['Number of WT Binders'] = main_df_f[wt_hla_list].apply(count_number_of_binders, axis=1)
    main_df_f['Number of MUT Binders'] = main_df_f[mut_hla_list].apply(count_number_of_binders, axis=1)
    print('finished counting binders')
    ### adding peptide mean entropy no bias
    # main_df_f['Peptide mean entropy'] = main_df_f.apply(lambda x: calc_pepe_mean_entropy(x['Peptide'], x['protein']),
    #                                                     axis=1)
    # print('finished entropy')
    ### adding monthly data

    ##### get max month #######
    main_df_f = main_df_f.apply(find_number_of_monthly_occurrences, axis=1)
    # main_df_f = parallelize_dataframe(df=main_df_f, func=find_number_of_monthly_occurrences, n_cores=32)

    print('finished the monthly part ......')

    main_df_f['Hydropathicity Delta'] = main_df_f[['protein', 'mutations position', 'mutations aa']].apply(
        phisical_change,
        phis='H', axis=1)
    main_df_f['Bulkiness Delta'] = main_df_f[['protein', 'mutations position', 'mutations aa']].apply(phisical_change,
                                                                                                      phis='B',
                                                                                                      axis=1)

    main_df_f['anchor state'] = main_df_f[['mutations position', 'Peptide', 'start pep', 'end pep']].apply(
        row_anchor_paprams,
        axis=1)

    print('finished the resi part .....')
    main_df_f['mutation ID'] = main_df_f[['protein', 'mutations position', 'mutations aa']].apply(phisical_change,
                                                                                                  id=True,
                                                                                                  axis=1)

    ################ self for upper resi #############
    main_df_f['self WT'] = main_df_f['Peptide'].apply(is_upper_pep_in_genome)
    main_df_f['self MUT'] = main_df_f['MUT peptides'].apply(is_upper_pep_in_genome)
    print('finished the self part .....')
    # main_df_f['Delta binders'] = main_df_f['Number of MUT Binders'] - main_df_f['Number of WT Binders']
    # no_bind_df = main_df_f[main_df_f['Number of WT Binders'] == 0].copy()
    # no_bind_df['nothing to something'] = no_bind_df['Number of MUT Binders'].apply(lambda x: True if x > 0 else False)
    #
    # bind_df = main_df_f[main_df_f['Number of MUT Binders'] == 0].copy()
    # bind_df['something to nothing'] = bind_df['Number of WT Binders'].apply(lambda x: True if x > 0 else False)



    # -------------------------------ESCAPE SCORE------------------------------- <editor-fold>



    def smart_delta(score_row):
        my_new_exposure = 0
        my_new_delta = 0
        my_hla_list = []
        my_hla_exposure_list = []
        global allele_high_dict
        for allele in allele_high_dict.keys():

            if score_row[allele] > 2 and score_row['MUT ' + allele] < 2:  # find if an escape happend
                my_new_delta = my_new_delta + ((score_row[allele] - score_row['MUT ' + allele])) * allele_high_dict[
                    allele]
                my_hla_list.append(allele)

            if score_row[allele] < 2 and score_row['MUT ' + allele] > 2:  # find if there was exposure
                my_new_exposure = my_new_exposure + (score_row['MUT ' + allele] - score_row[allele]) * \
                                  allele_high_dict[allele]
                my_hla_exposure_list.append(allele)

        score_row['escape score prevalent alleles'] = my_new_delta if len(my_hla_list) != 0 else None
        score_row['prevalent alleles seen escape'] = my_hla_list
        score_row['exposure score prevalent alleles'] = my_new_exposure if len(my_hla_exposure_list) != 0 else None
        score_row['prevalent alleles seen exposure'] = my_hla_exposure_list
        return score_row


    main_df_f = main_df_f.apply(smart_delta, axis=1)

    print('finished smart delta ....')
    return main_df_f



# ---------------------------MIAN------------------------------ </editor-fold>


# ------------concat all pritien dfs together to get a main df------------ <editor-fold>

text_files = [f for f in os.listdir(path_with_data) if f.endswith('final_df.xlsx')]
all_prot_df = pd.DataFrame()
count = 0
for file in text_files:
    # if file[:-14] != 'E':
    #     continue
    temp_df = pd.read_excel(path_with_data + file, index_col=0)
    temp_df['protein'] = temp_df['Peptide'].apply(lambda x: file[:-14])
    print(file[:-14])

    # shift column 'Name' to first position
    first_column = temp_df.pop('protein')
    temp_df.insert(0, 'protein', first_column)
    if count == 0:
        all_prot_df = temp_df.copy()
        count += 1
        continue
    all_prot_df = pd.concat([all_prot_df, temp_df])

all_prot_df['Delta'] = all_prot_df['MUT mean'] - all_prot_df['mean']
all_prot_df.to_csv(path_to_save + 'main_df.csv')
main_df = all_prot_df

main_df = main_df[main_df['start pep'] != -1]
main_df = main_df[main_df['mutations aa'] != '-']
main_df = main_df[main_df['mutations aa'] != 'C']



allele_data_omst['Allele'] = allele_data_omst['Allele'].apply(get_allele_g_type)

# allele dta to dict
allele_high_dict = {}
for idx in allele_data_omst.index:
    key = allele_data_omst.at[idx, 'Allele']
    allele_high_dict[key] = allele_data_omst.at[idx, 'Freq']

# main_df = parallelize_dataframe(df=main_df, func=main_hla_ranking, n_cores=32)
df_split = np.array_split(main_df, 32)


with concurrent.futures.ProcessPoolExecutor() as executor:
    df_list = executor.map(main_hla_ranking, df_split)
df_list_finish = [df for df in df_list]
print(len(df_list_finish))

main_df = pd.concat(df_list_finish)
main_df.to_csv(path_to_save + 'END_GAME_12_22.csv')

path_to_nas = '/run/user/1000/gvfs/afp-volume:host=HERTZ-LAB-NAS.local,user=sinai,volume=students/Sinai/GISAID_nov_2022/'
if os.path.exists(path_to_nas):
    print('found path')
    main_df.to_csv(path_to_nas + 'END_GAME_12_22.csv')

# ----------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------ </editor-fold>
endtime = time.time()
seconds = endtime - starttime
hours = seconds // (60 * 60)
seconds %= (60 * 60)
minutes = seconds // 60
seconds %= 60
print('end of script....')
print('time:')
print("%02i:%02i:%02i" % (hours, minutes, seconds))

# plot protien
mail_function(subject='sars hla ranking ', content='done....')

