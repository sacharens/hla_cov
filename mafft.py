"""


by: Sinai Sacharen

"""
import subprocess

import pandas as pd

# ------------------------""" imports """----------------------------------- <editor-fold>

# from hla_cov.Tool_File import *
# from hla_cov.pygmail import *
from Tool_File import *

# ------------------------------------------------------------------------ </editor-fold>


# -----------------------""" paths """----------------------------------- <editor-fold>

path_with_data = '/home/sacharen/Downloads/'
path_to_save = '/home/sacharen/Documents/' + '/hla_ee/'
nas_path = '/run/user/1000/gvfs/afp-volume:host=HERTZ-LAB-NAS.local,user=sinai,volume=students/Sinai/GISAID_nov_2022/'
if os.path.exists(path_to_save) and os.path.exists(path_with_data) and os.path.exists(nas_path):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')

# ------------------------------------------------------------------------ </editor-fold>


# -------------------------------""" global variables """------------------ <editor-fold>


ref_df = pd.read_excel('/home/sacharen/Documents/find_ref/new_ref_df.xlsx')
ref_df.drop(['Unnamed: 0'], axis=1, inplace=True)
ref_df.set_index(['GISAID name'], inplace=True)


# ------------------------------------------------------------------------ </editor-fold>


# ---------------------------""" MAIN """--------------------------------- <editor-fold>
file_og_name = 'allprot1108'

cmd = 'tar -xvf ' + path_with_data + file_og_name + '.tar.xz --directory ' + path_to_save
os.system(cmd)

fasta_file_gisaid = file_og_name + '.fasta'

print(fasta_file_gisaid)
# P_list = ['N','NSP8','NSP9']
p_list = ['Spike', 'NSP10', 'NSP9', 'NSP8', 'NSP7', 'NSP6', 'NSP5', 'NSP4', 'NSP3', 'NSP2', 'NSP1', 'NS7b', 'NSP11',
          'NSP16',
          'NSP15', 'NSP14', 'NSP13', 'NSP12', 'N', 'NS8', 'NS7a', 'NS6', 'M', 'E', 'NS3']

count_file_seq = 0
human = 0
cour = 0
number_of_seq_dict = {}

with open(path_to_save + file_og_name + '/' + fasta_file_gisaid, "r", encoding='latin1') as input_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        count_file_seq += 1
        # print(record.seq)
        P = re.findall(r'\w+', record.id)[0]
        if P == 'NS9a' or P == 'NS9b' or P == 'NS9c':
            continue
        if 'Human' in record.description:
            human += 1
        else:
            continue
        ref_seq_f = ref_df.loc[P, 'sequence']
        if P not in number_of_seq_dict.keys():
            number_of_seq_dict[P] = [1, 0]
        else:
            number_of_seq_dict[P][0] += 1

        range_val = len(ref_seq_f) * 0.025
        low_limit = len(ref_seq_f) - range_val
        hight_limit = len(ref_seq_f) + range_val
        if record.seq.find('X') == -1 and low_limit < len(record.seq) < hight_limit:
            cour += 1
            number_of_seq_dict[P][1] += 1

            with open(path_to_save + P + '_no_x.fa', 'a') as handle:
                SeqIO.write(record, handle, 'fasta')

        if count_file_seq % 777777 == 0:
            print(count_file_seq)

print('here we at ..........')

file1 = open(path_to_save + "log_file.txt", "a")  # append mode
file1.write("total number of sequences:" + str(count_file_seq) + " \n")
file1.write("total number of sequences after curation :" + str(cour) + " \n")
file1.write("total number of HUMAN sequences :" + str(human) + " \n")
file1.write("--------------------------------------------------------------------------------" + '\n')
for key, val in number_of_seq_dict.items():
    file1.write("total number of sequences in protein " + key + " before curation:" + str(val[0]) + " \n")
    file1.write("total number of sequences in protein " + key + " after :" + str(val[1]) + " \n")
    file1.write("--------------------------------------------------------------------------------" + '\n' + '\n')
file1.close()

# -------------------------------------MAFFT---------------------------------- </editor-fold>


for P in p_list:
    all_seq_list = []
    input_handle = open(path_to_save + P + '_no_x' + '.fa', "r")
    for record in SeqIO.parse(input_handle, "fasta"):
        all_seq_list.append(record)

    input_handle.close()

    ##  here we eak with the fct that the file is to big ##
    third = round(len(all_seq_list) / 3)
    random.shuffle(all_seq_list)

    list_for_save_A = all_seq_list[:third]
    list_for_save_B = all_seq_list[third:third * 2]
    list_for_save_C = all_seq_list[third * 2:]

    # discarding x
    SeqIO.write(list_for_save_A, path_to_save + P + '_no_x_A' + '.fa', 'fasta')
    SeqIO.write(list_for_save_B, path_to_save + P + '_no_x_B' + '.fa', 'fasta')
    SeqIO.write(list_for_save_C, path_to_save + P + '_no_x_C' + '.fa', 'fasta')

    for f in ['A', 'B', 'C']:
        input_file = P + '_no_x_' + f + '.fa'
        output_file = path_to_save + P + "_maffet_alinged_" + f + ".fa"
        cd_c = 'cd ' + path_to_save + '; '
        command = cd_c + "mafft --6merpair --keeplength --thread -1 --anysymbol --addfragments  " + input_file + " /home/sacharen/Documents/Maffet_projet/" + P + "_fasta.fa > " + output_file
        sed_c = "'s/*//g'"
        comm = 'sed -i ' + sed_c + ' ' + path_to_save + input_file
        subprocess.run('%s' % 'cd', shell=True)

        subprocess.run('%s' % comm, shell=True)

        try:
            subprocess.run('%s' % command, shell=True)
        except subprocess.CalledProcessError:
            continue

    cat_com = 'cat ' + path_to_save + P + '_maffet_alinged_A' + '.fa ' + path_to_save + P + '_maffet_alinged_B' + '.fa ' + path_to_save + P + \
              '_maffet_alinged_C' + '.fa > ' + path_to_save + P + '_maffet_alinged.fa'

    subprocess.run('%s' % cat_com, shell=True)
    #### remove the  helper files ####
    rm1 = 'rm -v ' + path_to_save + P + '_no_x_A.fa ' + path_to_save + P + '_no_x_B.fa ' + path_to_save + P + '_no_x_C.fa'
    subprocess.run('%s' % rm1, shell=True)
    rm2 = 'rm -v ' + path_to_save + P + '_maffet_alinged_A.fa ' + path_to_save + P + '_maffet_alinged_B.fa ' + path_to_save + P \
          + '_maffet_alinged_C.fa'
    subprocess.run('%s' % rm2, shell=True)

    print(P,
          'done-----------------------------------------------------------------------------------------------------')

# -------------------------------""" get all dates """------------------ <editor-fold>

date_list = pd.date_range(start="3/3/2020", end=datetime.datetime.today()).tolist()

new_date_list = []
for i in date_list:
    month_year = str(i.month) + '-' + str(i.year)
    if month_year not in new_date_list:
        new_date_list.append(month_year)

date_list = pd.to_datetime(new_date_list).strftime("%Y-%b").tolist()
date_list.pop(-1)


# ------------------------------------------------------------------------ </editor-fold>


def all_mutation_info_for_protein(P_f, path_to_save_f):

    print('started '+P_f+'..............')
    global ref_df
    global date_list
    ref_sequence = ref_df.loc[P_f, 'sequence']

    monthly_df = pd.DataFrame()
    monthly_count_df = pd.DataFrame()
    count = 0
    count_c = 0
    count1 = 0

    list_of_recoreds = []


    with open(path_to_save_f + P_f + '_maffet_alinged.fa', "r", encoding='latin1') as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            # print(record)
            count1 += 1
            list_of_recoreds.append(record)

    if count1 != len(list_of_recoreds):
        print('number of seq donot add up')
        # exit()
    ####################################################  PFM and entropy
    p, c = simple_list_PFM(list_of_recoreds, ref_sequence=ref_sequence)
    p.to_excel(path_to_save_f + P_f + '_PFM.xlsx')
    c.to_excel(path_to_save_f + P_f + '_COUNT.xlsx')
    try:
        p.to_excel(nas_path + P_f + '_PFM.xlsx')
    except OSError:
        print('didnot save monthly pfm')
    simple_entropy_function(path_to_save=path_to_save_f, save_file_name=P_f + '_entropy_matrix', pfm_df=p)
    print(P_f)
    print('PFM done!')

    # -----------------------------------entropy no bias---------------------- <editor-fold>

    # this part we creat all the monthly PFM and count matrix and-
    # the entropy with no bias, whic meens  a random sampel of equal size for evary month is taken and then
    # a pfm is calculated

    list_for_save = []
    protien_dict = {}
    for month in date_list:
        protien_dict[month] = []

    for rec in list_of_recoreds:
        try:
            dt = find_date(rec.description)
        except NameError:
            continue
        # print(type(dt))
        NoneType = type(None)
        if type(dt) != NoneType:
            dt = dt.strftime("%Y-%b")
            if dt in date_list:
                protien_dict[dt].append(rec)

    recored_temp_list = []
    for key in protien_dict.keys():
        recored_temp_list = protien_dict[key]

        ##### this part is for monthly count ####
        p, c = simple_list_PFM(recored_temp_list, ref_sequence=ref_sequence)
        c['date'] = c['A'].apply(lambda x: key)
        c['protein'] = c['A'].apply(lambda x: P_f)
        if count_c == 0:
            count_c += 1
            monthly_count_df = c.copy()
        else:
            monthly_count_df = pd.concat([monthly_count_df, c])

        ##### this part is for monthly PFM ####
        p, c = simple_list_PFM(recored_temp_list, ref_sequence=ref_sequence)
        p['date'] = p['A'].apply(lambda x: key)
        p['protein'] = p['A'].apply(lambda x: P_f)
        if count == 0:
            count += 1
            monthly_df = p.copy()
        else:
            monthly_df = pd.concat([monthly_df, p])

        random.shuffle(recored_temp_list)
        list_for_save = list_for_save + recored_temp_list[:3901]

    monthly_df.to_csv(path_to_save_f + P_f + '_monthly_pfm.csv')
    monthly_count_df.to_csv(path_to_save_f + P_f + '_monthly_count.csv')
    print(P_f, 'monthly PFM done!!!!!!!!!!!!!')

    ####################################################  PFM and entropy no bias



    if len(list_for_save) * len(date_list) != 3900 * len(date_list):
        print('number of seq donot add up in rntropy no bias')
        print(len(list_for_save) * len(date_list), 3900 * len(date_list))
        # exit()

    p, c = simple_list_PFM(list_for_save, ref_sequence=ref_sequence)
    p.to_excel(path_to_save_f + P_f + '_entropy no bias PFM.xlsx')
    c.to_excel(path_to_save_f + P_f + '_entropy no bias COUNT.xlsx')
    simple_entropy_function(path_to_save=path_to_save_f, save_file_name=P_f + '_entropy no bias entropy_matrix', pfm_df=p)

    print(P_f, 'entropy no bias PFM done!')


##################################   threding ##################################

# list_save_path = [path_to_save] * len(p_list)


# with concurrent.futures.ProcessPoolExecutor(max_workers=len(p_list)) as executor:
#     executor.map(all_mutation_info_for_protein,p_list, list_save_path)

for prot in p_list:
    all_mutation_info_for_protein(P_f=prot, path_to_save_f=path_to_save)


# concat all monthlay freq and count dfs

all_monthly_df = pd.DataFrame()
all_monthly_count_df = pd.DataFrame()
count = 0
for P in p_list:
    p_all = pd.read_csv(path_to_save + P + '_monthly_pfm.csv', sep=',', index_col=0)
    c_all = pd.read_csv(path_to_save + P + '_monthly_count.csv', sep=',', index_col=0)

    if count == 0:
        count += 1
        all_monthly_df = p_all.copy()
        all_monthly_count_df = c_all.copy()
    else:
        all_monthly_df = pd.concat([all_monthly_df, p_all])
        all_monthly_count_df = pd.concat([all_monthly_count_df, c_all])

all_monthly_df.to_csv(path_to_save + 'monthly_df_all_prot.csv')
all_monthly_count_df.to_csv(path_to_save + 'monthly_count_df_all_prot.csv')
try:
    all_monthly_df.to_csv(nas_path + 'monthly_count_df_all_prot.csv')
except OSError:
    print('did not save to NAS monthly df')
# ------------------------------------------------------------------------ </editor-fold>
mail_function(subject='gisaid pipline', content='done....')
