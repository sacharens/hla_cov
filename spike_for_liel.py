import time

from Tool_File import *

path_with_data = '/home/sacharen/Documents/' + '/hla_ee/'
path_to_save = '/home/sacharen/Documents/' + '/spike_11_22/'
nas_path = '/run/user/1000/gvfs/afp-volume:host=HERTZ-LAB-NAS.local,user=sinai,volume=students/Sinai/spike_11_22/'
if os.path.exists(path_to_save) and os.path.exists(path_with_data) and os.path.exists(nas_path):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')

########################
### file name       ###
#######################
fasta_file_gisaid = 'Spike_maffet_alinged.fa'

ref_df = pd.read_excel('/home/sacharen/Documents/find_ref/new_ref_df.xlsx')
ref_df.drop(['Unnamed: 0'], axis=1, inplace=True)
ref_df.set_index(['protein'], inplace=True)
ref_sequence = ref_df.loc['surface glycoprotein', 'sequence']
#
# clean_the_non_common(file_and_path=path_to_save + fasta_file_gisaid, common_length=len(ref_sequence) + 1
#                      , save_file_name_and_path=path_to_save + 'only_common.fa')

# when we will wont it to be up to date
date_list = pd.date_range(start="12/12/2019", end=datetime.datetime.today()).tolist()

# date_list = pd.date_range(start="12/12/2019",end="10/12/2020")
new_date_list = []
for i in date_list:
    month_year = str(i.month) + '-' + str(i.year)
    if month_year not in new_date_list:
        new_date_list.append(month_year)

date_list = pd.to_datetime(new_date_list).strftime("%Y-%b").tolist()
date_list.pop(-1)
date_list1 = date_list.copy()



# ----------------------------------------------------------------------------------------------------------------------
# dividing  the data by month
# ----------------------------------------------------------------------------------------------------------------------
print('dividing  the data by month ...')
start = time.time()
# creating a file for every date
for month in date_list:
    try:
        os.mkdir(path_to_save + month)
    except OSError as error:
        print(error)

for month in date_list:
    print(month)
    list_for_save = []
    input_handle = open(path_to_save + fasta_file_gisaid, "r", encoding='latin1')
    for record in SeqIO.parse(input_handle, "fasta"):
        # find the date of the sequence
        dt = find_date(record.description)
        # print(type(dt))
        NoneType = type(None)
        if type(dt) != NoneType:
            dt = dt.strftime("%Y-%b")
            # print(month, dt)
            if dt == month:
                list_for_save.append(record)
        else:
            continue

    input_handle.close()

    SeqIO.write(list_for_save, path_to_save + '/' + month + '/' + month + '_sequences.fa', 'fasta')
    # finding the number of repetition
print('finished dividing  the data by month ...')
print('time :', (time.time() - start) / 60)

# ----------------------------------------------------------------------------------------------------------------------
#  for all countries in the fasta file find contribution per month
# ----------------------------------------------------------------------------------------------------------------------

print('country contribution per month ...')


def contribution_per_month(date_list, path):
    df = pd.DataFrame(columns=date_list)
    df['full name'] = ''
    for month in date_list:
        print(month)
        input_handle = open(path + month + '/' + month + '_sequences.fa', "r", encoding='latin1')
        for record in SeqIO.parse(input_handle, "fasta"):
            location = find_seq_location(record.description)
            if location in df.index:
                df.at[location, month] = df.at[location, month] + 1
            else:
                df.at[location, month] = 1
                df.at[location, 'full name'] = record.description
            df.fillna(0, inplace=True)

        input_handle.close()
    return df


df = contribution_per_month(date_list=date_list, path=path_to_save)
df['sum'] = df.sum(axis=1)
# df.loc['sum'] = df.sum()
df.to_excel(path_to_save + 'location df.xlsx')

df = pd.read_excel(path_to_save + 'location df.xlsx', index_col=0)
####

# ----------------------------------------------------------------------------------------------------------------------
#  taking only countries that contributed more than 100 sequences
# ----------------------------------------------------------------------------------------------------------------------

df = df[df['sum'] > 10]
list_of_countries = list(df.index)
if list_of_countries[-1] == 'sum':
    list_of_countries.pop(-1)

path_and_file = path_to_save + fasta_file_gisaid


# ----------------------------------------------------------------------------------------------------------------------
#  country contributions
# ----------------------------------------------------------------------------------------------------------------------


def thread_for_country(date_list, path_and_file, country, save_path_perant, ref_sequence):
    # for country in list_of_countries:
    print('starting threading  {}....'.format(country))  # --------------------------------------------------------start

    try:
        os.mkdir(save_path_perant + country)
    except OSError as error:
        print(error)

    save_path = save_path_perant + country + '/'
    print(save_path)

    #########temp edit
    # all_files_in_folder = os.listdir(save_path)
    # if 'the_df.xlsx' in all_files_in_folder:
    #     return
    ##########

    ##########################  creating fasta file with seq from country  ################################
    list_for_save = []
    input_handle = open(path_and_file, "r", encoding='latin1')
    for record in SeqIO.parse(input_handle, "fasta"):
        # find the date of the sequence
        location = find_seq_location(record.description)
        if location == country:
            list_for_save.append(record)

    input_handle.close()
    SeqIO.write(list_for_save, save_path + country + '_sequences.fa', 'fasta')

    column_list = ['position', 'mutated aa']
    # column_list.extend(date_list)
    the_mane_df = pd.DataFrame(columns=column_list)

    for month in date_list:
        print(month)
        list_for_save = []
        input_handle = open(save_path + country + '_sequences.fa', "r", encoding='latin1')
        for record in SeqIO.parse(input_handle, "fasta"):
            # find the date of the sequence
            dt = find_date(record.description)
            # print(type(dt))
            NoneType = type(None)
            if type(dt) != NoneType:
                dt = dt.strftime("%Y-%b")
                # print(month, dt)
                if dt == month:
                    list_for_save.append(record)
            else:
                continue

        input_handle.close()
        if len(list_for_save) == 0:
            continue

        try:
            os.mkdir(save_path + month)
        except OSError as error:
            print(error)
        num_of_seq_in_month = len(list_for_save)

        ####################################################################
        ####### if the month has less than 10 sequences we dont care ######

        if num_of_seq_in_month < 10:
            continue
        ####################################################################

        SeqIO.write(list_for_save, save_path + month + '/' + month + '_sequences.fa', 'fasta')

        pfm, count_matrix = simple_PFM(file_and_path=save_path + month + '/' + month + '_sequences.fa',
                                       ref_sequence=ref_sequence)
        mutation_df = simple_mutation_finder(count_df=count_matrix, protein='spike', ref_sequence=ref_sequence)

        pfm.to_excel(save_path + month + '/' + month + '_PFM.xlsx')
        mutation_df.to_excel(save_path + month + '/' + month + '_mutation.xlsx')

        mutation_df.reset_index(inplace=True)
        mutation_df.drop(mutation_df[mutation_df['number of occurrences'] == 0].index, inplace=True)
        mutation_df['number of occurrences'] = (mutation_df['number of occurrences'] / num_of_seq_in_month) * 100
        mutation_df['mutation aa'] = mutation_df['new_index_no'].apply(lambda x: x[-1:])
        filtered_df = mutation_df[mutation_df['number of occurrences'] > 1].copy(
            deep=True)  # droping all sequense that are less the 1 percent
        num_of_seq_in_month = str(num_of_seq_in_month)
        filtered_df['amino acid position'] = filtered_df['amino acid position'].apply(lambda x: int(x))
        filtered_df.reset_index(inplace=True, drop=True)
        the_mane_df[month + ' ' + num_of_seq_in_month] = 0
        filtered_df.to_excel(save_path + month + '/' + month + '_filtered.xlsx')
        for idx in filtered_df.index:
            a = int(filtered_df.loc[idx, 'amino acid position'])
            temp_list = list(the_mane_df['position'])
            if a in temp_list:

                # if filtered_df.at[idx,'amino acid position'] in the_mane_df['position']:
                the_mane_df.reset_index(inplace=True, drop=True)
                for mane_idx in the_mane_df.index:

                    b = int(the_mane_df.loc[mane_idx, 'position'])
                    x = the_mane_df.loc[mane_idx, 'mutated aa']
                    y = filtered_df.loc[idx, 'mutation aa']

                    if x != y and b == a:
                        if temp_list.count(a) > 1:
                            continue
                        else:
                            new_mane_idx = len(the_mane_df.index) + 1
                            the_mane_df.loc[new_mane_idx, 'mutated aa'] = filtered_df.loc[idx, 'mutation aa']
                            the_mane_df.loc[new_mane_idx, 'position'] = filtered_df.loc[idx, 'amino acid position']
                            the_mane_df.loc[new_mane_idx, month + ' ' + num_of_seq_in_month] = filtered_df.loc[
                                idx, 'number of occurrences']
                            break

                    if x == y and b == a:
                        the_mane_df.loc[mane_idx, month + ' ' + num_of_seq_in_month] = filtered_df.loc[
                            idx, 'number of occurrences']
                        break

                    # if the_mane_df.at[mane_idx, 'mutated aa'] == filtered_df.at[idx, 'mutation aa'] and the_mane_df.at[mane_idx, 'position'] == filtered_df.at[idx, 'amino acid position']:
                    #     the_mane_df.at[mane_idx, month] = filtered_df.at[idx, 'number of occurrences']
                    #     break
                    # else:
                    #     if the_mane_df.at[mane_idx, 'position'] == filtered_df.at[idx, 'amino acid position']and the_mane_df.at[mane_idx, 'mutated aa'] != filtered_df.at[idx, 'mutation aa']:
                    #         new_mane_idx = len(the_mane_df.index) + 1
                    #         the_mane_df.at[new_mane_idx, 'mutated aa'] = filtered_df.at[idx, 'mutation aa']
                    #         the_mane_df.at[new_mane_idx, 'position'] = filtered_df.at[idx, 'amino acid position']
                    #         the_mane_df.at[new_mane_idx, month] = filtered_df.at[idx, 'number of occurrences']
                    #         break


            else:
                the_mane_df.reset_index(inplace=True, drop=True)
                new_mane_idx = len(the_mane_df.index) + 1
                the_mane_df.loc[new_mane_idx, 'mutated aa'] = filtered_df.loc[idx, 'mutation aa']
                the_mane_df.loc[new_mane_idx, 'position'] = filtered_df.loc[idx, 'amino acid position']
                the_mane_df.loc[new_mane_idx, month + ' ' + num_of_seq_in_month] = filtered_df.loc[
                    idx, 'number of occurrences']

            the_mane_df.fillna(0, inplace=True)
        the_mane_df.sort_values(by=['position'], inplace=True)
    the_mane_df.to_excel(save_path + 'the_df.xlsx')

    print(
        'threading {} ended....'.format(country))  # ----------------------------------------------------------------end


try:
    os.mkdir(path_to_save + 'country')
except OSError as error:
    print(error)
save_path = path_to_save + 'country' + '/'

country = list_of_countries
path_and_file = path_to_save + fasta_file_gisaid
path_and_file = [path_and_file] * len(country)
date_list = [date_list] * len(country)
ref_sequence = [ref_sequence] * len(country)
save_path_perant = [save_path] * len(country)

start = time.time()
with concurrent.futures.ThreadPoolExecutor() as executor:
    executor.map(thread_for_country, date_list, path_and_file, country, save_path_perant, ref_sequence)
print('all threads time :', (time.time() - start) / 60)

# ----------------------------------------------------------------------------------------------------------------------
#  creating the mane df with all the country's contributions
# ----------------------------------------------------------------------------------------------------------------------

print('creating the mane df with all the countrys contributions .... ')

path_to_save = path_to_save + 'country/'  # must end with a /
col = ['position', 'mutated aa', 'country'] + date_list1
all_country_df = pd.DataFrame(columns=col).astype(float)

# date_list= [ '2020-Mar', '2020-Jul', '2020-Aug','2020-Sep', '2020-Oct', '2020-Nov', '2020-Dec']
#####
# df = pd.read_excel(path_to_save + 'Iceland/the_df.xlsx', index_col=0)
####


# dictionary of number of countries per month
country_per_month_dict = {}
for _ in date_list1:
    country_per_month_dict[_] = 0

start = time.time()
all_files_in_folder = os.listdir(path_to_save)
for c in all_files_in_folder:
    df = pd.read_excel(path_to_save + c + '/the_df.xlsx', index_col=0)
    df.reset_index(inplace=True, drop=True)
    df['country'] = df['position'].apply(lambda x: c)
    # dictionary of number of countries per month
    list_c = list(df.columns)
    list_c.remove('position')
    list_c.remove('mutated aa')
    list_c.remove('country')
    temp_col_dict = {}
    for temp_col in list_c:
        temp_col_dict[temp_col] = temp_col.split(' ')[0]

    for k in temp_col_dict.values():
        country_per_month_dict[k] += 1

    df.rename(columns=temp_col_dict, inplace=True)
    # for df_idx in df.index:
    #     mutated_aa=df.at[df_idx,'mutated aa']
    #     position =df.at[df_idx,'position']
    # df = df.replace(0,None)
    all_country_df = all_country_df.append(df)
    # all_country_df = pd.concat([all_country_df, df ])
all_country_df = all_country_df[col]
print('========================================================================')

all_country_df.to_excel(save_path + 'all_country_df.xlsx')
# all_country_df.fillna(None,inplace=True)
# temp_df = all_country_df.replace(0,None)

temp_df = all_country_df.groupby(['position', 'mutated aa'], as_index=False).mean()
temp_df.to_excel(save_path + 'mean_country_df.xlsx')

n_df = all_country_df.groupby(['position', 'mutated aa'], as_index=False).sum()
col_list_n_df = list(n_df.columns)
col_list_n_df.remove('position')
col_list_n_df.remove('mutated aa')
for col in col_list_n_df:
    print(col)
    print(country_per_month_dict[col])
    n_df[col] = n_df[col] / (country_per_month_dict[col])

n_df.to_excel(save_path + 'normalised_by_month_country_df.xlsx')

print('all the normalisation time :', (time.time() - start) / 60)
