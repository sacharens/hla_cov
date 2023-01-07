"""
ensemble hla script

this scrip will take in a list of peptides and the alleles that are predicted to bind to them.
then the script will invoke the following predictores to make sure they are sound.
* mix mhc pred
* mhc flurry

the output will be the how many of the alleles are predicted to bind from the given set

by: Sinai Sacharen

"""

# -------------------------""" imports """------------------------------ <editor-fold>
import sys
import subprocess
import os
import time
import pandas as pd
from Bio import SeqIO
from scipy.stats import hmean
import concurrent.futures

from Bio.Seq import Seq
import re
import io
import math
import ast
from sklearn import preprocessing
# from mhcflurry import Class1AffinityPredictor
# predictor = Class1AffinityPredictor.load()


import numpy as np
import ast

""" paths """
# ------------------------------------------------------------------------ <editor-fold>
# C:\Users\User\Documents\COVID-19\peptide_valatadion



path_to_tools = '/home/sacharen/'
path_with_data = '/home/sacharen/Documents/ansambel_hla/'
path_to_save = path_with_data+'protensemble/'
if os.path.exists(path_to_tools) and os.path.exists(path_to_save) and os.path.exists(path_with_data):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')
    exit()

# ------------------------------------------------------------------------ </editor-fold>


# -----------------------------import files--------------------------------- <editor-fold>

prevalent_allel_file = path_with_data + 'prevalent_alleles'
# full_alele_file = pd.read_csv(path_with_data + 'all_netMHCpan_alleles.txt', sep='\s+', header=None)
path_with_all_protien_fa = '/home/sacharen/Documents/find_ref/sars_protien_fasta/'
ref_df = pd.read_excel('/home/sacharen/Documents/find_ref/new_ref_df.xlsx')
###############   the HLA acording to most frequent  from paper ###################


# hla_list = ['HLA-A01:01','HLA-A02:01','HLA-A02:05','HLA-A03:01','HLA-A11:01','HLA-A23:01','HLA-A24:02','HLA-A25:01']



# ------------------------------------------------------------------------ <editor-fold>

# get all 9 mers from fasta refrence

#
# pep_list = ['NALRITFGG','ALRITFGGP','LRITFGGPS','RITFGGPSD','ITFGGPSDS','TFGGPSDST','FGGPSDSTG','GGPSDSTGS','GPSDSTGSN','PSDSTGSNQ','SDSTGSNQN','DSTGSNQNG','STGSNQNGG','TGSNQNGGA','GSNQNGGAR','SNQNGGARS','NQNGGARSK','QNGGARSKQ']
#

# ----------------------------------FUNCTIONS-------------------------------------- <editor-fold>


def read_in_allels(prealent_allel_file_and_path, list_format = False, str_format = False):
    """
    tis function will read in the file and return all the allelels in it in the following format:
    'HLA-B40:01,HLA-B58:01,HLA-B15:01'
    so the netHmcpan predicter can read it in
    :param prealent_allel_file_and_path:
    :return: string
    """
    with open(prealent_allel_file_and_path) as file_in:
        lines = ''
        list_lines = []
        for line in file_in:
            temp_line = line.rstrip() + ','
            if len(temp_line) < 2:
                continue
            else:
                lines += temp_line
                list_lines.append(temp_line[:-1])

        if list_format:
            return list_lines
        if str_format:
            return lines[:-1]


def all_available_kmers_f_seq(seq, k):
    """
    Gets a string. Returns a list with all kmers (motifs) contained in the string,
    @param seq: string
    @param k: int. the size of motifs to extract (kmers)
    @return: string

    """
    kmer_list = []
    for start in range(len(seq) - k + 1):
        kmer = seq[start:start+k]
        kmer_list.append(kmer)

    return kmer_list





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


def majority_voting(ensemble_df):
    """
    this function will take in the ensemble df and will apply majority voting to it
    :param ensemble_df: the ensemble df
    :return: 2 columns  1- if its a majority vote 2- the geometric mean of the scores of the majority vote
    """
    # see if 2 or more of the scores are smaller than 2
    # if so then it is a majority vote
    # if not then it is not a majority vote
    # if it is a majority vote then take the geometric mean of the scores
    # if not then take the mean of the scores
    count = 0
    score_rank_list = []
    for score in ['mixMHC Rank', 'mhcflurry_presentation_percentile', 'netMHCpan Rank']:
        pred_score = ensemble_df[score]
        if pred_score <= 2:
            count += 1
            score_rank_list.append(score)
    if count > 2:
        ensemble_df['majority vote'] = 1
        ensemble_df['majority vote score'] = hmean(ensemble_df[score_rank_list])
    else:
        ensemble_df['majority vote'] = 0
        ensemble_df['majority vote score'] = hmean(ensemble_df[
                                                       ['mixMHC Rank', 'mhcflurry_presentation_percentile',
                                                        'netMHCpan Rank']])
    return ensemble_df








def creat_pep_file_from_list(pep_list_f, file_path_and_name_f):
    """
    the net MHC predictor needs to get a peptide file as an input so this function will creat that file
    from a given list of peptides
    :param pep_list_f:
    :param file_path_and_name_f:
    """
    epfile = open(file_path_and_name_f, 'w')

    # epetope = 'QYIKWPWYI'
    for p_f in pep_list_f:
        s = p_f + "\n"

        # Writing a string to file
        epfile.write(s)

    # Closing file
    epfile.close()


def make_mix_hla_fomat(hla_list_f):
    """
    this function will take in a list of alleles and will returen one
    that has the right format for th mixpred predictor
    A2501,B0801,B1801  (str)
    """

    hla_str = ''
    for hla in hla_list_f:
        hla_str += hla[4:].replace(':', '') + ','

    return hla_str[:-1]




def mix_mhc_pred(pep_list,hla_list):

    HLA_str = make_mix_hla_fomat(hla_list)
    temp_pep_file = path_with_data + "mix_pep_file.txt"
    creat_pep_file_from_list(pep_list, temp_pep_file)
    output_file = path_with_data + 'mix_pred_out.csv'
    command = path_to_tools+'MixMHCpred/' + 'MixMHCpred -i ' + temp_pep_file + ' -o ' + output_file + ' -a ' + HLA_str
    subprocess.run(command, shell=True)

    mix_pred_df = pd.read_csv(output_file, sep='\s+', comment='#')
    mix_pred_df.set_index('Peptide', inplace=True)
    mix_pred_df.drop('%Rank_bestAllele', inplace=True, axis=1)
    mix_pred_df = mix_pred_df.filter(regex='^%Ra', axis=1)
    pd.melt(mix_pred_df,id_vars=[])

    # turn the hla back to regular format:
    col_dict = {}
    for col in mix_pred_df.columns:
        col_dict[col] = 'HLA-'+col[6:7]+col[7:9]+':'+col[9:11]

    mix_pred_df.rename(columns=col_dict,inplace=True)
    mix_pred_df.reset_index(inplace=True)
    mix_pred_df_melt = pd.melt(mix_pred_df, id_vars=['Peptide'], value_vars=mix_pred_df.columns[1:],
                               var_name='allele', value_name='Rank')

    mix_pred_df_melt.rename(columns={'allele':'HLA','Rank':'mixMHC Rank'},inplace = True)

    return mix_pred_df_melt


def mhc_flurry(pep_list,hla_list):
    dict_input = {}
    for h in hla_list:
        dict_input[h] = pep_list

    input_csv = pd.DataFrame.from_dict(dict_input)
    input_csv = pd.melt(input_csv,value_vars=hla_list,var_name='allele', value_name='peptide')
    input_csv.set_index('allele',drop=True,inplace=True)
    input_csv.to_csv(path_with_data + "flurry_input_file.csv")
    temp_pep_file = path_with_data + "flurry_input_file.csv"
    output_file = path_with_data + 'flurry_pred_out.csv'
    # touch_comm = 'touch '+output_file
    # subprocess.run(touch_comm, shell=True)
    command = 'mhcflurry-predict '+ temp_pep_file + ' > '+ output_file
    subprocess.run(command, shell=True)

    flurry_pred_df = pd.read_csv(output_file, sep=',',skiprows=4)
    return_df = flurry_pred_df[['allele','peptide','mhcflurry_presentation_percentile']].copy()
    return_df.rename(columns={'allele':'HLA','peptide':'Peptide'},inplace = True)
    return return_df


def net_mhc_pan(pep_list,hla_list):
    hla_str = ','.join(hla_list)
    temp_pep_file = path_with_data + "net_pep_file.txt"
    creat_pep_file_from_list(pep_list, temp_pep_file)
    output_file = path_with_data + 'net_mhc_out.csv'
    command = path_to_tools+'netMHCpan-4.1/' + "netMHCpan " + " -a " + hla_str + " -p " + temp_pep_file + " -xls " + " -xlsfile " + output_file

    subprocess.check_output('%s' % command, shell=True)

    predict_df = pd.read_csv(output_file, sep='\t', header=[0, 1])
    new_col_dict = {}
    current_hla = ''
    for col in predict_df.columns:
        print(col)
        if col[0].startswith('HLA'):
            current_hla = col[0]
        if current_hla == '':
            new_col_dict[col] = col[1]
        if col[1] == 'EL_Rank':
            new_col_dict[col] = current_hla

    predict_df.columns = predict_df.columns.to_flat_index()
    predict_df.rename(columns=new_col_dict, inplace=True)
    for col in predict_df.columns:
        if type(col) is tuple:
            predict_df.drop(col, axis=1, inplace=True)

    predict_df.reset_index(inplace=True)
    predict_df = pd.melt(predict_df, id_vars=['Peptide'], value_vars=hla_list,
                               var_name='HLA', value_name='Rank')
    predict_df.rename(columns={'Rank':'netMHCpan Rank'},inplace = True)

    return predict_df
# ------------------------------------------------------------------------ </editor-fold>

HLA_list = read_in_allels(prevalent_allel_file,list_format=True)

for p in ['NSP10', 'NSP9', 'NSP8', 'NSP7', 'NSP6', 'NSP5', 'NSP4', 'NSP3', 'NSP2', 'NSP1', 'NS7b', 'NSP11', 'NSP16',
          'NSP15', 'NSP14'
    , 'NSP13', 'NSP12', 'N', 'NS8', 'NS7a', 'NS6', 'M', 'E', 'NS3', 'Spike']:

    print('processing protein '+p+'...........................')
    file_and_path_with_fasta = path_with_all_protien_fa + p + "_fasta.fa"
    with open(file_and_path_with_fasta, "r", encoding='latin1') as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            seq_str = str(record.seq)
            PEP_list = all_available_kmers_f_seq(seq=seq_str, k=9)

    mix_df = mix_mhc_pred(pep_list=PEP_list, hla_list=HLA_list)
    flurry_df = mhc_flurry(pep_list=PEP_list, hla_list=HLA_list)
    net_df = net_mhc_pan(pep_list=PEP_list, hla_list=HLA_list)
    ensemble_df = pd.concat([mix_df, flurry_df,net_df], axis=1)


    # drop duplicates columns from ensemble_df

    ensemble_df = ensemble_df.loc[:, ~ensemble_df.columns.duplicated()]

    ensemble_df_parallel = parallelize_dataframe(df=ensemble_df, func=majority_voting, n_cores=8)

    ensemble_df_parallel.to_excel(path_to_save + p + '_ensemble_df.xlsx')
