"""
validation script.

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
path_to_save = path_with_data
if os.path.exists(path_to_tools) and os.path.exists(path_to_save) and os.path.exists(path_with_data):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')
    exit()

# ------------------------------------------------------------------------ </editor-fold>


# ------------------------------------------------------------------------ </editor-fold>


###############   the HLA acording to most frequent  from paper ###################


hla_list = ['HLA-A01:01','HLA-A02:01','HLA-A02:05','HLA-A03:01','HLA-A11:01','HLA-A23:01','HLA-A24:02','HLA-A25:01']



# ------------------------------------------------------------------------ </editor-fold>

# get all 9 mers from fasta refrence


pep_list = ['NALRITFGG','ALRITFGGP','LRITFGGPS','RITFGGPSD','ITFGGPSDS','TFGGPSDST','FGGPSDSTG','GGPSDSTGS','GPSDSTGSN','PSDSTGSNQ','SDSTGSNQN','DSTGSNQNG','STGSNQNGG','TGSNQNGGA','GSNQNGGAR','SNQNGGARS','NQNGGARSK','QNGGARSKQ']


# ------------------------------------------------------------------------ </editor-fold>













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




def mix_mhc_pred(pep_list):
    global hla_list
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
                               var_name='HLA', value_name='Rank')

    return mix_pred_df_melt


def mhc_flurry(pep_list):
    global hla_list


    dict_input = {}
    # The
    # input
    # CSV
    # file is expected
    # to
    # contain
    # columns “allele”, “peptide”, and, optionally, “n_flank”, and “c_flank”.

    output_file = path_with_data + 'flurry_pred_out.csv'
    command = 'mhcflurry-predict'+ INPUT.csv + ' –out '+output_file
    subprocess.run(command, shell=True)

    flurry_pred_df = pd.read_csv(output_file, sep='\s+', comment='#')




def net_mhc_pan(pep_list):

    global hla_list

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


    predict_df.rename(columns=col_dict,inplace=True)
    predict_df.reset_index(inplace=True)
    mix_pred_df_melt = pd.melt(predict_df, id_vars=['Peptide'], value_vars=predict_df.columns[:],
                               var_name='HLA', value_name='Rank')

    return predict_df
# ------------------------------------------------------------------------ </editor-fold>
