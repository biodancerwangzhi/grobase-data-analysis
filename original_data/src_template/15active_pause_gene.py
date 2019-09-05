# -*- coding:utf-8 -*-
import os, sys

import numpy as np
from scipy.misc import factorial
import scipy.stats as stats
import csv

def read_csv_as_ndarray(csvDir):
    csv_file = open(csvDir)
    csv_reader_list = csv.reader(csv_file)
    tmp_list = []
    for ln in csv_reader_list:
        tmp_list.append(ln)
    csv_ndarray = np.array(tmp_list)
    return csv_ndarray

def read_csv_col_ndarray(csvDir, col, row_start):
    return np.loadtxt(csvDir,delimiter=',',unpack=True, skiprows=row_start, usecols=col)

def read_csv_row_ndarray(csvDir, row, col_start):
    rowData = open(csvDir)
    newLine = rowData.readline()
    for i in range(row):
        newLine =rowData.readline()
    return np.array(newLine.strip('\n').split(',')[col_start:])

def get_list_val_index_dic(list_a):
    dic_a = {}
    for i in range(len(list_a)):
        val = list_a[i]
        dic_a.setdefault(val, i)
    return dic_a 
    
def get_poisson_probability(k, _lambda):
    l = _lambda   # float
    k = int(k)      # int
    if k <= 170:    # factorial(k) over inf
        p_value = l ** k * np.exp(-l) / factorial(k)
    else:
        #print k, _lambda
        p_value = 0
    return p_value

def get_bulk_poisson_probability(k_list, _lambda, folder_val):
    l = _lambda*folder_val
    k_list = k_list*folder_val
    k_list = k_list.astype(int)
    k_list = k_list.tolist()
    l_list = [l]*len(k_list)
    return map(get_poisson_probability, k_list, l_list)

def get_density_list(csvDir, sample_col):
    length_array = read_csv_col_ndarray(csvDir, 1, 1)
    count_array = read_csv_col_ndarray(csvDir, sample_col, 1)
    return count_array/length_array


def get_active_gene(gene_body_file, intergenic_file, output_file):
    gene_body_sample_name_list = open(gene_body_file).readline().strip().split(',')
    gene_body_sample_val_index_dict = get_list_val_index_dic(gene_body_sample_name_list)
    intergenic_sample_name_list = open(intergenic_file).readline().strip().split(',')
    intergenic_sample_val_index_dict = get_list_val_index_dic(intergenic_sample_name_list)
    #intergenic_total_length = np.sum(read_csv_col_ndarray(intergenic_file, 1, 1))

    # gene name
    old_column = np.loadtxt(gene_body_file, delimiter=',',unpack=True, skiprows=1, usecols=0, dtype='string')
    
    sample_list = ['ensembl_id']
    for sample_name in gene_body_sample_name_list:
        if sample_name in ['ensembl_id', 'length', 'hugo_id', 'biotype']:
            continue
        sample_list.append(sample_name)
        gene_body_density = get_density_list(gene_body_file, gene_body_sample_val_index_dict[sample_name])
        intergenic_density = get_density_list(intergenic_file, intergenic_sample_val_index_dict[sample_name])
        intergenic_density = [x for x in intergenic_density if x>0]
        # 去掉异常值
        intergenic_density = sorted(intergenic_density)[:int(len(intergenic_density)*0.95)]
        back_val = np.mean(intergenic_density)
        p_value= get_bulk_poisson_probability(gene_body_density, back_val, 1000)   
        old_column = np.vstack([old_column, p_value])
    
    old_column_T = old_column.T
    old_column_T = np.vstack([np.array(sample_list), old_column_T])
    np.savetxt(output_file, old_column_T, delimiter=',', fmt='%s')

def get_pi(promoter_file, gene_body_file, output_file):
    promoter_sample_name_list = open(promoter_file).readline().strip().split(',')
    promoter_sample_val_index_dict = get_list_val_index_dic(promoter_sample_name_list) 
    gene_body_sample_name_list = open(gene_body_file).readline().strip().split(',')
    gene_body_sample_val_index_dict = get_list_val_index_dic(gene_body_sample_name_list)

    # gene name
    old_column = np.loadtxt(gene_body_file, delimiter=',',unpack=True, skiprows=1, usecols=0, dtype='string')

    sample_list = ['ensembl_id']
    for sample_name in gene_body_sample_name_list:
        if sample_name in ['ensembl_id', 'length', 'hugo_id', 'biotype']:
            continue
        sample_list.append(sample_name)    
        gene_body_density = get_density_list(gene_body_file, gene_body_sample_val_index_dict[sample_name])
        promoter_density = get_density_list(promoter_file, promoter_sample_val_index_dict[sample_name])
        pi = []
        for x, y in zip(promoter_density, gene_body_density):
            if y == 0:
                pi.append(0)
            else:
                pi.append(float(x)/y)
        pi = np.array(pi)
        old_column = np.vstack([old_column, pi])
    old_column_T = old_column.T
    old_column_T = np.vstack([np.array(sample_list), old_column_T])
    np.savetxt(output_file, old_column_T, delimiter=',', fmt='%s')     

def fisher_test(gene_body_count, promoter_proximal_peak_count, float_lambda, gene_body_length, promoter_window_length):
    list_a = [int(gene_body_count), int(promoter_proximal_peak_count)]
    list_b = [float_lambda* gene_body_length, float_lambda*promoter_window_length]
    oddsratio, pvalue = stats.fisher_exact([list_a, list_b])  
    return pvalue 

def get_fisher_test(promoter_windows_file, gene_body_file,intergenic_file, output_file):
    gene_body_sample_name_list = open(gene_body_file).readline().strip().split(',')
    gene_body_sample_val_index_dict = get_list_val_index_dic(gene_body_sample_name_list)
    intergenic_sample_name_list = open(intergenic_file).readline().strip().split(',')
    intergenic_sample_val_index_dict = get_list_val_index_dic(intergenic_sample_name_list)
    # intergenic_total_length = np.sum(read_csv_col_ndarray(intergenic_file, 1, 1))        
    promoter_sample_name_list = open(promoter_windows_file).readline().strip().split(',')
    promoter_sample_val_index_dict = get_list_val_index_dic(promoter_sample_name_list)    

    old_column = np.loadtxt(gene_body_file, delimiter=',',unpack=True, skiprows=1, usecols=0, dtype='string')
    
    body_length = np.loadtxt(gene_body_file, delimiter=',',unpack=True, skiprows=1, usecols=1)

    promoter_length = []
    promoter_windows_length = np.loadtxt(promoter_windows_file, delimiter=',',unpack=True, skiprows=1, usecols=1) 
    int_total = len(promoter_windows_length) 
    int_num =2
    while int_num< int_total-1:
        promoter_length.append(promoter_windows_length[int_num])
        int_num +=100

    sample_list = ['ensembl_id']
    for sample_name in gene_body_sample_name_list:
        if sample_name in ['ensembl_id', 'length', 'hugo_id', 'biotype']:
            continue
        sample_list.append(sample_name)
        gene_body_density = get_density_list(gene_body_file, gene_body_sample_val_index_dict[sample_name])
        intergenic_density = get_density_list(intergenic_file, intergenic_sample_val_index_dict[sample_name])
        intergenic_density = [x for x in intergenic_density if x>0]
        intergenic_density = sorted(intergenic_density)[:int(len(intergenic_density)*0.95)]
        back_val = np.mean(intergenic_density)

        promoter_windows_density = get_density_list(promoter_windows_file, promoter_sample_val_index_dict[sample_name])

        int_total = len(promoter_windows_density)
        max_promoter_window_density = []
        # feature count row from 2
        int_num =2
        while int_num< int_total-1:
            int_max_count = np.max(promoter_windows_density[int_num:int_num+100])
            max_promoter_window_density.append(int_max_count)
            int_num +=100
        
        gene_body_count = gene_body_density * body_length
        promoter_count = np.array(max_promoter_window_density) * promoter_length

        l_list = [back_val]*len(gene_body_count)
        fisher_pvalue = map(fisher_test, gene_body_count, promoter_count, l_list, body_length, promoter_length)
        old_column = np.vstack([old_column, fisher_pvalue])
    
    old_column_T = old_column.T
    old_column_T = np.vstack([np.array(sample_list), old_column_T])
    np.savetxt(output_file, old_column_T, delimiter=',', fmt='%s')        



def main():
    gene_body_file = '../results/gene_body_feature_count_matrix.csv'
    intergenic_file = '../results/intergenic_feature_count_matrix.csv'
    active_gene_output_file = '../results/active_gene_p_value.csv'
    get_active_gene(gene_body_file, intergenic_file, active_gene_output_file)

    promoter_file = '../results/promoter_feature_count_matrix.csv'
    pi_ouput_file = '../results/pi_gene_fc.csv'
    get_pi(promoter_file, gene_body_file, pi_ouput_file)

    promoter_windows_file = '../results/promoter_windows_feature_count_matrix.csv'
    fisher_output_file = '../results/fisher_pvalue_gene.csv'
    get_fisher_test(promoter_windows_file, gene_body_file,intergenic_file, fisher_output_file)

if __name__ == '__main__':
    main()
