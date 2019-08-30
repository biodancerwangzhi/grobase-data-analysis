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



def get_active_gene():
    pass

def get_fisher_test(arr1, arr2):
    pass  

def get_pi(promoter, gene_body):
    pass


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

def get_list_val_index_dic(list_a):
    dic_a = {}
    for i in range(len(list_a)):
        val = list_a[i]
        dic_a.setdefault(val, i)
    return dic_a 

def main():
    c = open('./specie_active_gene_num.txt', 'w') 
    for specie_dir in ['human', 'mouse' 'ara', 'cae', 'dro',  'mac', 'man',  'pan', 'pla', 'rat','sac', 'sch', 'zea']:
        ############### active gene
        print specie_dir
        c.write('########specie_dir:'+specie_dir+'\n')
        gene_body_file = '../%s/gene_body_feature_count_matrix.csv'%specie_dir
        intergenic_file = '../%s/intergenic_feature_count_matrix.csv'%specie_dir
        promoter_file = '../%s/promoter_feature_count_matrix.csv'%specie_dir

        gene_body_sample_name_list = open(gene_body_file).readline().strip().split(',')
        gene_body_sample_val_index_dict = get_list_val_index_dic(gene_body_sample_name_list)

        intergenic_sample_name_list = open(intergenic_file).readline().strip().split(',')
        intergenic_sample_val_index_dict = get_list_val_index_dic(intergenic_sample_name_list)

        promoter_sample_name_list = open(promoter_file).readline().strip().split(',')
        promoter_sample_val_index_dict = get_list_val_index_dic(promoter_sample_name_list)  

        intergenic_total_length = np.sum(read_csv_col_ndarray(intergenic_file, 1, 1))

        for sample_name in promoter_sample_name_list:
            if sample_name in ['ensembl_id', 'length', 'hugo_id', 'biotype']:
                continue
            gene_body_density = get_density_list(gene_body_file, gene_body_sample_val_index_dict[sample_name])
            intergenic_density = get_density_list(intergenic_file, intergenic_sample_val_index_dict[sample_name])
            # promoter_density = get_density_list(promoter_file, promoter_sample_val_index_dict[sample_name])

            #import matplotlib.pyplot as plt
        
            # intergenic_density = [x for x in intergenic_density if x>0 and x <1]
            intergenic_density = [x for x in intergenic_density if x>0]
            intergenic_density = sorted(intergenic_density)[:int(len(intergenic_density)*0.95)]

            back_val = np.mean(intergenic_density)
            # gene_body_density = [x for x in gene_body_density if x>2*back_val]

            # if 0.001 > back_val or back_val > 0.01:
                # print back_val
            

            p_value= get_bulk_poisson_probability(gene_body_density, back_val, 1000)   
            # promoter_density = [x for x in promoter_density if x <1 and x>0]
            c.write(sample_name + '\t' + str(len([x for x in p_value if x<0.001])) + '\n')  

    c.close()

if __name__ == '__main__':
    main()
