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

def get_gene_window_dict(file, windows_num = 390):
    dict_a = {}
    list_ln = open(file).readlines()	
    int_total = len(list_ln)
    int_num =2
    while int_num< int_total-1:
        list_gene_windows = list_ln[int_num:int_num+windows_num]
        gene_name = '_'.join(list_gene_windows[0].split('_')[:-1])
        # 
        list_windows_density = [float(i.strip('\n').split('\t')[-1])/(int(i.strip('\n').split('\t')[-2])) for i in list_gene_windows]
		# int_max_count = max(list_bin_count)
        dict_a.setdefault(gene_name, list_windows_density)
        int_num +=windows_num
    return dict_a

def get_active_gene():
    pass

def get_fisher_test(arr1, arr2):
    pass  

def get_pi(promoter, gene_body):
    pass

'''
def get_gsm_lambda_dict(dir_file):
    dict_a ={}
    list_file = os.listdir(dir_file)
    for file_name in list_file:
        if file_name.endswith('summary'):
            continue
        total_length = 0
        total_count = 0
        for ln in open(dir_file+file_name):
            if ln[0] == '#' or ln.split('\t')[0] == 'Geneid':
                continue  
            total_length += int(ln.split('\t')[-2])
            total_count += int(ln.split('\t')[-1].strip('\n'))
        float_background_density = float(total_count)/total_length
        dict_a.setdefault(file_name.split('_')[0], float_background_density)
    return dict_a    
# dict_gsm_lambda = get_gsm_lambda_dict('../results/intergenic_region_feature_count/')

def get_gsm_lambda(str_gsm_id):
    return dict_gsm_lambda[str_gsm_id]

mer_30_grch38 = 2484365113  
mer_50_grch38 = 2701171919
def get_gsm_lambda(str_gsm_id):
    return 0.03*dict_total_mapped_reads[str_gsm_id]/mer_30_grch38/2
    #return 0.04


def get_gsm_lambda(str_gsm_id):  
    print dict_total_mapped_reads[str_gsm_id]
    print float(dict_total_mapped_reads[str_gsm_id])/3200000000/2
    return float(dict_total_mapped_reads[str_gsm_id])/3200000000/2.0
'''    

def get_poisson_probability(k, _lambda):
    l = _lambda
    try:
        p_value = l ** k * np.exp(-l) / factorial(k)
        #print p_value
    except:
        # print  np.exp(-l) / factorial(k)
        p_value = float("inf")
    return p_value

# transcribed or not
def get_gsm_gene_poisson_probability_dict(str_file_name,float_lambda):
    dict_a={}
    num = 0
    for ln in open(str_file_name):
        if ln[0] == '#' or ln.split('\t')[0] == 'Geneid':
            continue   
        list_ln = ln.split('\t')
        feature_type = list_ln[0].split('_')[-1]
        gene_name = '_'.join(list_ln[0].split('\t')[0].split('_')[0:-1])
        if feature_type == 'tes':
            int_count = int(list_ln[-1].strip('\n'))
            int_length = int(list_ln[-2])
            gene_body_density = float(int_count)/int_length
            p_value = get_poisson_probability(gene_body_density*1000, float_lambda*1000)   
            if p_value< 0.01:
                num+=1
            dict_a.setdefault(gene_name, p_value)
    return dict_a
#dict_gene_p_value = get_gsm_gene_poisson_probability_dict('../results/onek2_tes_feature_count/GSM2296633_count.txt') 
#print dict_gene_p_value  

# get max bin in tss+- 1k。{gene name： max value}
def get_gene_promoter_proximal_peak_dict(str_file_name):
	dict_a = {}
	# file rows number
	int_total = 8597942  # change for different species
	list_ln = open(str_file_name).readlines()	
	# feature count row from 2
	int_num =2
	while int_num< int_total-1:
		list_ln_gene = list_ln[int_num:int_num+390]
        # str_gene_name = list_ln_gene[0].split('\t')[0].split('_')[0]  bug because gene name have more than one '_'
		str_gene_name = '_'.join(list_ln_gene[0].split('\t')[0].split('_')[0:-1])
		list_bin_count = [int(i.strip('\n').split('\t')[-1])for i in list_ln_gene]
		int_max_count = max(list_bin_count)
		dict_a.setdefault(str_gene_name, int_max_count)
		int_num +=390
	return dict_a
#print get_gene_promoter_proximal_peak_dict('../results/promoter_proximal_peak_feature_count/GSM2545173_gene_count.txt')

def get_gene_body_count_length_dict(str_file_name):
    dict_a={}
    for ln in open(str_file_name):
        if ln[0] == '#' or ln.split('\t')[0] == 'Geneid':
            continue   
        list_ln = ln.split('\t')
        feature_type = list_ln[0].split('_')[-1]
        gene_name = '_'.join(list_ln[0].split('_')[0:-1])
        if feature_type == 'tes':
            int_count = int(list_ln[-1].strip('\n'))
            int_length = int(list_ln[-2])
            dict_a.setdefault(gene_name, [int_count, int_length])
    return dict_a  
#print get_gene_body_density_dict('../results/onek2_tes_feature_count/GSM2296633_count.txt')        

def get_fisher_exact_test_dict(gsm_id, float_lambda):
    dict_a={}
    dict_gene_promoter_proximal_peak = get_gene_promoter_proximal_peak_dict('../results/promoter_proximal_peak_feature_count/%s_gene_count.txt'%gsm_id)
    #print dict_gene_promoter_proximal_peak.keys()[:10]
    dict_gene_body_list_count_length = get_gene_body_count_length_dict('../results/onek2_tes_feature_count/%s_count.txt'%gsm_id) 
    #print dict_gene_body_list_count_length.keys()[:10]    
    #float_lambda
    for gene_name in dict_gene_body_list_count_length.keys():
        # promoter proximal may be drop, because tss-1k < 0
        try:
            list_count_length = dict_gene_body_list_count_length[gene_name]
            gene_body_count = list_count_length[0]
            gene_body_length = list_count_length[1]
            promoter_proximal_peak_count = dict_gene_promoter_proximal_peak[gene_name]
            list_a = [gene_body_count, promoter_proximal_peak_count]
            list_b = [float_lambda* gene_body_length, float_lambda*50]
            oddsratio, pvalue = stats.fisher_exact([list_a, list_b])
            dict_a.setdefault(gene_name, pvalue)
        except:
            pass        
    return dict_a   

def get_pi_dict(gsm_id):    
    dict_a={}
    dict_gene_promoter_proximal_peak = get_gene_promoter_proximal_peak_dict('../results/promoter_proximal_peak_feature_count/%s_gene_count.txt'%gsm_id)
    dict_gene_body_list_count_length = get_gene_body_count_length_dict('../results/onek2_tes_feature_count/%s_count.txt'%gsm_id)  
    #float_lambda
    for gene_name in dict_gene_body_list_count_length.keys():
        try:
            list_count_length = dict_gene_body_list_count_length[gene_name]
            gene_body_count = list_count_length[0]
            gene_body_length = list_count_length[1]
            promoter_proximal_peak_count = dict_gene_promoter_proximal_peak[gene_name]
            if gene_body_count == 0:
                float_pi = float("-inf")
            else:
                float_pi = float(promoter_proximal_peak_count)/50/gene_body_count *gene_body_length
                # print 'float_pi', float_pi
            
            dict_a.setdefault(gene_name, [float_pi, promoter_proximal_peak_count, gene_body_count])
        except:
            pass        
    return dict_a  
 
import json 
'''
def main():   
    json_dir = '../results/active_pause_json/'
    output_dir = '../results/active_pause_list/'
    os.system('mkdir %s'%output_dir)
    for gsm_id in dict_total_mapped_reads.keys():
        print '###############'
        print gsm_id
        float_lambda = get_gsm_lambda(gsm_id)
        print 'background rate', float_lambda
        
        with open(json_dir+'fisher_test_'+gsm_id+'.json', 'r') as f:
            dict_fisher_test = json.load(fp=f)   
        list_fisher_test = [y for x, y in zip(dict_fisher_test.values(), dict_fisher_test.keys()) if x <= 0.01]

        with open(json_dir+'pi_'+gsm_id+'.json', 'r') as g:
            dict_pi = json.load(fp=g)      
        list_pi = [y for x,y in zip(dict_pi.values(), dict_pi.keys()) if x >=2]
        
        c = open(output_dir+'paused_'+gsm_id+'.txt', 'w')
        for pause_gene_name in set(list_fisher_test) & set(list_pi):
            c.write(pause_gene_name+'\n')
        c.close()

        with open(json_dir+'active_'+gsm_id+'.json', 'r') as h:
            dict_active = json.load(fp=h)     
        
        list_active = [y for x,y in zip(dict_active.values(), dict_active.keys()) if x <=0.01]
        c = open(output_dir+'active_'+gsm_id+'.txt', 'w')        
        for active_gene_name in list_active:
            c.write(active_gene_name+'\n')
        c.close()
main()
    

def main():
    output_dir = '../results/active_pause_json/'
    os.system('mkdir %s'%output_dir)
    for gsm_id in dict_total_mapped_reads.keys():
        print '###############'
        print gsm_id
		# use intergenic region density as lambda
        float_lambda = get_gsm_lambda(gsm_id)
        #print 'background rate', float_lambda    
    
        #print len(dict_pi.keys()), 'PI over 2'
        #print len(dict_active.keys()), 'actived gene'
        
		# compute all gene fisher exact value to judge pause or not
        dict_fisher_test = get_fisher_exact_test_dict(gsm_id, float_lambda)
        with open(output_dir+'fisher_test_'+gsm_id+'.json', 'w') as f:
            json.dump(dict_fisher_test, f) 
        print len(dict_fisher_test.keys()), 'paused gene'
        dict_pi = get_pi_dict(gsm_id)
        with open(output_dir+'pi_'+gsm_id+'.json', 'w') as g:
            json.dump(dict_pi, g)         
        print len(dict_pi.keys()), 'PI over 2'
		
		# compute active poisson value
        dict_active = get_gsm_gene_poisson_probability_dict('../results/onek2_tes_feature_count/%s_count.txt'%gsm_id, float_lambda)
        with open(output_dir+'active_'+gsm_id+'.json', 'w') as h:
            json.dump(dict_active, h)            
        print len(dict_active.keys()), 'actived gene' 

main()
'''


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
    # active gene
    
    gene_body_file = '../human/gene_body_feature_count_matrix.csv'
    intergenic_file = '../human/intergenic_feature_count_matrix.csv'
    promoter_file = '../human/promoter_feature_count_matrix.csv'
    gene_body_sample_name_list = open(gene_body_file).readline().strip().split(',')
    gene_body_sample_val_index_dict = get_list_val_index_dic(gene_body_sample_name_list)

    intergenic_sample_name_list = open(intergenic_file).readline().strip().split(',')
    intergenic_sample_val_index_dict = get_list_val_index_dic(intergenic_sample_name_list)

    promoter_sample_name_list = open(promoter_file).readline().strip().split(',')
    promoter_sample_val_index_dict = get_list_val_index_dic(promoter_sample_name_list)    

    intergenic_total_length = np.sum(read_csv_col_ndarray(intergenic_file, 1, 1))

    # ???????????prmoter? gene body? intergenic??density??????  
    sample_name = 'GSM2552215'

    gene_body_density = get_density_list(gene_body_file, gene_body_sample_val_index_dict[sample_name])
    
    intergenic_density = get_density_list(intergenic_file, intergenic_sample_val_index_dict[sample_name])

    promoter_density = get_density_list(promoter_file, promoter_sample_val_index_dict[sample_name])

    import matplotlib.pyplot as plt
    data = [promoter_density, gene_body_density, intergenic_density]
    fig, ax = plt.subplots()
    ax.set_title('Multiple Samples with Different sizes')
    ax.boxplot(data)

    plt.show()

    '''
    sorted_intergenic_ndarray = sorted(intergenic_ndarray)
    sorted_intergenic_ndarray = sorted_intergenic_ndarray[int(len(intergenic_ndarray)*0.1): int(len(intergenic_ndarray)*0.9)]
    intergenic_total_count = np.sum(sorted_intergenic_ndarray)
    
    
    intergenic_density = intergenic_total_count / intergenic_total_length
    print intergenic_density
    print max(gene_body_density)
    print get_poisson_probability(gene_body_density[3], intergenic_density)
    # print intergenic_density
    '''

    '''
    promoter_windows_file = ''
    promoter_windows_name_list =read_csv_col_ndarray(promoter_windows_file, 0, 1)
    sample_list= open(promoter_windows_file).readline().strip('\n').split(',')[1:]

    # pause gene PI



    # pause gene fisher test
    '''
main()
