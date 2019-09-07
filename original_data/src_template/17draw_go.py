import os

import numpy as np



def read_csv_col_ndarray(csvDir, col, row_start):
    return np.loadtxt(csvDir,delimiter=',',unpack=True, skiprows=row_start, usecols=col)

def get_list_val_index_dic(list_a):
    dic_a = {}
    for i in range(len(list_a)):
        val = list_a[i]
        dic_a.setdefault(val, i)
    return dic_a 

def get_set_classification(set_total_gene, set_active_gene, set_paused_gene):
    set_active_paused = set_active_gene & set_paused_gene
    set_active_unpaused = set_active_gene & (set_total_gene - set_paused_gene)
    set_unactive_paused = (set_total_gene -set_active_gene) & set_paused_gene
    set_unactive_unpaused = (set_total_gene -set_active_gene) & (set_total_gene - set_paused_gene)
    return set_active_paused, set_active_unpaused, set_unactive_paused, set_unactive_unpaused

def write_gene_to_file(set_gene, file_name):
    with open(file_name, 'w') as f:
        for i in set_gene:
            f.write(i+'\n')

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
importr('org.At.tair.db')
importr('clusterProfiler')

def go_analysis(set_gene_file, output):
    DE_list = robjects.StrVector(list(set_gene_file))

    rcode_bp = 'enrichGO(gene= %s,keyType="TAIR",OrgDb =%s ,ont ="%s", \
    pAdjustMethod = "BH",pvalueCutoff = 0.01, qvalueCutoff = 0.05, \
    readable = TRUE)'%(DE_list.r_repr(), 'org.At.tair.db', 'BP')
    rcode_cc = 'enrichGO(gene= %s,keyType="TAIR",OrgDb =%s ,ont ="%s", \
    pAdjustMethod = "BH",pvalueCutoff = 0.01, qvalueCutoff = 0.05, \
    readable = TRUE)'%(DE_list.r_repr(), 'org.At.tair.db', 'CC') 
    rcode_mf = 'enrichGO(gene= %s,keyType="TAIR",OrgDb =%s ,ont ="%s", \
    pAdjustMethod = "BH",pvalueCutoff = 0.01, qvalueCutoff = 0.05, \
    readable = TRUE)'%(DE_list.r_repr(), 'org.At.tair.db', 'MF')   

    ego_bp = robjects.r(rcode_bp)
    ego_cc = robjects.r(rcode_cc)
    ego_mf = robjects.r(rcode_mf)

    c = open(output + '_bp.txt', 'w')
    c.write(str(ego_bp.r_repr()))
    c.close()
    c = open(output + '_cc.txt', 'w')
    c.write(str(ego_cc.r_repr()))
    c.close()
    c = open(output + '_mf.txt', 'w')
    c.write(str(ego_mf.r_repr()))
    c.close()

def main():
    output_dir = '../results/active_pause_go/'
    os.system('mkdir ' + output_dir)
    gene_file = '../results/active_gene_p_value.csv'

    gene_list = np.loadtxt(gene_file, delimiter=',',unpack=True, skiprows=1, usecols=0, dtype='string')
    set_total_gene = set(gene_list)

    active_sample_name_list = open(gene_file).readline().strip().split(',')
    active_sample_val_index_dict = get_list_val_index_dic(active_sample_name_list)

    pi_file = '../results/pi_gene_fc.csv'
    pi_sample_name_list = open(pi_file).readline().strip().split(',')
    pi_sample_val_index_dict = get_list_val_index_dic(pi_sample_name_list)   
    fisher_file = '../results/fisher_pvalue_gene.csv'
    fisher_sample_name_list = open(fisher_file).readline().strip().split(',')
    fisher_sample_val_index_dict = get_list_val_index_dic(fisher_sample_name_list)   

    for sample_name in active_sample_name_list:
        if sample_name in ['ensembl_id', 'length', 'hugo_id', 'biotype']:
            continue
        else:
            active_gene_score = read_csv_col_ndarray(gene_file, active_sample_val_index_dict[sample_name], 1)
            active_gene_set = set(gene_list[active_gene_score < 0.001])
            # pi pause
            paused_gene_pi = read_csv_col_ndarray(pi_file, pi_sample_val_index_dict[sample_name], 1)
            paused_gene_pi_set = set(gene_list[paused_gene_pi > 4])            
            set_active_pause, set_active_unpause, set_unactive_pause, set_unactive_unpause = get_set_classification(set_total_gene, active_gene_set, paused_gene_pi_set)
            go_analysis(set_active_pause, output_dir + 'pi_'+ sample_name+ '_active_pause')
            go_analysis(set_active_unpause, output_dir + 'pi_'+ sample_name+ '_active_unpause')
            go_analysis(set_unactive_pause, output_dir + 'pi_'+ sample_name+ '_unactive_pause')
            go_analysis(set_unactive_unpause, output_dir + 'pi_'+ sample_name+ '_unactive_unpause')
            # fisher pause
            paused_gene_fisher = read_csv_col_ndarray(fisher_file, fisher_sample_val_index_dict[sample_name], 1)
            paused_gene_fisher_set = set(gene_list[paused_gene_fisher < 0.001])            
            set_active_pause, set_active_unpause, set_unactive_pause, set_unactive_unpause = get_set_classification(set_total_gene, active_gene_set, paused_gene_fisher_set)
            go_analysis(set_active_pause, output_dir +'fisher_'+ sample_name+ '_active_pause')
            go_analysis(set_active_unpause, output_dir +'fisher_'+ sample_name+ '_active_unpause')
            go_analysis(set_unactive_pause, output_dir +'fisher_'+ sample_name+ '_unactive_pause')
            go_analysis(set_unactive_unpause, output_dir +'fisher_'+ sample_name+ '_unactive_unpause')

if __name__ == '__main__':
    main()