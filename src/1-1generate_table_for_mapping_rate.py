# -*- coding: utf-8 -*-
import json
import os

def extract_unique_multiple_align(bam_stat_dir):
    split_label = '############'   
    unique_multiple_total_dict={}
    split_log_list=open(bam_stat_dir).read().split(split_label)   
    for split_log_item in split_log_list[1:]:
        gsm_id = split_log_item.split('.bam')[0].strip()
        split_log_line_list=split_log_item.split('\n')
        for split_log_line in split_log_line_list:
            if split_log_line.startswith('Total records:'):
                total_reads=int(split_log_line.split('Total records:                          ')[1].split('\n')[0])
            if split_log_line.startswith('mapq < mapq_cut (non-unique):'):
                multi_reads=int(split_log_line.split('mapq < mapq_cut (non-unique):           ')[1].split('\n')[0])
            if split_log_line.startswith('mapq >= mapq_cut (unique):'):
                unique_reads=int(split_log_line.split('mapq >= mapq_cut (unique):              ')[1].split('\n')[0])
        unique_multiple_total_dict[gsm_id]= [(unique_reads+multi_reads)/float(total_reads), float(unique_reads)/(unique_reads + multi_reads)]
    return unique_multiple_total_dict    

def extract_feature_count_align(feature_count_log_dir):
    list_feature_map_rate = []
    list_gsm_id = []
    for ln in open(feature_count_log_dir):
        if ln.startswith('||    Successfully assigned alignments : '):
            list_feature_map_rate.append(ln.split(':')[1].split(')')[0])
        if ln.startswith('||             Output file : GSM'):
            list_gsm_id.append(ln.split(':')[1].split('.')[0])
    dict_gsm_feature_mapRate = {}        
    for gsm_id, map_rate in zip(list_gsm_id, list_feature_map_rate):
        if 'gene' in gsm_id:
            gsm_id = gsm_id.split('_')[0].strip()
            map_rate = map_rate.split('(')[1]
            dict_gsm_feature_mapRate.setdefault(gsm_id, map_rate)
    return dict_gsm_feature_mapRate

#dir_specie = '../primary_results/'
c = open('./1map_rate.csv', 'w')
first_ln = 'gsm_id!unique_map_rate!map_rate!featureCount_map_rate\n'
c.write(first_ln)

#list_specie_file = os.listdir(dir_specie)
list_specie_file = ['ara', 'cae', 'dro', 'human', 'mac', 'man','mouse','pan','pla', 'rat', 'sac','sch', 'zea']
dir_specie = './'
for specie_file in list_specie_file:
    dir_specie_bam_stat_log_dir = dir_specie + specie_file + '/logs/10bam_stat.log'
    dict_bam_stat_map_rate = extract_unique_multiple_align(dir_specie_bam_stat_log_dir)
    #dir_feature_count_log = dir_specie + specie_file + '/logs/11exons_gene_featureCount.log'
    #dict_feature_count_map_rate = extract_feature_count_align(dir_feature_count_log)
    
    for gsm_id in dict_bam_stat_map_rate.keys():
        map_unique_map_list = dict_bam_stat_map_rate[gsm_id]
        #feature_map_rate = dict_feature_count_map_rate[gsm_id]
        feature_count_log_dir='/data1/szx/%s/results/exons_gene_featureCount/%s_gene_count.txt.summary'%(specie_file, gsm_id)
        try:
            for ln in open(feature_count_log_dir):
                if ln.startswith('Assigned'):
                    #print ln.split('a')
                    assigned_num = float(ln.split('\t')[1].strip())
                if ln.startswith('Unassigned_NoFeatures'):
                    #print ln.split('a')
                    unassigned_num = float(ln.split('\t')[1].strip())
            feature_map_rate = assigned_num / (assigned_num + unassigned_num)
            new_ln = gsm_id + '!' + str(round(float(map_unique_map_list[0]), 3)*100) + '!'+ str(round(float(map_unique_map_list[1]), 3)*100) + '!'+ str(round(float(feature_map_rate), 3)*100) + '\n'
            c.write(new_ln)
        except:
            print gsm_id

c.close()