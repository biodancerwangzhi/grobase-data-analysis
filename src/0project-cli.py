# -*- coding:utf-8 -*-
import os
import csv
import re
import json

f=open('./common.json')
content = f.read()
j=json.loads(content)

table_header_list = j['table_header_list']
dict_organism_abbr = j['dict_organism_abbr']
dict_abbr_index = j['dict_abbr_index']
dict_abbr_annotation = j['dict_abbr_annotation']

# 1 先创建所有物种的文件夹，和下级的文件夹
gsm_detail_dir = '../original_data/gsm_detail_results.csv'
set_specie_dir = set()
set_organism = set()
dict_gsm_srr_list = {}
for gsm_detail_ln in open(gsm_detail_dir):
    gsm_field_list = gsm_detail_ln.split('!')
    specie_name = gsm_field_list[4]
    specie_abbr = dict_organism_abbr[specie_name]
    set_organism.add(specie_name)
    set_specie_dir.add(specie_abbr)
    gsm_id = gsm_field_list[0]
    srr_field = gsm_field_list[2]
    srr_id_list = []
    if '##' in srr_field:
        for srr_id in srr_field.split('##'):
            srr_id = srr_id.strip('#').strip('"')
            srr_id_list.append(srr_id)
    else:
        srr_id_list.append(srr_field.strip('#').strip('"'))
    dict_gsm_srr_list.setdefault(gsm_id, srr_id_list)

os.system('mkdir ../results/')
os.system('mkdir ../results/bulk_specie_results/')
for specie_dir in set_specie_dir:
    
    os.system('mkdir ../results/bulk_specie_results/%s'%specie_dir)
    os.system('mkdir ../results/bulk_specie_results/%s/original_data'%specie_dir)
    os.system('mkdir ../results/bulk_specie_results/%s/original_data/original_fastq'%specie_dir)
    os.system('mkdir ../results/bulk_specie_results/%s/results'%specie_dir)
    os.system('mkdir ../results/bulk_specie_results/%s/logs'%specie_dir)
    os.system('mkdir ../results/bulk_specie_results/%s/src'%specie_dir)
    
    # 2 拷贝src template
    os.system('cp -r -f ../original_data/src_template/* ../results/bulk_specie_results/%s/src'%specie_dir)
    # 3 改变模板中的注释文件  
    align_bwa = '../results/bulk_specie_results/%s/src/5align_bwa.py'%specie_dir
    strFile = open(align_bwa).read().replace('{{bwa_index}}', dict_abbr_index[specie_dir])
    c = open(align_bwa, 'w')
    c.write(strFile)
    c.close()

    gene_exon_feature_count = '../results/bulk_specie_results/%s/src/11exons_gene_featureCount.py'%specie_dir
    strFile = open(gene_exon_feature_count).read().replace('{{annotation_gtf}}', dict_abbr_annotation[specie_dir])
    c = open(gene_exon_feature_count, 'w')
    c.write(strFile)
    c.close()

    ##### 后面还有些文件需要加到src模板，并且动态修改它们的注释文件名。  



# 4 将每个物种的gsm 和 srr id 放进每个物种文件夹的original_data。  
# 将相应的fastq mv进相应的目录中。这步不要急，先print一下，然后再mv。等其他的都运行完再执行。

dict_specie_gsm_list = {}
dict_gsm_srr = {}
dict_gsm_srr_list = {}
for gsm_detail_ln in open(gsm_detail_dir):
    gsm_field_list = gsm_detail_ln.split('!')
    specie_name = gsm_field_list[4]
    gsm_id = gsm_field_list[0].strip('"').strip()
    srr_id = gsm_field_list[2]
    dict_gsm_srr.setdefault(gsm_id, srr_id)
    if specie_name in dict_specie_gsm_list.keys():
        dict_specie_gsm_list[specie_name].append(gsm_id)
    else:
        dict_specie_gsm_list.setdefault(specie_name, [gsm_id])

#print dict_gsm_srr
for specie_name in dict_specie_gsm_list.keys():
    specie_abbr = dict_organism_abbr[specie_name]
    gsm_id_list = dict_specie_gsm_list[specie_name]
    # print gsm_id_list
    gsm_txt_file = open('../results/bulk_specie_results/%s/original_data/gsm.txt'%specie_abbr, 'w')
    srr_txt_file = open('../results/bulk_specie_results/%s/original_data/srr.txt'%specie_abbr, 'w')
    for gsm_id in gsm_id_list:

        gsm_id=gsm_id.strip('"').strip()
        try:
            srr_id = dict_gsm_srr[gsm_id]
        except:
            print gsm_id
        gsm_txt_file.write(gsm_id + '\n')
        srr_txt_file.write(srr_id + '\n')

        #os.system('mv ../original_data/original_fastq/%s* ../results/bulk_specie_results/%s/original_data/original_fastq/'%(srr_id, specie_abbr))
    gsm_txt_file.close()
    srr_txt_file.close()