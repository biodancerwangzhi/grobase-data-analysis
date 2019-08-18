# -*- coding: utf-8 -*-
import os, json

f=open('./common.json')
content = f.read()
j=json.loads(content)

dict_abbr_index = j['dict_abbr_index']
set_organism = dict_abbr_index.values()
dict_abbr_annotation = j['dict_abbr_annotation']
set_annotation = dict_abbr_annotation.values()

# 1 判断是否有基因组， index，annotation。
# 当没有的时候输出一下，后面再手工下载一下。
genome_dir_list = os.listdir('/data1/szx/metadata/genome')
index_name_list = os.listdir('/data1/szx/metadata/index')
annotation_list = os.listdir('/data1/szx/metadata/gene_annotation/original_gtf')

for organism_index in set_organism:
    organism_name = organism_index.replace('_bwa_index', '')
    if organism_index not in index_name_list:
        print organism_index + ' dont have genome index.'
    if organism_name not in genome_dir_list:
        print organism_name + ' dont have genome， we need download it.'

for annotation in set_annotation:
    if annotation not in annotation_list:
        print annotation + ' dont have annotation, download it.'

# 2 检测一下新下载的sra的文件名称和文件大小。  
# 下次下载新的sra时再用。

# 3 检测一下总表上的gsm和所有已经分析完的gsm的id比较。  
# have analysised 
analysised_gsm_list = []
# 现在的分析结果都存在这些文件夹下。
dir1 = '/tmp/szx/project/grobase/original_computation/primary_results/'
dir2 = '/data1/szx/'
dir_list1=os.listdir(dir1)
dir_list2=os.listdir(dir2)
    
for dir_name in dir_list1:
    print dir1 + dir_name
    if os.path.isdir(dir1 + dir_name):
        if os.path.exists(dir1 + dir_name + '/results/merge_bam'):
            analysised_gsm_list.extend(os.listdir(dir1 + dir_name + '/results/merge_bam'))

for dir_name in dir_list2:
    if os.path.isdir(dir2 + dir_name):
        if os.path.exists(dir2 + dir_name + '/results/merge_bam'):
            analysised_gsm_list.extend(os.listdir(dir2 + dir_name + '/results/merge_bam'))

analysised_gsm_list = [x for x in analysised_gsm_list if not x.endswith('.bai')]
print len(analysised_gsm_list)
# need to analysis
gsm_detail_dir = '../original_data/gsm_detail_results.csv'
# gsm_detail_dir = './gsm_detail_results.csv'

for gsm_detail_ln in open(gsm_detail_dir):
    gsm_field_list = gsm_detail_ln.split('!')
    specie_name = gsm_field_list[4]
    gsm_id = gsm_field_list[0]
    if gsm_id + '.bam' not in analysised_gsm_list:
        print specie_name, gsm_id
