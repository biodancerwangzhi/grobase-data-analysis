# -*- coding: utf-8 -*-
import os, json

f=open('./common.json')
content = f.read()
j=json.loads(content)

dict_abbr_index = j['dict_abbr_index']
set_organism = dict_abbr_index.values()
dict_abbr_annotation = j['dict_abbr_annotation']
set_annotation = dict_abbr_annotation.values()
dict_organism_abbr = j['dict_organism_abbr']

# 1 将fq文件移到对应的文件夹。  
# todo: 未来有fq文件时再写入


# 2 创造基因组并把faidx拷到文件夹，并重命名成一个通用名，以便后面代码不需要改
for (specie_name, specie_genome) in zip(dict_abbr_index.keys(), dict_abbr_index.values()):
    specie_genome = specie_genome.replace('_bwa_index', '')
    genome_dir = '/data1/szx/metadata/genome/' + specie_genome
    faidx_dir = '/data1/szx/metadata/genome/' + specie_genome + '.fai'
    os.system('cp %s ../results/bulk_specie_results/%s/src/genome.fai'%(faidx_dir,specie_name))