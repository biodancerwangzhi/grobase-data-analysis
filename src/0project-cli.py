# -*- coding:utf-8 -*-
import os
import csv
import re

table_header_list = ['gsm_id', 'srx_id', 'srr_id', 'gse', 'organism', 'seq_type', 
'source', 'tissue', 'bio_type', 'paper_title', 'pmid', 'organization', 'date_time']
dict_organism_abbr = {
    'Arabidopsis thalian': 'ara', 
    'Caenorhabditis elegans': 'cae', 
    'Drosophila melanogaster': 'dro', 
    'Homo sapiens': 'human', 
    'Macaca mul': 'mac',
    'Manihot esculen': 'man',
    'Mus musculus': 'mouse', 
    'Pan troglodytes': 'pan',
    'Plasmodium falciparum': 'pla', 
    'Rattus norvegicus': 'rat', 
    'Saccharomyces cerevisiae': 'sac',
    'Schizosaccharomyces pombe': 'sch',
    'Zea mays': 'zea'
}

dict_abbr_index = {
    'ara':'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa_bwa_index', 
    'cae':'Caenorhabditis_elegans.WBcel235.dna.toplevel.fa_bwa_index', 
    'dro':'Drosophila_melanogaster.BDGP6.dna.toplevel.fa_bwa_index', 
    'human':'Homo_sapiens.GRCh38.dna.primary_assembly.fa_bwa_index', 
    'mac':'Macaca_mulatta.Mmul_8.0.1.dna.toplevel.fa.gz_bwa_index',
    'man':'Manihot_esculenta.Manihot_esculenta_v6.dna.toplevel.fa_bwa_index',
    'mouse':'Mus_musculus.GRCm38.dna.primary_assembly.fa_bwa_index', 
    'pan':'Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa.gz_bwa_index',
    'pla':'Plasmodium_falciparum.EPr1.dna.toplevel.fa_bwa_index', 
    'rat':'Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa_bwa_index', 
    'sac':'Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa_bwa_index',
    'sch':'Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz_bwa_index',
    'zea':'Zea_mays.B73_RefGen_v4.dna.toplevel.fa_bwa_index'
}



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

    ##### 后面还有些文件需要加到src模板，并且动态修改它们的注释文件名。  


    

