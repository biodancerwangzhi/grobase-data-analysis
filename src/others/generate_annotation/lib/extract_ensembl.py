# -*- coding: utf-8 -*-
#########################
# Ensembl annotation needs to be divided according to different analysis pipelines.
# feature was extracted based on gene or exon, protein or non-conding.
# gtf can be trans to bed, and the chr from a number to chr+number, gene name become gene_biotype+'_'+gene_name, score become gene length.
#########################
import os

# file_ensembl_annotation = '/stornext/HPCScratch/home/wang.zhi/reference/gene/gene_annotation/Homo_sapiens.GRCh38.93.gtf'

'''
# 从gene，transcript，exon，intron这些feature中选择分析内容属于基因的哪部分。  
def extract_genefeature(gtf_file, gene_feature_list, output_file):
    c = open(output_file, 'w')
    for ln_gtf in open(gtf_file):
        if ln_gtf[0] == '#':
            continue
        if ln_gtf.split('\t')[2] in gene_feature_list:
            c.write(ln_gtf)
    c.close()
#extract_genefeature(file_ensembl_annotation, ['gene'])

# 选择像蛋白编码基因，非编码的基因分类。  
def extract_genetype(gtf_file, gene_type_list, output_file):
    c = open(output_file, 'w')
    for ln_gtf in open(gtf_file):
        if ln_gtf[0] == '#':
            continue   
        if ln_gtf.split('gene_biotype "')[1].split('"')[0] in gene_type_list:
            c.write(ln_gtf)
    c.close()
#extract_genetype('Homo_sapiens.GRCh38.93_gene.gtf', ['protein_coding'], output_file) 
#extract_genefeature(file_ensembl_annotation, ['gene', 'exon']) 
#extract_genetype('Homo_sapiens.GRCh38.93_gene_exon.gtf', ['protein_coding'])
'''

# 上面代码可以抽象成，任何抽取feature的方案，都是基于提取gtf的ln，然后根据条件提取合适的ln。  
# 经过这一抽象，条件就可以轻易组合如抽取蛋白编码基因的exon。 


def condition_filter(ln_gtf, condition_dict):
    if 'gene_feature_list' in condition_dict.keys():
        if ln_gtf.split('\t')[2] not in condition_dict['gene_feature_list']:
            return False
    
    if 'gene_type_list' in condition_dict.keys():
        if ln_gtf.split('gene_biotype "')[1].split('"')[0] not in condition_dict['gene_type_list']:
            return False
    return ln_gtf

def extract_condition_ln(gtf_file, output_file, condition_dict):
    c = open(output_file, 'w')
    for ln_gtf in open(gtf_file):
        if ln_gtf[0] == '#':
            continue   
        ln_gtf = condition_filter(ln_gtf, condition_dict)
        if ln_gtf:
            c.write(ln_gtf)
    c.close()

# example1: all gtf_protein_coding 
# condtion_dict = {'gene_type_list': ['protein_coding']}
# extract_condition_ln('Homo.gtf', 'Homo__protein.gtf', condition_dict)

# example 2: gtf_protein_coding exons
# condtion_dict = {'gene_type_list': ['protein_coding'], 'gene_feature_list': ['exon']}
# extract_condition_ln('Homo.gtf' 'Homo__protein.gtf', condition_dict,)
'''
#an application to get unreplicated sorted exon based on protein coding gene
list_gtf = open('Homo_sapiens.GRCh38.93_gene_exon_protein_coding.gtf').readlines()
dic_gene_exon = {}
for ln_gtf in list_gtf:
    list_ln_gtf = ln_gtf.split('\t')
    # only consider the gene only exist in one postion in genomes.
    str_gene_name = ln_gtf.split('gene_name "')[1].split('"')[0]
    if str_gene_name in dic_gene_exon.keys():
        dic_gene_exon[str_gene_name]['postion'].append([int(list_ln_gtf[3]), int(list_ln_gtf[4])])
        dic_gene_exon[str_gene_name]['type'].append(list_ln_gtf[2])
    else:
        dic_gene_exon.setdefault(str_gene_name, {'type': [list_ln_gtf[2]], 'chrom': list_ln_gtf[0], 'strand': list_ln_gtf[6], 'postion':[[int(list_ln_gtf[3]), int(list_ln_gtf[4])]]})


c1 = open('tmp_exon.bed', 'w')
d1 = open('tmp_gene.bed', 'w')
for gene_name, dict_gene_exon in zip(dic_gene_exon.keys(), dic_gene_exon.values()):
    str_chrom = 'chr' + dict_gene_exon['chrom']
    str_strand = dict_gene_exon['strand']
    gene_start = dict_gene_exon['postion'][0][0]
    gene_end = dict_gene_exon['postion'][0][1]  
    list_gene_postion =[y for x,y in zip(dict_gene_exon['type'], dict_gene_exon['postion']) if x == 'gene']
    if len(list_gene_postion) != 1:
        print gene_name, 'exists on genomics more than once'
        continue
    else:
        list_gene_postion = [y for x,y in zip(dict_gene_exon['type'], dict_gene_exon['postion']) if x == 'gene'][0]
        int_gene_length = int(list_gene_postion[1]) - int(list_gene_postion[0])
        dict_gene_exon['type'].remove('gene')                
        dict_gene_exon['postion'].remove(list_gene_postion)
        ln_bd_gene = str_chrom + '\t' + str(list_gene_postion[0]) + '\t' + str(list_gene_postion[1]) + '\t' + gene_name + '\t'+str(int_gene_length)+'\t' + str_strand + '\n'
        d1.write(ln_bd_gene)
    for exon_postion in sorted(dict_gene_exon['postion']):
        exon_start = exon_postion[0]
        exon_end = exon_postion[1]
        ln_bd_exon = str_chrom + '\t' + str(exon_start) + '\t' + str(exon_end) + '\t' + gene_name + '\t' +str(int_gene_length) + '\t' + str_strand + '\n'
        c1.write(ln_bd_exon)

d1.close()
c1.close() 


os.system('bedtools sort -i tmp_exon.bed > tmp_sort_exon.bed')
os.system('bedtools merge -i tmp_sort_exon.bed -c 4 -o distinct> tmp_merge_sort_exon.bed')

c = open('./merge_sort_exon_without_intersect.bed', 'w')
dic_gene_strand = {}
for ln_gene in open('./tmp_gene.bed'):
	gene_name = ln_gene.split('\t')[0]
	strand = ln_gene.split('\t')[-1].strip('\n')
	dic_gene_strand.setdefault(gene_name, strand)
	
for ln_exon in open('tmp_merge_sort_exon.bed'):
	exon_label = ln_exon.split('\t')[-1].strip('\n')
	set_exon = set(exon_label.split(','))
	if len(set_exon) == 1:
		str_exon = list(set_exon)[0]
		list_ln = ln_exon.split('\t')
		new_ln = list_ln[0]+'\t'+list_ln[1]+'\t'+list_ln[2]+'\t'+str_exon+'\t'+'0'+'\t'+dic_gene_strand[list_ln[0]]+'\n'
		c.write(new_ln)
c.close()
os.system('rm ./tmp*')
'''

