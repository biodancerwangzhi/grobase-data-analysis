import os
import pandas as pd

dir_gene_exon_feature_count = '../results/exons_gene_featureCount/'
dir_promoter_gene_body_feature_count = '../results/promoter_gene_body_featureCount/'
dir_promoter_gene_body_windows_feature_count = '../results/promoter_geneBody_windows_featureCount/'
dir_intergenic_feature_count = '../results/intergenic_region_featureCount/'

list_exon_gene_feature_count = [x for x in os.listdir(dir_gene_exon_feature_count) if x.endswith('.txt')]
list_exon_feature_count = [x for x in list_exon_gene_feature_count if 'exon' in x]
list_gene_feature_count = [x for x in list_exon_gene_feature_count if 'gene' in x]

list_promoter_gene_body_feature_count = [x for x in os.listdir(dir_promoter_gene_body_feature_count) if x.endswith('.txt')]
list_promoter_feature_count = [x for x in list_promoter_gene_body_feature_count if 'promoter' in x]
list_gene_body_feature_count = [x for x in list_promoter_gene_body_feature_count if 'gene_body' in x]

list_promoter_gene_body_windows_feature_count = [x for x in os.listdir(dir_promoter_gene_body_windows_feature_count) if x.endswith('.txt')]
list_promoter_windows_feature_count = [x for x in list_promoter_gene_body_windows_feature_count if 'promoter' in x]
list_gene_body_windows_feature_count = [x for x in list_promoter_gene_body_windows_feature_count if 'gene_body' in x]

list_intergenic_feature_count = [x for x in os.listdir(dir_intergenic_feature_count) if x.endswith('.txt')]

exon_ensembl_name_list = [y.split('\t')[0].strip('\n') for y in [x for x in open(dir_gene_exon_feature_count+ list_exon_feature_count[0])][2:]]
exon_ensembl_column = pd.Series(exon_ensembl_name_list, name= 'ensembl_id')
gene_ensembl_name_list = [y.split('\t')[0].strip('\n') for y in [x for x in open(dir_gene_exon_feature_count+ list_gene_feature_count[0])][2:]]
gene_ensembl_column = pd.Series(gene_ensembl_name_list, name= 'ensembl_id')

promoter_ensembl_name_list = [y.split('\t')[0].strip('\n').replace('_promoter', '') for y in [x for x in open(dir_promoter_gene_body_feature_count+ list_promoter_feature_count[0])][2:]]
promoter_ensembl_column = pd.Series(promoter_ensembl_name_list, name= 'ensembl_id')
gene_body_ensembl_name_list = [y.split('\t')[0].strip('\n').replace('_gene_body', '') for y in [x for x in open(dir_promoter_gene_body_feature_count+ list_gene_body_feature_count[0])][2:]]
gene_body_ensembl_column = pd.Series(gene_body_ensembl_name_list, name= 'ensembl_id')

promoter_windows_ensembl_name_list = [y.split('\t')[0].strip('\n').replace('_promoter', '') for y in [x for x in open(dir_promoter_gene_body_windows_feature_count+ list_promoter_windows_feature_count[0])][2:]]
promoter_windows_ensembl_column = pd.Series(promoter_windows_ensembl_name_list, name= 'ensembl_id')
gene_body_windows_ensembl_name_list = [y.split('\t')[0].strip('\n').replace('_gene_body', '') for y in [x for x in open(dir_promoter_gene_body_windows_feature_count+ list_gene_body_windows_feature_count[0])][2:]]
gene_body_windows_ensembl_column = pd.Series(gene_body_windows_ensembl_name_list, name= 'ensembl_id')

intergenic_ensembl_name_list = [y.split('\t')[0].strip('\n').replace('_intergenic', '') for y in [x for x in open(dir_intergenic_feature_count+ list_intergenic_feature_count[0])][2:]]
intergenic_ensembl_column = pd.Series(intergenic_ensembl_name_list, name= 'ensembl_id')

exon_length_list = [y.split('\t')[-2].strip('\n') for y in [x for x in open(dir_gene_exon_feature_count+ list_exon_feature_count[0])][2:]]
exon_length_column = pd.Series(exon_length_list, name= 'length')

gene_length_list = [y.split('\t')[-2].strip('\n') for y in [x for x in open(dir_gene_exon_feature_count+ list_gene_feature_count[0])][2:]]
gene_length_column = pd.Series(gene_length_list, name= 'length')

promoter_length_list = [y.split('\t')[-2].strip('\n').replace('_promoter', '') for y in [x for x in open(dir_promoter_gene_body_feature_count+ list_promoter_feature_count[0])][2:]]
promoter_length_column = pd.Series(promoter_length_list, name= 'length')

gene_body_length_list = [y.split('\t')[-2].strip('\n').replace('_gene_body', '') for y in [x for x in open(dir_promoter_gene_body_feature_count+ list_gene_body_feature_count[0])][2:]]
gene_body_length_column = pd.Series(gene_body_length_list, name= 'length')

promoter_windows_length_list = [y.split('\t')[-2].strip('\n').replace('_promoter', '') for y in [x for x in open(dir_promoter_gene_body_windows_feature_count+ list_promoter_windows_feature_count[0])][2:]]
promoter_windows_length_column = pd.Series(promoter_windows_length_list, name= 'length')

gene_body_windows_length_list = [y.split('\t')[-2].strip('\n').replace('_gene_body', '') for y in [x for x in open(dir_promoter_gene_body_windows_feature_count+ list_gene_body_windows_feature_count[0])][2:]]
gene_body_windows_length_column = pd.Series(gene_body_windows_length_list, name= 'length')

intergenic_length_list = [y.split('\t')[-2].strip('\n').replace('_intergenic', '') for y in [x for x in open(dir_intergenic_feature_count+ list_intergenic_feature_count[0])][2:]]
intergenic_length_column = pd.Series(intergenic_length_list, name= 'length')

annotation = '/data1/szx/metadata/gene_annotation/original_gtf/{{gtf_file}}'
dict_gsm_hugo_biotype = {}
for ln in open(annotation):
    if ln[0] != '#' and ln.split('\t')[2] == 'gene':
        try:
            ensembl_id = ln.split('gene_id "')[1].split('"')[0]
            hugo_id = ''
            biotype = ''
            biotype = ln.split('gene_biotype "')[1].split('"')[0]
            hugo_id = ln.split('gene_name "')[1].split('"')[0]
            dict_gsm_hugo_biotype.setdefault(ensembl_id, [hugo_id, biotype])
        except:
            dict_gsm_hugo_biotype.setdefault(ensembl_id, ['no', 'no'])
    
old_column_exon = [exon_ensembl_column]
hugo_column_exon = pd.Series([dict_gsm_hugo_biotype[x][0] for x in exon_ensembl_name_list], name= 'hugo_id')
biotype_column_exon = pd.Series([dict_gsm_hugo_biotype[x][1] for x in exon_ensembl_name_list], name= 'biotype')
old_column_exon.append(hugo_column_exon)
old_column_exon.append(exon_length_column)
old_column_exon.append(biotype_column_exon)

old_column_gene = [gene_ensembl_column]
hugo_column_gene = pd.Series([dict_gsm_hugo_biotype[x][0] for x in gene_ensembl_name_list], name= 'hugo_id')
biotype_column_gene = pd.Series([dict_gsm_hugo_biotype[x][1] for x in gene_ensembl_name_list], name= 'biotype')
old_column_gene.append(hugo_column_gene)
old_column_gene.append(gene_length_column)
old_column_gene.append(biotype_column_gene)

old_column_promoter = [promoter_ensembl_column]
hugo_column_promoter = pd.Series([dict_gsm_hugo_biotype[x][0] for x in promoter_ensembl_name_list], name= 'hugo_id')
biotype_column_promoter = pd.Series([dict_gsm_hugo_biotype[x][1] for x in promoter_ensembl_name_list], name= 'biotype')
old_column_promoter.append(promoter_length_column)
old_column_promoter.append(hugo_column_promoter)
old_column_promoter.append(biotype_column_promoter)

old_column_gene_body = [gene_body_ensembl_column]
hugo_column_gene_body = pd.Series([dict_gsm_hugo_biotype[x][0] for x in gene_body_ensembl_name_list], name= 'hugo_id')
biotype_column_gene_body = pd.Series([dict_gsm_hugo_biotype[x][1] for x in gene_body_ensembl_name_list], name= 'biotype')
old_column_gene_body.append(gene_body_length_column)
old_column_gene_body.append(hugo_column_gene_body)
old_column_gene_body.append(biotype_column_gene_body)

old_column_promoter_windows = [promoter_windows_ensembl_column]
old_column_promoter_windows.append(promoter_windows_length_column)

old_column_gene_body_windows = [gene_body_windows_ensembl_column]
old_column_gene_body_windows.append(gene_body_windows_length_column)

old_column_intergenic = [intergenic_ensembl_column]
old_column_intergenic.append(intergenic_length_column)

for gsm_feature_count in list_exon_feature_count:
    gsmId = gsm_feature_count.split('_')[0]
    exp_list = [int(y.split('\t')[-1].strip('\n')) for y in [x for x in open(dir_gene_exon_feature_count+ gsm_feature_count)][2:]]
    new_column = pd.Series(exp_list, name= gsmId)
    old_column_exon.append(new_column)

for gsm_feature_count in list_gene_feature_count:
    gsmId = gsm_feature_count.split('_')[0]
    exp_list = [int(y.split('\t')[-1].strip('\n')) for y in [x for x in open(dir_gene_exon_feature_count+ gsm_feature_count)][2:]]
    new_column = pd.Series(exp_list, name= gsmId)
    old_column_gene.append(new_column)

for gsm_feature_count in list_promoter_feature_count:
    gsmId = gsm_feature_count.split('_')[0]
    exp_list = [int(y.split('\t')[-1].strip('\n')) for y in [x for x in open(dir_promoter_gene_body_feature_count+ gsm_feature_count)][2:]]
    new_column = pd.Series(exp_list, name= gsmId)
    old_column_promoter.append(new_column)

for gsm_feature_count in list_gene_body_feature_count:
    gsmId = gsm_feature_count.split('_')[0]
    exp_list = [int(y.split('\t')[-1].strip('\n')) for y in [x for x in open(dir_promoter_gene_body_feature_count+ gsm_feature_count)][2:]]
    new_column = pd.Series(exp_list, name= gsmId)
    old_column_gene_body.append(new_column)

for gsm_feature_count in list_promoter_windows_feature_count:
    gsmId = gsm_feature_count.split('_')[0]
    exp_list = [int(y.split('\t')[-1].strip('\n')) for y in [x for x in open(dir_promoter_gene_body_windows_feature_count+ gsm_feature_count)][2:]]
    new_column = pd.Series(exp_list, name= gsmId)
    old_column_promoter_windows.append(new_column)

for gsm_feature_count in list_gene_body_windows_feature_count:
    gsmId = gsm_feature_count.split('_')[0]
    exp_list = [int(y.split('\t')[-1].strip('\n')) for y in [x for x in open(dir_promoter_gene_body_windows_feature_count+ gsm_feature_count)][2:]]
    new_column = pd.Series(exp_list, name= gsmId)
    old_column_gene_body_windows.append(new_column)

for gsm_feature_count in list_intergenic_feature_count:
    gsmId = gsm_feature_count.split('_')[0]
    exp_list = [int(y.split('\t')[-1].strip('\n')) for y in [x for x in open(dir_intergenic_feature_count+ gsm_feature_count)][2:]]
    new_column = pd.Series(exp_list, name= gsmId)
    old_column_intergenic.append(new_column)
  
save_exon = pd.concat(old_column_exon, axis=1)  
save_exon.to_csv('../results/exon_feature_count_matrix.csv', index=False, sep=',')

save_gene = pd.concat(old_column_gene, axis=1)  
save_gene.to_csv('../results/gene_feature_count_matrix.csv', index=False, sep=',')

save_promoter = pd.concat(old_column_promoter, axis=1)  
save_promoter.to_csv('../results/promoter_feature_count_matrix.csv', index=False, sep=',')

save_gene_body = pd.concat(old_column_gene_body, axis=1)  
save_gene_body.to_csv('../results/gene_body_feature_count_matrix.csv', index=False, sep=',')

save_promoter_windows = pd.concat(old_column_promoter_windows, axis=1)  
save_promoter_windows.to_csv('../results/promoter_windows_feature_count_matrix.csv', index=False, sep=',')

save_gene_body_windows = pd.concat(old_column_gene_body_windows, axis=1)  
save_gene_body_windows.to_csv('../results/gene_body_windows_feature_count_matrix.csv', index=False, sep=',')

save_intergenic = pd.concat(old_column_intergenic, axis=1)  
save_intergenic.to_csv('../results/intergenic_feature_count_matrix.csv', index=False, sep=',')