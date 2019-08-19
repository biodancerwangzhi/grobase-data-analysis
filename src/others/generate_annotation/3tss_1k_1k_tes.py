import os

from lib import generate_region_from_gtf


output_dir = '/data1/szx/metadata/gene_annotation/tss_1k_tes_1k/'
os.system('mkdir %s'%output_dir)
dir_protein_gene = '/data1/szx/metadata/gene_annotation/protein_gene/'

condition_dict={
    'tss_tes': ['tss', 'tes'],
    'left_right_length': [1000, 1000]
}

for file in os.listdir(dir_protein_gene):
    file_dir = dir_protein_gene + file
    output = output_dir + file.replace('.gtf', '__tss_1k_tes_1k.gtf')
    generate_region_from_gtf(file_dir, output, condition_dict)