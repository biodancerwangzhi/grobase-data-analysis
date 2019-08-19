import os

import lib
annotation_dir = '/data1/szx/metadata/gene_annotation/original_gtf/'
output_dir = '/data1/szx/metadata/gene_annotation/protein_gene/'
os.system('mkdir %s'%output_dir)

condition_dict = {
    'gene_feature_list': ['gene'],
    'gene_type_list': ['protein_coding']
}

for file in os.listdir(annotation_dir):
    file_dir = annotation_dir + file
    output = output_dir + file.replace('.gtf', '__protein_coding__gene.gtf')
    lib.extract_condition_ln(file_dir, output, condition_dict)
