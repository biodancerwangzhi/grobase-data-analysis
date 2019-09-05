# -*- coding:utf-8 -*-
import os

from lib import generate_slide_windows
from lib import generate_slide_windows_with_same_bin_number

#generate_slide_windows('Homo_sapiens.GRCh38.93_gene_protein_coding_tss_1k.gtf', 'Homo_sapiens.GRCh38.93_gene_protein_coding_tss_1k_windows50_slide5.gtf', 50, 5)            
# generate_slide_windows_with_same_bin_number
annotation_dir = '/data1/szx/metadata/gene_annotation/'
output_dir = annotation_dir + 'pause_gene_windows/'
os.system('mkdir %s'%output_dir)

dir_protein_gene = '/data1/szx/metadata/gene_annotation/protein_gene/'
for file in os.listdir(dir_protein_gene):
    file_dir = dir_protein_gene + file
    promoter_windows_output = output_dir + file.replace('.gtf', '__promoter_windows.gtf')
    geneBody_windows_output = output_dir + file.replace('.gtf', '__geneBody_windows.gtf')
    # promoter windows
    generate_slide_windows_with_same_bin_number(annotation_dir + 'tss_updown_1k/' + file.replace('.gtf', '__1k_tss_1k.gtf'), promoter_windows_output, 50, 100)
    # gene body 小窗
    generate_slide_windows_with_same_bin_number(annotation_dir + 'tss_1k_tes_1k/' + file.replace('.gtf', '__tss_1k_tes_1k.gtf'), geneBody_windows_output, 1000, 390)