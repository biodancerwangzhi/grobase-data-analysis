# -*- coding:utf-8 -*-
'''
    compare promoter and gene body density, to find out pausing genes
'''

import os
os.system('mkdir ../results/promoter_gene_body_featureCount/')
# protein.1k_tss_1k
promoter_annotation_file = '/data1/szx/metadata/gene_annotation/tss_updown_1k/{{updown_1k_tss}}'
# protein.tss_1k_tes_1k
gene_body_annotation_file = '/data1/szx/metadata/gene_annotation/tss_1k_tes_1k/{{tss_1k_tes_1k}}'

pair_end_sam = []

sam_list = os.listdir('../results/merge_bam') 
for bam_file in sam_list:
    if bam_file.endswith('.bam'):
        popen_content = os.popen('samtools view -H ../results/merge_bam/%s'%bam_file).read()
        if '../results/clean_reads/pair_end' in popen_content:
            pair_end_sam.append(bam_file)

for sam_file in sam_list:
    if '.bai' in sam_file:
        continue
    if sam_file in pair_end_sam:
        os.system('featureCounts -p -T 40 -t gene -s 1 -O -a %s \
-o ../results/promoter_gene_body_featureCount/%s_promoter.txt -M ../results/merge_bam/%s'%(promoter_annotation_file, sam_file.strip('.bam'), sam_file))    
    else:
        os.system('featureCounts -T 40 -t gene -s 1 -O -a %s \
-o ../results/promoter_gene_body_featureCount/%s_promoter.txt -M ../results/merge_bam/%s'%(promoter_annotation_file, sam_file.strip('.bam'), sam_file))


for sam_file in sam_list:
    if '.bai' in sam_file:
        continue
    if sam_file in pair_end_sam:
        os.system('featureCounts -p -T 40 -t gene -s 1 -O -a %s \
-o ../results/promoter_gene_body_featureCount/%s_gene_body.txt -M ../results/merge_bam/%s'%(gene_body_annotation_file, sam_file.strip('.bam'), sam_file))    
    else:
        os.system('featureCounts -T 40 -t gene -s 1 -O -a %s \
-o ../results/promoter_gene_body_featureCount/%s_gene_body.txt -M ../results/merge_bam/%s'%(gene_body_annotation_file, sam_file.strip('.bam'), sam_file))
