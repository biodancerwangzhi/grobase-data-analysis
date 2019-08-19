# -*- coding:utf-8 -*-
'''
    比较promoter区和gene体的密度，从而计算是否暂停。
'''

import os
os.system('mkdir ../results/promoter_gene_body_featureCount/')
# protein.1k_tss_1k
promoter_annotation_file = '/data1/szx/metadata/gene_annotation/human/Homo_sapiens.protein.1k_tss_1k.gtf'
# protein.tss_1k_tes
gene_body_annotation_file = '/data1/szx/metadata/gene_annotation/human/Homo_sapiens.protein.tss_1k_tes.gtf'

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
