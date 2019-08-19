# -*- coding:utf-8 -*-

'''
    In theory, gro-seq exon and intron density are follow average distribution
    using gene region density to represent intron, generate exon/intron as qc
'''

import os
os.system('mkdir ../results/exons_gene_featureCount')
annotation_file = '/data1/szx/metadata/gene_annotation/original_gtf/{{annotation_gtf}}'

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
        os.system('featureCounts -p -T 40 -t gene -s 1 -a %s -O \
-o ../results/exons_gene_featureCount/%s_gene_count.txt -M ../results/merge_bam/%s'%(annotation_file, sam_file.strip('.bam'), sam_file))    
    else:
        os.system('featureCounts -T 40 -t gene -s 1 -a %s -O \
-o ../results/exons_gene_featureCount/%s_gene_count.txt -M ../results/merge_bam/%s'%(annotation_file, sam_file.strip('.bam'), sam_file))


for sam_file in sam_list:
    if '.bai' in sam_file:
        continue
    if sam_file in pair_end_sam:
        os.system('featureCounts -p -T 40 -t exon -s 1 -a %s -O \
-o ../results/exons_gene_featureCount/%s_exon_count.txt -M ../results/merge_bam/%s'%(annotation_file, sam_file.strip('.bam'), sam_file))    
    else:
        os.system('featureCounts -T 40 -t exon -s 1 -a %s -O \
-o ../results/exons_gene_featureCount/%s_exon_count.txt -M ../results/merge_bam/%s'%(annotation_file, sam_file.strip('.bam'), sam_file))
