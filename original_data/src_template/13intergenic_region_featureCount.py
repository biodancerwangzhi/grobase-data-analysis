import os
os.system('mkdir ../results/intergenic_region_featureCount/')
# intergenic_region_file
intergenic_region_file = '/data1/szx/metadata/gene_annotation/intergenic_region/{{intergenic_region}}'

pair_end_sam = []
sam_list = os.listdir('../results/merge_bam') 
for bam_file in sam_list:
    if bam_file.endswith('.bam'):
        popen_content = os.popen('samtools view -H ../results/merge_bam/%s'%bam_file).read()
        if '../results/clean_reads/pair_end' in popen_content:
            pair_end_sam.append(bam_file)

sam_list = os.listdir('../results/merge_bam') 

for sam_file in sam_list:
    if '.bai' in sam_file:
        continue
    if sam_file in pair_end_sam:
        os.system('featureCounts -p -T 40 -t gene -s 1 -a %s \
-o ../results/intergenic_region_featureCount/%s_count.txt -M ../results/merge_bam/%s'%(intergenic_region_file, sam_file.strip('.bam'), sam_file))    
    else:
        os.system('featureCounts -T 40 -t gene -s 1 -a %s \
-o ../results/intergenic_region_featureCount/%s_count.txt -M ../results/merge_bam/%s'%(intergenic_region_file, sam_file.strip('.bam'), sam_file))