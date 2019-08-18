import os

align_file_list = os.listdir('../results/bwa_align')
os.system('mkdir ../results/bwa_align_sorted_bam/')
result_list = [x.split('_')[1].strip('.bam') for x in os.listdir('../results/bwa_align_sorted_bam')]
align_file_list = [x for x in align_file_list if x.split('.')[0] not in result_list]

print len(align_file_list)

for align_file in align_file_list:
    align_file = align_file.strip('\n')
    os.system('samtools sort -m 3G --threads 40 ../results/bwa_align/%s \
-o ../results/bwa_align_sorted_bam/sorted_%s.bam'%(align_file,align_file.split('.')[0]))
