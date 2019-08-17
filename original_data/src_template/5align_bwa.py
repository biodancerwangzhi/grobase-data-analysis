import os

trim_fastq_list=os.listdir('../results/clean_reads/')
finished_fastq_list = [x.strip('.sam') for x in os.listdir('../results/bwa_align')]
# os.system('mkdir ../results/bwa_align')
try:
    trim_fastq_list.remove('pair_end')
except:
    pass
for fastq_name in trim_fastq_list:
    if fastq_name not in finished_fastq_list:
        os.system('bwa mem /data1/szx/metadata/index/{{bwa_index}}/bwa_index -t 40 \
../results/clean_reads/%s > ../results/bwa_align/%s.sam'%(fastq_name,fastq_name.strip('.fastq')))

set_pair_file_name = set()
for pair_file in os.listdir('../results/clean_reads/pair_end'):
    set_pair_file_name.add(pair_file.split('_')[0])
for pair_fastq__file_name in set_pair_file_name:
    os.system('bwa mem /data1/szx/metadata/index/{{bwa_index}}/bwa_index -t 40 \
../results/clean_reads/pair_end/%s_1_trim_paired.fastq.gz ../results/clean_reads/pair_end/%s_2_trim_paired.fastq.gz \
> ../results/bwa_align/%s.fastq.gz.sam'%(pair_fastq__file_name,pair_fastq__file_name,pair_fastq__file_name))
