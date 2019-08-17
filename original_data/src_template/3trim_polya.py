import os
fastq_dir='../original_data/original_fastq/'
output_dir = '../results/reads_after_trim_polya/'
os.system('mkdir %s'%output_dir)

fastq_file_list=os.listdir(fastq_dir)

for i in fastq_file_list:
    fastq_file=fastq_dir+i
    output_file=output_dir+i
    os.system('cutadapt -j 30 -a "A{100}" -o %s %s'%(output_file,fastq_file))