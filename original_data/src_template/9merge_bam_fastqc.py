import os
import multiprocessing

def run_fastqc(fastq_dir, output_dir, fastq_set, n):
    name = multiprocessing.current_process().name
    print name,'starting'
    for fastq_file in fastq_set:
        os.system('fastqc -t 40 %s%s -o %s'%(fastq_dir, fastq_file, output_dir))
    return

def run_fastqc_with_multithreads(fastq_dir, output_dir, fastq_file_list):
    sample_num = len(fastq_file_list)
    if sample_num < 40:
        thread_num = sample_num    
        each_group_num = 1 
    else:
        thread_num = 40
        each_group_num = sample_num/thread_num +1        
    processes_pool = list()

    for n in range(thread_num):
        try:
            split_bam = fastq_file_list[n*each_group_num: n*each_group_num + each_group_num]
        except:
            split_bam = fastq_file_list[n*each_group_num:]
            
        p = multiprocessing.Process(target = run_fastqc, args=(fastq_dir, output_dir, split_bam,n,))
        #p = multiprocessing.Process(target = print_gsmid_list, args=(split_bam,n,))
        p.start()
        processes_pool.append(p)
    
    for p in processes_pool:
        p.join()
    print 'Process end.'

bam_dir='../results/merge_bam/'
output_dir = '../results/merge_bam_fastqc/'
os.system('mkdir %s'%output_dir)

bam_file_list=[x for x in os.listdir(bam_dir) if x.endswith('.bam')]
run_fastqc_with_multithreads(bam_dir, output_dir, bam_file_list)



