import os
import multiprocessing

def run_bamindex(bam_dir, bam_set, n):
    name = multiprocessing.current_process().name
    print name,'starting'
    for align_file_name in bam_set:
        os.system('samtools index %s%s \
    %s%s.bai'%(bam_dir, align_file_name,bam_dir, align_file_name))        
    return


def run_bamindex_with_multithreads(bam_dir, align_file_list):
    sample_num = len(align_file_list)
    if sample_num < 40:
        thread_num = sample_num    
        each_group_num = 1 
    else:
        thread_num = 40
        each_group_num = sample_num/thread_num +1        
    processes_pool = list()

    for n in range(thread_num):
        try:
            split_bam = align_file_list[n*each_group_num: n*each_group_num + each_group_num]
        except:
            split_bam = align_file_list[n*each_group_num:]
            
        p = multiprocessing.Process(target = run_bamindex, args=(bam_dir, split_bam,n,))
        #p = multiprocessing.Process(target = print_gsmid_list, args=(split_bam,n,))
        p.start()
        processes_pool.append(p)
    
    for p in processes_pool:
        p.join()
    print 'Process end.'

bam_dir = '../results/merge_bam/'
align_file_list = os.listdir('../results/merge_bam')
run_bamindex_with_multithreads(bam_dir, align_file_list)