import os
import multiprocessing

def run_trimmomatic_se(fastq_dir,out_dir, gsm_id_set, n):
    name = multiprocessing.current_process().name
    print name,'starting'
    for fastq_file in gsm_id_set:
        os.system('java -Xmx2G -jar /home/szx/miniconda2/share/trimmomatic-0.38-1/trimmomatic.jar SE %s%s %s%s LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30 -threads 10'%(fastq_dir,fastq_file, out_dir,fastq_file))    
    return

def run_trimmomatic_with_multithreads_se(fastq_dir, out_dir, gsm_id_list):
    sample_num = len(gsm_id_list)
    if sample_num < 30:
        thread_num = sample_num    
        each_group_num = 1 
    else:
        thread_num = 30
        each_group_num = sample_num/thread_num +1        
    processes_pool = list()
    for n in range(thread_num):
        try:
            split_bam = gsm_id_list[n*each_group_num: n*each_group_num + each_group_num]
        except:
            split_bam = gsm_id_list[n*each_group_num:]
            
        p = multiprocessing.Process(target = run_trimmomatic_se, args=(fastq_dir, out_dir, split_bam,n,))
        #p = multiprocessing.Process(target = print_gsmid_list, args=(split_bam,n,))
        p.start()
        processes_pool.append(p)
    
    for p in processes_pool:
        p.join()
    print 'Process end.'


def main():
    fastq_dir='../results/reads_after_trim_polya/'
    out_dir='../results/clean_reads/'

    fastq_file_list=os.listdir(fastq_dir)

    os.system('mkdir %s'%out_dir)
    os.system('mkdir  %spair_end'%fastq_dir)
    os.system('mkdir  %spair_end'%out_dir)

    list_se_fastq_file = []
    for fastq_file in fastq_file_list:
        fastq_file=fastq_file.strip('\n')
        file_name=fastq_file.split('.')[0]
        if  file_name[-2:] in ['_1','_2']:
            os.system('mv %s%s %spair_end/%s'%(fastq_dir,fastq_file,fastq_dir,fastq_file))
            #print 'double_end  ',file_name    
        else:   
            list_se_fastq_file.append(fastq_file)

    if list_se_fastq_file != []:
        run_trimmomatic_with_multithreads_se(fastq_dir, out_dir, list_se_fastq_file)
                    
    pair_fastq_list=os.listdir('%spair_end'%fastq_dir)

    if pair_fastq_list == []:
        print 'no pair end reads'
        pass
    else:
        set_pair_set = set()
        for pair_fastq in pair_fastq_list:
            pair_file_name=pair_fastq.split('.')[0][:-2]
            print pair_file_name
            set_pair_set.add(pair_file_name)
        for pair_file_file_name in set_pair_set:
            os.system('java -Xmx100G -jar /home/szx/miniconda2/share/trimmomatic-0.38-1/trimmomatic.jar PE -threads 35 \
            ../results/reads_after_trim_polya/pair_end/%s_1.fastq.gz ../results/reads_after_trim_polya/pair_end/%s_2.fastq.gz \
            ../results/clean_reads/pair_end/%s_1_trim_paired.fastq.gz ../results/clean_reads/pair_end/%s_1_trim_unpaired.fastq.gz \
            ../results/clean_reads/pair_end/%s_2_trim_paired.fastq.gz ../results/clean_reads/pair_end/%s_2_trim_unpaired.fastq.gz \
            LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30 \
            '%(pair_file_file_name,pair_file_file_name,pair_file_file_name,pair_file_file_name,pair_file_file_name,pair_file_file_name))       		

if __name__ == '__main__':
    main()