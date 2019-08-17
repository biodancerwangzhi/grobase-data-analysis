import os

os.system('mkdir ../results/merge_bam')
a_list = open('../original_data/gsm.txt').readlines()

b_list = open('../original_data/srr.txt').readlines()

for i, j in zip(a_list, b_list):
    if j == '\r\n':
        print i, 'not have srr'
    else:
        if '##' in j:
            split_j = j.strip('\r\n').split('##')
            new_command = 'samtools merge ../results/merge_bam/%s.bam'%i.strip('\r\n')
            for split_j_item in split_j:
                new_command += ' ../results/bwa_align_sorted_bam/sorted_%s.bam'%split_j_item.strip('#').strip('"')
            #print '%s'%new_command
            os.system('%s'%new_command)    
        else:
            os.system('mv ../results/bwa_align_sorted_bam/sorted_%s.bam ../results/merge_bam/%s.bam'%(j.strip('\r\n'),i.strip('\r\n')))