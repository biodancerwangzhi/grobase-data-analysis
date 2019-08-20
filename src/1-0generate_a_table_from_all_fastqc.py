import os
dir_specie = './'
c = open('./1fastqc_summary.csv', 'w')
first_ln = 'gsm_id!fastqc_pass_num!sequence_length!total_sequence!gc_rate\n'
c.write(first_ln)

list_specie_file = os.listdir(dir_specie)
list_specie_file = ['ara', 'cae', 'dro', 'human', 'mac', 'man','mouse','pan','pla', 'rat', 'sac','sch', 'zea']
for specie_file in list_specie_file:
    dir_specie_fastqc = dir_specie + specie_file + '/results/merge_bam_fastqc/'
    list_fastq_sample = os.listdir(dir_specie_fastqc)
    for sample_name in list_fastq_sample:
        if sample_name.endswith('.html'):
            gsm_id = sample_name.split('_')[0]
            dir_fastqc = dir_specie_fastqc + sample_name
            str_fastqc = open('%s'%dir_fastqc).read()

            fastqc_pass_num = str_fastqc.count('[PASS]')    
            sequence_length = str_fastqc.split('<td>Sequence length</td><td>')[1].split('</td>')[0] 
            total_sequence = str_fastqc.split('<td>Total Sequences</td><td>')[1].split('</td>')[0]     
            gc_rate = str_fastqc.split('<td>%GC</td><td>')[1].split('</td>')[0]      
            new_ln = gsm_id + '!' + str(fastqc_pass_num)+ '!' + sequence_length+'!'+total_sequence+'!'+gc_rate+'\n'
            c.write(new_ln)

c.close()