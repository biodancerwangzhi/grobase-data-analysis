import os
dir_specie = './'

c = open('./pair_end.txt', 'w')
specie_dir = ['ara', 'cae', 'dro', 'human', 'mac', 'man','mouse','pan','pla', 'rat', 'sac','sch', 'zea']
for specie_file in specie_dir:
    dir_bam = dir_specie + specie_file + '/results/merge_bam/'
    for bam_file in os.listdir(dir_bam):
        if bam_file.endswith('.bai'):
            continue
        else:
            gsm_id = bam_file.strip('.bam')
            popen_content = os.popen('samtools view -H %s'%dir_bam + bam_file).read()
            if '../results/clean_reads/pair_end' in popen_content:
                c.write(gsm_id +'\t' + 'true\n')
            else:
                c.write(gsm_id +'\t' + 'false\n')

c.close()