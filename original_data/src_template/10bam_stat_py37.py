#  run 'nohup python 10bam_stat_py37.py 1>../logs/10bam_stat.log 2>&1 &'  to get log file
import os

list_bam_bai = os.listdir('../results/merge_bam')
list_bam = [i for i in list_bam_bai if i.endswith('.bam')]
for str_bam in list_bam:
	os.system('echo "############ %s\n"'%str_bam) 
	os.system('bam_stat.py -i ../results/merge_bam/%s'%str_bam)
