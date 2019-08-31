# -*- coding: utf-8 -*-
'''
配置分几个层面， 
基因范围
    tss附近还是test附近。
    上游长度和下游长度。 

有没有windows
    windows长度  
    slide长度
'''

condtion_dict = {
    'tss_tes': ['tss', 'tes'], 
    'upper_downstream': {
        'tss_upper':0,
        'tss_down': 0,
        'tes_upper': 0,
        'tes_down': 0
    }
}

# 将tss上下游范围， 抽象成左右点 和上下游长度
def parse_tss_tes_conditon(condition_dict):
    if condition_dict['tss_tes'] == ['tss']:
        left_point = 'tss'
        right_point = 'tss'      
    if condition_dict['tss_tes'] == ['tes']:
        left_point = 'tes'
        right_point = 'tes'    
    if condition_dict['tss_tes'] in [['tss', 'tes'], ['tes', 'tss']]:
        left_point = 'tss'
        right_point = 'tes'
    return left_point, right_point

def generate_region_from_gtf(file, output_file, condition_dict):  
    if 'tss_tes' in condition_dict.keys():   
        left_point, right_point= parse_tss_tes_conditon(condition_dict)
        [left_length, right_length] = condition_dict['left_right_length']  
    '''    
    if 'windows_length' in condition_dict.keys():
        windows_length = condition_dict['windows_length']
        slide_length = condition_dict['slide_length']
    '''
    c = open(output_file, 'w')
    for ln_gene in open(file):
        if ln_gene[0] == '#':
            continue
        list_gene = ln_gene.strip('\n').split('\t')
        strand = list_gene[6]
        start = list_gene[3]
        end = list_gene[4]
        length = int(end) - int(start)
        left = start
        right = end

        if 'tss_tes' in condition_dict.keys(): 
            length_threshold = 0
            # length_threshold is the parameter for  1k tss 1k, some gene less than 2kb.
            if 'length_threshold' in condition_dict.keys():   
                length_threshold = condition_dict['length_threshold']   
                if length < length_threshold:   
                    left_length = int(length/2) if(left_length>=0) else int(length/2)*-1
                    right_length = int(length/2) if (right_length>=0) else int(length/2)*-1

            if strand == '+':
                if left_point == 'tss':
                    left = int(start) + left_length
                else:
                    left = int(end) + left_length
                if right_point == 'tss': 
                    right = int(start) + right_length
                else:
                    right = int(end) + right_length
            else:
                # 左边的点用来确定左边，负链
                if left_point == 'tss':
                    left = int(start) - right_length
                else:
                    left = int(end) - right_length
                if right_point == 'tss':
                    right = int(start) - left_length
                else:
                    right = int(end) - left_length

        if left <= 0:
            left = 1

        new_ln = ln_gene.split('\tgene\t')[0]+'\tgene\t' + str(left)+'\t'+ str(right)+'\t.\t'+strand+'\t.\tgene_id' +ln_gene.split('.\tgene_id')[1]
        c.write(new_ln)
    c.close()                 

def generate_tss_region_from_gtf(file, output_file, distance=1000):
    c = open(output_file, 'w')
    for ln_gene in open(file):
        if ln_gene[0] == '#':
            continue
        list_gene = ln_gene.strip('\n').split('\t')
        strand = list_gene[6]
        start = list_gene[3]
        end = list_gene[4]
        if strand == '+':
            tss = int(start)
        else:
            tss = int(end)
        tss_left_1k = tss - distance
        tss_rigth_1k = tss + distance
        if tss_left_1k <0:
            continue         
        new_ln = ln_gene.split('\tgene\t')[0]+'\tgene\t' + str(tss_left_1k)+'\t'+ str(tss_rigth_1k)+'\t.\t'+strand+'\t.\tgene_id' +ln_gene.split('.\tgene_id')[1]
        c.write(new_ln)
    c.close()    

def generate_range_region_from_gtf(file, output_file, distance=1000):
    c = open(output_file, 'w')
    for ln_gene in open(file):
        if ln_gene[0] == '#':
            continue
        list_gene = ln_gene.strip('\n').split('\t')
        strand = list_gene[6]
        start = list_gene[3]
        end = list_gene[4]
        start_left_1k = int(start) - distance
        end_rigth_1k = int(end) + distance
        if start_left_1k <0:
            continue         
        new_ln = ln_gene.split('\tgene\t')[0]+'\tgene\t' + str(start_left_1k)+'\t'+ str(end_rigth_1k)+'\t.\t'+strand+'\t.\tgene_id' +ln_gene.split('.\tgene_id')[1]
        c.write(new_ln)
    c.close() 
#generate_tss_region_from_gtf('./Homo_sapiens.GRCh38.93_gene_protein_coding.gtf', 'Homo_sapiens.GRCh38.93_gene_protein_coding_tss_1k.gtf')        
#generate_tss_region_from_gtf('./Homo_sapiens.GRCh38.93_gene_protein_coding.gtf', 'Homo_sapiens.GRCh38.93_gene_protein_coding_tss_3k.gtf', 3000)        

def generate_slide_windows(file, output_file, bin_length, slide_length = 50):
    c = open(output_file, 'w')
    for ln_gene in open(file):
        if ln_gene[0] == '#':
            continue
        list_gene = ln_gene.strip('\n').split('\t')
        strand = list_gene[6]
        start = int(list_gene[3])
        end = int(list_gene[4]) 
        length = end - start
        bin_number = (length-bin_length)/slide_length
        list_type = ln_gene.split('gene_id')[1].split(';')

        for bin_num in range(bin_number):
            new_start = start + bin_num * slide_length
            new_end = new_start + bin_length
            new_ln = ln_gene.split('\tgene\t')[0]+'\tgene\t' + \
			str(new_start) +'\t'+ str(new_end) +'\t.\t'+strand+'\t.\tgene_id "' + \
			ln_gene.split('gene_id "')[1].split('"')[0]+'_'+str(bin_num)+'"; '
            for bio_type in list_type[1:-1]:
		new_ln = new_ln+bio_type+'; '
            new_ln = new_ln + '\n'
            c.write(new_ln)
    c.close()    
#generate_slide_windows('Homo_sapiens.GRCh38.93_gene_protein_coding_tss_3k.gtf', 'Homo_sapiens.GRCh38.93_gene_protein_coding_tss_3k_windows50_slide5.gtf', 50, 5)            
#generate_slide_windows('Homo_sapiens.GRCh38.93_gene_protein_coding_tss_1k.gtf', 'Homo_sapiens.GRCh38.93_gene_protein_coding_tss_1k_windows50_slide5.gtf', 50, 5)            

import math
def generate_slide_windows_with_same_bin_number(file, output_file, bin_length, bin_number):
    c = open(output_file, 'w')
    for ln_gene in open(file):
        if ln_gene[0] == '#':
            continue
        list_gene = ln_gene.strip('\n').split('\t')
        strand = list_gene[6]
        start = int(list_gene[3])
        end = int(list_gene[4]) 
        length = end - start
        if length <bin_length:
            bin_length = length-1
        slide_length = int(math.ceil((length-bin_length)/ float(bin_number)))
        list_type = ln_gene.split('gene_id')[1].split(';')

        for bin_num in range(bin_number):
            new_start = start + bin_num * slide_length
            new_end = new_start + bin_length
            new_ln = ln_gene.split('\tgene\t')[0]+'\tgene\t' + \
			str(new_start) +'\t'+ str(new_end) +'\t.\t'+strand+'\t.\tgene_id "' + \
			ln_gene.split('gene_id "')[1].split('"')[0]+'_'+str(bin_num)+'"; '
            for bio_type in list_type[1:-1]:
		new_ln = new_ln+bio_type+'; '
            new_ln = new_ln + '\n'
            c.write(new_ln)
    c.close()   		 

def generate_tss_1k_tes(file, output_file):
    c = open(output_file, 'w')
    for ln_gene in open(file):
        if ln_gene[0] == '#':
            continue
        list_gene = ln_gene.strip('\n').split('\t')
        strand = list_gene[6]
        start = int(list_gene[3])
        end = int(list_gene[4]) 
        length = end - start   
        list_type = ln_gene.split('gene_id')[1].split(';')
        if length > 3000:
            if strand == '+':
                tss = start
                tes = end
                tss_downstream_1k = tss + 1000
                # new_ln_tss = ln_gene.split('\tgene\t')[0]+'\tgene\t' + str(tss) +'\t'+ str(tss_downstream_1k) +'\t.\t'+strand+'\t.\tgene_id "' +ln_gene.split('gene_id "')[1].split('"')[0]+'_tss'+'"; gene_version '+ln_gene.split('"; gene_version ')[1]
                new_ln_tes = ln_gene.split('\tgene\t')[0]+'\tgene\t' + str(tss_downstream_1k + 1) +'\t'+ str(tes) +'\t.\t'+strand+'\t.\tgene_id "' + \
                    ln_gene.split('gene_id "')[1].split('"')[0]+'_tes'+'"; '                
                for bio_type in list_type[1:-1]:
                    new_ln_tes = new_ln_tes+bio_type+'; '
                new_ln_tes = new_ln_tes + '\n'
                #c.write(new_ln_tss)
                c.write(new_ln_tes)
            else:
                tss = end
                tes = start
                tss_upstream_1k =tss  -1000
                #new_ln_tss = ln_gene.split('\tgene\t')[0]+'\tgene\t' + str(tss_upstream_1k)+'\t'+ str(tss)  +'\t.\t'+strand+'\t.\tgene_id "' +ln_gene.split('gene_id "')[1].split('"')[0]+'_tss'+'"; gene_version '+ln_gene.split('"; gene_version ')[1]
                new_ln_tes = ln_gene.split('\tgene\t')[0]+'\tgene\t' + str(tes)+'\t'+ str(tss_upstream_1k - 1)  +'\t.\t'+strand+'\t.\tgene_id "' + \
                    ln_gene.split('gene_id "')[1].split('"')[0]+'_tes'+'"; '
                for bio_type in list_type[1:-1]:
                    new_ln_tes = new_ln_tes+bio_type+'; '
                new_ln_tes = new_ln_tes + '\n'
                c.write(new_ln_tes)
                #c.write(new_ln_tss)    
    c.close()
        
#generate_tss_1k_tes('./Homo_sapiens.GRCh38.93_gene_protein_coding.gtf', 'Homo_sapiens.GRCh38.93_gene_protein_coding_tss_1k_tes.gtf')
#generate_tss_1k_tes('./Homo_sapiens.GRCh38.93_gene.gtf', 'Homo_sapiens.GRCh38.93_gene_tss_1k_tes.gtf')