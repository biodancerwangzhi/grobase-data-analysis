# -*- coding: utf-8 -*-
import os, json

f=open('../../common.json')
content = f.read()
j=json.loads(content)

dict_annotation_genome = j['dict_annotation_genome']


def bd4tobd6(ln):
    	list_a = ln.split('\t')
	new_ln = list_a[0]+'\t'+list_a[1]+'\t'+list_a[2]+'\t'+'#'+'\t'+'0'+'\t'+list_a[3]
	return new_ln

def gtf2bd(gtf_file, file_bd):
    c = open(file_bd, 'w')
    for ln_gtf in open(gtf_file):
        if ln_gtf[0] == '#':
            continue
        list_gtf = ln_gtf.split('\t')
        gene_length = str(int(list_gtf[4])-int(list_gtf[3]))
        try:
            gene_name = ln_gtf.split('gene_name "')[1].split('"')[0]
        except:
            gene_name = ln_gtf.split('gene_id "')[1].split('"')[0]
        new_ln = 'chr'+list_gtf[0]+'\t'+list_gtf[3]+'\t'+list_gtf[4]+'\t'+ \
ln_gtf.split('gene_biotype "')[1].split('"')[0]+'_'+gene_name+'\t'+gene_length+'\t'+list_gtf[6]+'\n'
        c.write(new_ln)
    c.close()

def generate_chr_size(genome_size_file):
    dict_chr_size={}
    for j in open(genome_size_file):
        chrom = 'chr'+j.split('\t')[0]
        chrom_size = j.split('\t')[1]
        dict_chr_size.setdefault(chrom, chrom_size)
    return dict_chr_size

def bed2gtf(infile, outfile):  
    inf  = open(infile, 'r')
    outf = open(outfile,'w')
    cont = 1
    for linea in inf:
        linea_split = linea.split()
        chrom = linea_split[0].lstrip('chr')
        ini_pos = int(linea_split[1])
        fin_pos = int(linea_split[2])
	strand = linea_split[-1].strip('\n')
	name = linea_split[3]
        #peak = linea_split[3]

        #outf.write(chrom + "\tpfurio\tpeak\t" + str(ini_pos) + "\t" + str(fin_pos) + '\t.\t+\t.\tpeak_id "' + peak + '";\n')
        outf.write(chrom + "\thavana\tgene\t" + str(ini_pos) + "\t" + str(fin_pos) + '\t.\t'+strand+'\t.\tgene_id "' +
                name + '";\n')
        cont += 1
    inf.close()
    outf.close()   

# get genomic intergenic region by subtract gene +- 10k
# replace genomeSize for different species with samtools faidx ../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
intergenic_dir = '/data1/szx/metadata/gene_annotation/intergenic_gene/'
os.system('mkdir %s'%intergenic_dir)

dir_annotation = '/data1/szx/metadata/gene_annotation/original_gtf/'
list_specie_file = os.listdir(dir_annotation)
for specie_file in list_specie_file:
    dir_specie = dir_annotation + specie_file
    dir_specie_bd = './tmp_'+ specie_file.replace('.gtf', '.bed')
    gtf2bd(dir_specie, dir_specie_bd)

    genome_size_file = '/data1/szx/metadata/genome/' + dict_annotation_genome[specie_file].replace('.fa_bwa_index', '.fa.fai')
    dict_chr_size = generate_chr_size(genome_size_file)
    
    c = open('tmp_expand_10k.bed','w')
    for i in open(dir_specie_bd):
        list_si = i.split('\t')
        if int(list_si[1])-10000 > 0:
            new_ln = list_si[0]+'\t'+str(int(list_si[1])-5000)+'\t'+str(int(list_si[2])+5000)+'\t'+list_si[3]+'\t'+list_si[4]+'\t'+list_si[5]
        else:
            new_ln = i
        c.write(new_ln)
    c.close()

    os.system('bedtools sort -i tmp_expand_10k.bed > tmp_sort.bed')
    os.system('bedtools merge -s -c 6 -o distinct -i tmp_sort.bed > tmp_merge.bed')
    os.system('bedtools sort -i tmp_merge.bed > tmp_merge_sort.bed')
    list_chr = []
    for i in open('tmp_merge_sort.bed'):
        chrom = i.split('\t')[0]
        if chrom in list_chr:
            continue
        else:
            list_chr.append(chrom)

    c = open('tmp_sort_genome_size.bed', 'w')
    for k in list_chr:
        if k in dict_chr_size.keys():
            new_ln_postive = k +'\t'+'1'+'\t'+ dict_chr_size[k].strip('\n') +'\t'+ '+'+'\n'
            new_ln_negative = k +'\t'+'1'+'\t'+ dict_chr_size[k].strip('\n') +'\t'+'-'+'\n'
            c.write(new_ln_postive)
            c.write(new_ln_negative)
    c.close()

    c = open('tmp_merge_sort.bd6','w')
    for ln in open('tmp_merge_sort.bed'):
        new_ln = bd4tobd6(ln)
        c.write(new_ln)
    c.close()

    c = open('tmp_sort_genome_size.bd6','w')
    for ln in open('tmp_sort_genome_size.bed'):
        new_ln = bd4tobd6(ln)
        c.write(new_ln)
    c.close()

    os.system('bedtools subtract -s -a tmp_sort_genome_size.bd6 -b tmp_merge_sort.bd6 >tmp_subtract.bed')

    bed2gtf('tmp_subtract.bed', 'tmp_intergenic_region.gtf')


    intergenic_file = intergenic_dir + specie_file.replace('.gtf', '')+'__intergenic.gtf'
    c = open(intergenic_file, 'w')

    num =0 
    total =0
    for ln in open('tmp_intergenic_region.gtf'):
        int_length = int(ln.split('\t')[4])-int(ln.split('\t')[3])
        num +=1
        new_ln = ln.strip('#";\n') +'"'+str(num)+'";\n'
        c.write(new_ln)
        total+=int_length

    print 'intergenic region ', total, 'bp'

    c.close()
    os.system('rm tmp*')
