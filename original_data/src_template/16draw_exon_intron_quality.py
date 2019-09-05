# -*- coding=utf-8 -*-
from matplotlib import pyplot as plt
import numpy as np
import itertools
import random
import os

import scipy.stats as st
from scipy.misc import factorial
from math import exp
import pandas as pd

from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures

import seaborn as sns
str_img_output = '../results/img/'
os.system('mkdir '+str_img_output)
plt.switch_backend('agg')

from scipy.stats import pearsonr
def drawExonIntronDensityScatter(x, y, sampleId, outputdir):
    fig=plt.figure()
    ax = fig.add_subplot(111)
    plt.title("Correlationship between exon and intron in %s"%sampleId)
    plt.xlabel("Exon/RPKM")
    plt.ylabel("Intron/RPKM")
    plt.autoscale(tight=True)
    #plt.grid()##plt.show()
    xPercentile = np.percentile(x,[0,25,50,75,100])
    xIQR = xPercentile[3] - xPercentile[1]
    xUpLimit = xPercentile[3]+xIQR*1.5
    xDownLimit = xPercentile[1]-xIQR*1.5   
    yPercentile = np.percentile(y,[0,25,50,75,100])
    yIQR = yPercentile[3] - yPercentile[1]
    yUpLimit = yPercentile[3]+yIQR*1.5  
    yDownLimit = yPercentile[1]-yIQR*1.5 
    x_original = x
    y_original = y
    x_a = [a for a,b in zip(x,y)if a < xUpLimit]
    y_a = [b for a,b in zip(x,y)if a < xUpLimit]
    x = [a for a,b in zip(x_a,y_a)if b < yUpLimit]
    y = [b for a,b in zip(x_a,y_a)if b < yUpLimit]  
    x_b = [a for a,b in zip(x,y)if a > xDownLimit]
    y_b = [b for a,b in zip(x,y)if a > xDownLimit]  
    x = [a for a,b in zip(x_b,y_b)if b > yDownLimit]
    y = [b for a,b in zip(x_b,y_b)if b > yDownLimit]
    if x == [] or y == []:
        pearson_r = round(pearsonr(x_original, y_original)[0], 2)
    else:
        pearson_r = round(pearsonr(x, y)[0], 2)      
    
    '''
    fp2 = np.polyfit(x,y,3)
    f2 = np.poly1d(fp2)

    fx = np.linspace(0,x[-1],1000)
    plt.plot(fx,f2(fx),linewidth=2,color='g')  
    
    length = len(x)
    datasets_X = np.array(x).reshape([length,1])

    datasets_Y = np.array(y)

    minX = min(datasets_X)
    maxX = max(datasets_X)
    X = np.arange(minX,maxX).reshape([-1,1]) 
    print X
    plt.scatter(datasets_X, datasets_Y, color = 'red') 
    
    regr = linear_model.LinearRegression()
    regr.fit(datasets_X, datasets_Y)
    #print regr.__dict__['_residues'], regr.__dict__['coef_']
    plt.plot(X,regr.predict(X),color = 'blue',linewidth=4)    
    '''
    #ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
    str_one_text_content = 'diagonal line'
    str_pearson_text_content = 'pearson R is ' + str(pearson_r)
    one_font_dict = {'family':'serif',
                 'color':'black',
                 'size':14} 

    pearson_font_dict = {'family':'serif',
                 'color':'blue',
                 'size':14}                  
    plt.text(0.2, 0.9, str_one_text_content, ha='center', va='center', transform=ax.transAxes, fontdict=one_font_dict)
    plt.text(0.8, 0.2, str_pearson_text_content, ha='center', va='center', transform=ax.transAxes, fontdict=pearson_font_dict)
    
    if x == [] or y == []:
        #g = sns.jointplot(x1, x2, kind="kde", size=7, space=0) 
        #sns.regplot(x=x_original, y=y_original, fit_reg=False, color='red')   
        ax.scatter(x_original, y_original, color='red', alpha = 0.3)
        diag_line, = ax.plot(ax.get_ylim(), ax.get_ylim(), ls="-", c=".3", color='black')
        diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3", color='blue')
    else:
        #sns.regplot(x=x, y=y, fit_reg=False, color='red')
        ax.scatter(x, y, color='#BADCDB',alpha = 0.3)
        diag_line, = ax.plot(ax.get_ylim(), ax.get_ylim(), ls="-", c=".3", color='black')
        diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3", color='blue')
    plt.savefig(outputdir)


def read_csv_col_ndarray(csvDir, col, row_start):
    return np.loadtxt(csvDir,delimiter=',',unpack=True, skiprows=row_start, usecols=col)

def get_list_val_index_dic(list_a):
    dic_a = {}
    for i in range(len(list_a)):
        val = list_a[i]
        dic_a.setdefault(val, i)
    return dic_a 

# find out the replicate for each treatment    
def generate_correalation_matrix(gene_file, list_sample):
    df = pd.DataFrame()

    gene_sample_name_list = open(gene_file).readline().strip().split(',')
    gene_sample_val_index_dict = get_list_val_index_dic(gene_sample_name_list)

    for sample_name in list_sample:
        gene_count = read_csv_col_ndarray(gene_file, gene_sample_val_index_dict[sample_name], 1)
        df[sample_name] = gene_count
    return df.corr()

from matplotlib.pyplot import cm
def drawHeatmapCorrealation(count_file, list_samples, sample_id, output_dir):   
    fig=plt.figure()
    ax = fig.add_subplot(111)
    matrix_correalation = generate_correalation_matrix(count_file, list_samples)
    matrix_correalation = np.array(matrix_correalation)
    matrix_correalation = np.round(matrix_correalation, 2)

    im = ax.imshow(matrix_correalation, interpolation='nearest', cmap=plt.get_cmap('Pastel2_r'))
    ax.set_title("Correlation matrix for %s technical replication"%sample_id)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(list_samples)))
    ax.set_yticks(np.arange(len(list_samples)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(list_samples)
    ax.set_yticklabels(list_samples)
    colors = []
    for sampleId in list_samples:
        if sampleId == sample_id:
            colors.append('red')
        else:
            colors.append('black')
    for xtick, color in zip(ax.get_xticklabels(), colors):
        xtick.set_color(color)
    for ytick, color in zip(ax.get_yticklabels(), colors):
        ytick.set_color(color)        
    
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right",rotation_mode="anchor")
    for i in range(len(list_samples)):
        for j in range(len(list_samples)):
            text = ax.text(j, i, matrix_correalation[i, j], ha="center", va="center", color="b")
    # Loop over data dimensions and create text annotations.
    plt.savefig(output_dir)

def get_list_val_index_dic(list_a):
    dic_a = {}
    for i in range(len(list_a)):
        val = list_a[i]
        dic_a.setdefault(val, i)
    return dic_a 

def read_csv_col_ndarray(csvDir, col, row_start):
    return np.loadtxt(csvDir,delimiter=',',unpack=True, skiprows=row_start, usecols=col)

def get_density_list(csvDir, sample_col):
    length_array = read_csv_col_ndarray(csvDir, 2, 1)
    count_array = read_csv_col_ndarray(csvDir, sample_col, 1)
    return count_array/length_array

def get_intron_density_list(csvGeneDir, csvExonDir, gene_col, exon_col):
    gene_length_array = read_csv_col_ndarray(csvGeneDir, 2, 1)
    exon_length_array = read_csv_col_ndarray(csvExonDir, 2, 1)
    gene_count_array = read_csv_col_ndarray(csvGeneDir, gene_col, 1)
    exon_count_array = read_csv_col_ndarray(csvExonDir, exon_col, 1)
    length_array = gene_length_array - exon_length_array
    count_array = gene_count_array - exon_count_array
    return count_array/length_array


def main():
    gene_featureCount = '../results/gene_feature_count_matrix.csv'
    exon_featureCount = '../results/exon_feature_count_matrix.csv'
    gene_sample_name_list = open(gene_featureCount).readline().strip().split(',')
    gene_sample_val_index_dict = get_list_val_index_dic(gene_sample_name_list)
    exon_sample_name_list = open(exon_featureCount).readline().strip().split(',')
    exon_sample_val_index_dict = get_list_val_index_dic(exon_sample_name_list)          
    
    os.system('mkdir ../results/img/')
    scatter_dir = '../results/img/exon_intron_scatter/'
    os.system('mkdir %s'%scatter_dir)
    
    for sampleId in gene_sample_name_list:
        if sampleId in ['ensembl_id', 'length', 'hugo_id', 'biotype']:
            continue
        exon_density = get_density_list(exon_featureCount, exon_sample_val_index_dict[sampleId])
        gene_density = get_density_list(gene_featureCount, gene_sample_val_index_dict[sampleId])
        #intron_density = get_intron_density_list(gene_featureCount, exon_featureCount, gene_sample_val_index_dict[sampleId], exon_sample_val_index_dict[sampleId])
        exonIntronScatterDir = "%sscatter_exon_intron_%s.png"%(scatter_dir, sampleId)
        drawExonIntronDensityScatter(exon_density*1000, gene_density*1000, sampleId, exonIntronScatterDir)
    

    gsm_detail_file = '/data1/szx/grobase-data-analysis/original_data/gsm_detail_results.csv'
    heatmap_dir = '../results/img/repeat_heatmap/'
    os.system('mkdir %s'%heatmap_dir)    
    dic_process_gsm = {}
    for ln in open(gsm_detail_file):
        ls = ln.split('!')
        gsm, gse, specie,tissue, treatment, = ls[0].strip('"'), ls[3], ls[4], ls[6], ls[9]
        process_key = gse+specie+tissue+treatment
        if process_key in dic_process_gsm.keys():
            dic_process_gsm[process_key].append(gsm)
        else:
            dic_process_gsm.setdefault(process_key, [gsm])

    exist_gsm = set([x.split('.')[0] for x in os.listdir('../results/merge_bam')])
    for list_samples in dic_process_gsm.values():
        if len(list_samples) > 1:
            for sample_id in list_samples:
                if sample_id in exist_gsm:
                    repeatCorScatterDir = '%srepeat_%s.png'%( heatmap_dir, sample_id)
                    drawHeatmapCorrealation(gene_featureCount, list_samples, sample_id, repeatCorScatterDir)  

if __name__ == '__main__':
    main()    