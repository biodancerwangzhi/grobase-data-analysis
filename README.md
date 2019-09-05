# grobase-data-analysis
数据批量分析，分两个层级，整个项目层级，单个物种层级。  
项目层级负责0构建脚手架

## 总层次
### 0构建项目脚手架(linux环境下)  
commit[地址](https://github.com/biodancerwangzhi/grobase-data-analysis/commit/a249c547fed542064f19c47e6e66fd28d41a844f/src/0project-cli.py)  
原始数据是一张具有物种对应的gsm id和srr id信息的总表。  
先创建一个bulk_specie_results文件夹，里面包含各个物种的文件夹。  
为每个物种文件夹建src,results，logs，original_data四个文件夹。  
需要拷贝每个物种需要使用的脚本模板进src，脚本根据物种名不同，使用不同的注释文件，动态修改一下模板。 
original_data有生成的每个物种的gsm id和srr id，用于后面的分析。之后就是在物种文件夹下批量运行脚本。  

### 0-1基于总表检验一些关键文件是否已经存在
commit[地址](https://github.com/biodancerwangzhi/grobase-data-analysis/commit/b77bea5b4b86b04621d5f4df929e41426808013d/src/0-1test_preprocess.py)
1检测注释文件。  
需要的注释文件，基因组文件，索引。一个不能少。
2检测原始sra文件。
检测新下载的sra是不是之前的分析结果中使用，如果使用文件大小是不是相同。  
3检验gsm id。
检验要分析的gsm id是否已经在产生的bam中。  

### 0-2移动关键文件  
commit[地址](https://github.com/biodancerwangzhi/grobase-data-analysis/commit/a49511e5cabab023547b88be301f98d0e41a6717/src/0-2mv_files.py)
1批量将sra变成fq，然后将fq文件移动到各种物种文件夹下。  
2将faidx文件移动到各种物种文件夹下。  

### 1-0统计所有bam文件的fastq 结果  
commit[地址]()
整理到一张excel表。   

### 1-1统计各种比对率  
commit[地址]()  
整理到一张execel表。  

## 物种层级  
### 0-0sra-bam的基础运算
commit[地址](https://github.com/biodancerwangzhi/grobase-data-analysis/commit/a249c547fed542064f19c47e6e66fd28d41a844f#diff-b6930e80c1c53fc9eb467f7a76d4132c/original_data/src_template/3trim_polya.py)
fq->bam  
包括常见的reads文件trim和碱基质控。比对，排序，建索引等。模板文件[3-8]  
bam stat和fastqc  
简单的比对统计和质控。模板文件[9-10]  

### 0-1各种feature count  
commit[地址](https://github.com/biodancerwangzhi/grobase-data-analysis/commit/c066a53c5a17c145e2e9b53066bc424da2440878/original_data/src_template/12promoter_geneBody_featureCount.py)  
需要各种注释文件生成，见others的注释生成部分。
1 gene和exon层面的featureCount。  
2 promoter和gene body的featureCount。  
3 promoter和gene body windows 的featureCount。  
4 intergenic region的featureCount。  

### 0-2提取express matrix  
commit[地址](https://github.com/biodancerwangzhi/grobase-data-analysis/commit/c066a53c5a17c145e2e9b53066bc424da2440878/original_data/src_template/14export_express_matrix.py )
上面生成的count需要输出到excel。  

### 0-3计算暂停基因  
commit[地址](https://github.com/biodancerwangzhi/grobase-data-analysis/commit/cefeead6b14b4aa5e09eea54b9eb772e586c76c8)
1 确定是否激活: 
    通过一些基础函数提取样本在express matrix的具体列。
    使用possion分布比较gene body的density和intergenic region（背景）的density。  

    这里还尝试了和RNA-Seq比一比激活基因的个数，
        从archs4下了所有rna-seq的表达矩阵，
        获取了具体gsm的count list，
        找到了网上的两种方法1 FPKM>1 2 zFPKM > -3，  
        但其中有很多问题，一个细胞系有很多样本选哪些比较，gro和rna-seq的gene索引找到表达值比较麻烦（懒得写），gro数据还要重算，而且比来比去也没多大意思对数据库没多大影响  

2 PI(folder change)模型
3 fisher test 模型 
    按公式计算
    整合到一个excel中

### 0-4各种数据的统计量分布情况。  
所有分析结果都要生成基于样本的统计。  
将分布的情况转化成百分比。  
设计图表，基于原始数据或百分比。  

### 0-5画质控图  
commit[地址]()

## others
### 0-0 annotation生成  
commit[地址](https://github.com/biodancerwangzhi/grobase-data-analysis/commit/c066a53c5a17c145e2e9b53066bc424da2440878/src/others/generate_annotation/4pause_gene_windows.py)
1 确定暂停基因需要promoter区和gene body区的reads密度，所以要构建gene的promoter区和gene body区。  
2 暂停模型比较复杂，还需要建一个滑窗去做统计分析。也会考虑背景值，intergenic region的density。  
