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
commit[地址]()
1将fq文件移动到各种物种文件夹下。  
2将faidx文件移动到各种物种文件夹下。  

## 物种层级  
