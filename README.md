# grobase-data-analysis
数据批量分析，分两个层级，整个项目层级，单个物种层级。  
项目层级负责0构建脚手架
## 0构建项目脚手架(linux环境下)  
脚本[地址]()  
原始数据是一张具有物种对应的gsm id和srr id信息的总表。  
先创建一个bulk_specie_results文件夹，里面包含各个物种的文件夹。  
为每个物种文件夹建src,results，logs，original_data四个文件夹。  
需要拷贝每个物种需要使用的脚本模板进src，脚本根据物种名不同，使用不同的注释文件，动态修改一下模板。 
original_data有每个物种的gsm id和srr id，用于后面的分析。之后就是在物种文件夹下批量运行脚本。  
