#包的安装
local({r <- getOption("repos")
r["CRAN"] <- "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"
options(repos=r)})

if (!requireNamespace("BiocManager", quietly=TRUE)){
  install.packages("BiocManager")
}
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

#查看project中有哪些数据类型
TCGAbiolinks:::getProjectSummary("TCGA-SKCM")


#查询数据  GDCquery()
query <- GDCquery(project = "TCGA-SKCM",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification",
                  workflow.type = "BCGSC miRNA Profiling")

## Not run:
query1 <- GDCquery(project = "TCGA-SKCM",
                   data.category = "Clinical")

query2 <- GDCquery(project = "TCGA-SKCM",
                  data.category = "Simple Nucleotide Variation",
                  data.type = "Masked Somatic Mutation")
View(getResults(query2))

query <- GDCquery(project = "TCGA-ACC",
                  data.category =  "Copy Number Variation",
                  data.type = "Masked Copy Number Segment",
                  sample.type = c("Primary solid Tumor"))
query.met <- GDCquery(project = c("TCGA-GBM","TCGA-LGG"),
                      legacy = TRUE,
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450")
query <- GDCquery(project = "TCGA-ACC",
                  data.category =  "Copy number variation",
                  legacy = TRUE,
                  file.type = "hg19.seg",
                  barcode = c("TCGA-OR-A5LR-01A-11D-A29H-01"))
query <-GDCquery(project = "TCGA-GBM",
                 data.category = "Gene expression",
                 data.type = "Gene expression quantification",
                 platform = "Illumina HiSeq", 
                 file.type  = "normalized_results",
                 experimental.strategy = "RNA-Seq",
                 barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
                 legacy = TRUE)



#下载数据  GDCdownload()
GDCdownload(query, method = "client", files.per.chunk = 10, directory="E:/download")
GDCdownload(query1, method = "api",directory="E:/download")
GDCdownload(query2, method = "api",directory="E:/download", files.per.chunk = 1)
GDCdownload(query2, method = "client", files.per.chunk = 10, directory="E:/download")

#保存整理数据 GDCprepare()
data <- GDCprepare(query = query,
                   save = TRUE,
                   directory ="E:/download",   #注意和GDCdownload设置的路径一致GDCprepare才可以找到下载的数据然后去处理。    
                   save.filename = "GBM.RData")  #存储一下，方便下载直接读取


data1 <- GDCprepare(query = query1,
                   save = TRUE,
                   directory ="E:/download", 
                   save.filename = "Clinical.RData")   
data2 <- GDCprepare(query = query2,
                   save = TRUE,
                   directory ="E:/download",  
                   save.filename = "Mutation.RData")

#基因ID转换
miRNA<-data.frame(data$miRNA_ID)
left_alias1<-merge(miRNA,change_2_2,by.x = "data.miRNA_ID",by.y = "alias_symbol",all.x=TRUE)
write.csv(left_alias3,file="D:/3melanoma/left_alias3.csv")


setwd("C:/Users/Desktop/circRNA分析结果/find_circ_result自己分析结果")#设置working directory
gene_symbol<-read.csv("SFTSV_24vscontrol_circBase_anno.csv",header=F,stringsAsFactors = F)[,11]#读取数据并提取含有gene_symbol的列
library(biomaRt)
mart <- useMart("ensembl","hsapiens_gene_ensembl")##小鼠选择mmusculus_gene_ensembl
gene_id<-getBM(attributes=c("external_gene_name","ensembl_gene_id"),filters = "external_gene_name",values = gene_symbol, mart = mart)#将输入的filters设置未external_gene_name(也就是gene_symbol),将输出的attributes设置为external_gene_name和emsembl_gene_id
write.table(gene_id,"SFTSV_24vscontrol_circBase_anno_gene_id.txt",row.names = F,col.names=F,quote=F)

#基因筛选
data$sum = apply(data[,substr(colnames(data),1,9) %in% "reads_per"],1,sum)


#MAD的计算：选择均小于0.5的代谢基因
median(abs(d-median(d)))*1.4826