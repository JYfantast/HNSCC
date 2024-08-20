
c("#2EC4B6","#E71D36")
c("#E71D36","#2EC4B6")

c("#0048A1","#AA1701")
c("#AA1701","#0048A1")

####1.1.TCGA-mRNA####

#加载R包
library(rjson)
library(tidyverse)

#设置路径 
setwd("E:/99other/31HNSCC/1data/TCGA/RNA")

#读入meta.data文件
json <- jsonlite::fromJSON("metadata.cart.2024-03-23.json")
#View(json)

#获取样本名称及文件名称
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
#sample_id[1:10]
file_sample <- data.frame(sample_id,file_name=json$file_name)  
#View(file_sample)

#获取counts所在位置
count_file <- list.files('gdc_download_20240323_130603.971589/',
                         pattern = '*.tsv',recursive = TRUE)
#count_file[1:10]

#获取每个文件名称
count_file_name <- strsplit(count_file,split='/')
#count_file_name[1:10]
count_file_name <- sapply(count_file_name,function(x){x[2]})
#count_file_name[1:10]

#构建一个空的数据框
matrix = data.frame(matrix(nrow=60660,ncol=0))

#逐个读取及合并
for (i in 1:length(count_file)){
  path = paste0('gdc_download_20240323_130603.971589//',count_file[i])   #Counts文件夹名
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  #3：unstranded,counts；4：stranded_first；5：stranded_second；6：tpm_unstranded；7：fpkm_unstranded；8：fpkm_uq_unstranded
  #data <- data[3]
  data <- data[6]
  #data <- data[7]
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}

#转化为gene_symbol
path = paste0('gdc_download_20240323_130603.971589//',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
#gene_name[1:10]
matrix0 <- cbind(gene_name,matrix)
#获取基因类型
gene_type <- data[-c(1:6),2]
#gene_type[1:10]
matrix0 <- cbind(gene_type,matrix0)

#将gene_name列去除重复的基因，保留基因最大表达量结果,min,mean
matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)
table(gene_name)

#保留mRNA
matrix0 <- subset(x = matrix0, gene_type == "protein_coding")
#table(gene_type)

#将gene_name列设为行名
rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-c(1,2)]

samenames = rownames(matrix0)

#转化为数值型
matrix0 = apply(matrix0,2,as.numeric)
matrix0 = as.data.frame(matrix0)
rownames(matrix0) = samenames

#创建mo.data
mo.data<-list()
mo.data$mRNA.TPM = matrix0

#保存
saveRDS(mo.data, file = "E:/99other/31HNSCC/1data/TCGA/results/mo.data1mRNA.rds")

dir.create('results')
mo.data$mRNA.TPM %>% {log2(1+.)} %>% write.table(gzfile("E:/99other/31HNSCC/1data/TCGA/results/tcga.mRNA.exp.all.txt.gz"),sep="\t",quote = F)


####1.2.TCGA-lncRNA####

#加载R包
library(rjson)
library(tidyverse)

#设置路径 
setwd("E:/99other/31HNSCC/1data/TCGA/RNA")

#读入meta.data文件
json <- jsonlite::fromJSON("metadata.cart.2024-03-23.json")
#View(json)

#获取样本名称及文件名称
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
#sample_id[1:10]
file_sample <- data.frame(sample_id,file_name=json$file_name)  
#View(file_sample)

#获取counts所在位置
count_file <- list.files('gdc_download_20240323_130603.971589/',
                         pattern = '*.tsv',recursive = TRUE)
#count_file[1:10]

#获取每个文件名称
count_file_name <- strsplit(count_file,split='/')
#count_file_name[1:10]
count_file_name <- sapply(count_file_name,function(x){x[2]})
#count_file_name[1:10]

#构建一个空的数据框
matrix = data.frame(matrix(nrow=60660,ncol=0))

#逐个读取及合并
for (i in 1:length(count_file)){
  path = paste0('gdc_download_20240323_130603.971589//',count_file[i])   #Counts文件夹名
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  #3：unstranded,counts；4：stranded_first；5：stranded_second；6：tpm_unstranded；7：fpkm_unstranded；8：fpkm_uq_unstranded
  #data <- data[3]
  data <- data[6]
  #data <- data[7]
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}

#转化为gene_symbol
path = paste0('gdc_download_20240323_130603.971589//',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
#gene_name[1:10]
matrix0 <- cbind(gene_name,matrix)
#获取基因类型
gene_type <- data[-c(1:6),2]
#gene_type[1:10]
matrix0 <- cbind(gene_type,matrix0)

#将gene_name列去除重复的基因，保留基因最大表达量结果,min,mean
matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)
#table(gene_name)

#保留mRNA
matrix0 <- subset(x = matrix0, gene_type == "lncRNA")
#table(gene_type)

#将gene_name列设为行名
rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-c(1,2)]
samenames = rownames(matrix0)

#转化为数值型
matrix0 = apply(matrix0,2,as.numeric)
matrix0 = as.data.frame(matrix0)
rownames(matrix0) = samenames
matrix0 = log2(matrix0 + 1)

#读取
mo.data <- readRDS("E:/99other/31HNSCC/1data/TCGA/results/mo.data1mRNA.rds")

mo.data$lncRNA.exp = matrix0

#保存
saveRDS(mo.data, file = "E:/99other/31HNSCC/1data/TCGA/results/mo.data2lncRNA.rds")


####1.3.TCGA-miRNA####

#加载R包
library(stringr)
library(tidyverse)

#设置路径
setwd("E:/99other/31HNSCC/1data/TCGA/miRNA")

#批量读入
count_files = dir("gdc_download_20240323_130711.898685/",pattern = "*isoforms.quantification.txt$",recursive = T)
exp = list()
for(i in 1:length(count_files)){
  exp[[i]] <- read.table(paste0("gdc_download_20240323_130711.898685/",count_files[[i]]),sep = "\t",header = T) %>%
    dplyr::select(c(6,4)) %>% ###count 3, RPM 4
    group_by(miRNA_region) %>% 
    summarise(reads_per_million_miRNA_mapped = sum(reads_per_million_miRNA_mapped)) 
}

#行数是不一样，合并
m = Reduce(function(x, y) merge(x, y, by= 'miRNA_region',all = T), exp)
#将NA值变为0
m[is.na(m)]=0
#行名转换
exp <- column_to_rownames(m,var = "miRNA_region")
#删除最后三行
exp = exp[-((nrow(exp)-2):nrow(exp)),]
#查看
#table(str_detect(rownames(exp),"mature,"))
#删除行名中的mature
rownames(exp) = str_remove(rownames(exp),"mature,")

#MIM开头的是mirBase数据库里的id，需要转换为大多以hsa-miR开头的成熟体id
#GDC数据库使用的mirBasev21版本的id，使用miRBaseVersions.db包更丝滑
library(miRBaseVersions.db)

#基因名字修改
mh <- select(miRBaseVersions.db,
             keys = rownames(exp),
             keytype = "MIMAT",
             columns = c("ACCESSION","NAME","VERSION"))
mh = mh[mh$VERSION=="21",]
mh = mh[match(rownames(exp),mh$ACCESSION),]
#identical(rownames(exp),mh$ACCESSION)
#table(!duplicated(mh$NAME))
rownames(exp) = mh$NAME

#读入
meta <- jsonlite::fromJSON("metadata.cart.2024-03-23.json")
ID = sapply(meta$associated_entities,
            function(x){x$entity_submitter_id})

#注意修改样本名字
file2id = data.frame(file_name = meta$file_id,
                     ID = ID)
count_files2 = stringr::str_split(count_files,"/",simplify = T)[,1]
file2id = file2id[match(count_files2,file2id$file_name),]
colnames(exp) = file2id$ID

exp = as.matrix(exp)

exp = log2(exp + 1)
exp = as.data.frame(exp)
#读取
mo.data <- readRDS("E:/99other/31HNSCC/1data/TCGA/results/mo.data2lncRNA.rds")

mo.data$miRMA.exp = exp

#保存
saveRDS(mo.data, file = "E:/99other/31HNSCC/1data/TCGA/results/mo.data3miRNA.rds")


####1.4.TCGA-methy####

#加载R包
library(rjson)
library(tidyverse)

#设置路径
setwd("E:/99other/31HNSCC/1data/TCGA/methy")

#读入meta.data文件
json <- jsonlite::fromJSON("metadata.cart.2024-03-23.json")

#获取样本名称及文件名称
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)  

length(rownames(file_sample))

file_sample1 = file_sample[c(1:300),]
file_sample2 = file_sample[c(301:580),]

write.table(file_sample1,'samples1.txt', sep="\t", quote=F, row.names = F)
write.table(file_sample2,'samples2.txt', sep="\t", quote=F, row.names = F)




#加载R包
library(stringr)
library(tidyverse)

#设置路径
setwd("E:/99other/31HNSCC/1data/TCGA/methy")

  #读入meta.data文件
json <- jsonlite::fromJSON("metadata.cart.2024-03-23.json")

#获取样本名称及文件名称
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)  

#获取矩阵所在位置
count_file <- list.files('1/',
                         pattern = '*level3betas.txt',recursive = TRUE)

#获取每个文件名称
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})

#构建一个空的数据框
matrix = data.frame(matrix(nrow=486427,ncol=0))

#逐个读取及合并
for (i in 1:length(count_file)){
  data<- read.delim(paste0('1//',count_file[i]),fill = TRUE,header = FALSE,row.names = 1)
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}

#删除任意带NA值的行
matrix = na.omit(matrix)
matrix = as.data.frame(matrix)

#读取
mo.data <- readRDS("E:/99other/31HNSCC/1data/TCGA/results/mo.data3miRNA.rds")

mo.data$meth.beta = matrix

#保存
saveRDS(mo.data, file = "E:/99other/31HNSCC/1data/TCGA/results/mo.data4methy.rds")

####1.5.TCGA-mutation####

#加载R包
library(maftools)
library(dplyr)

#目录
setwd("E:/99other/31HNSCC/1data/TCGA/mutation/gdc_download_20240323_131624.874215")
#合并所有数据
files <- list.files(pattern = '*.gz',recursive = TRUE)
all_mut <- data.frame()
for (file in files) {
  mut <- read.delim(file,skip = 7, header = T, fill = TRUE,sep = "\t")
  all_mut <- rbind(all_mut,mut)
}
#仅保留前12个字符
#all_mut$Tumor_Sample_Barcode = substr(all_mut$Tumor_Sample_Barcode,1,12)

#数据读入
all_mut <- read.maf(all_mut)

#数据整理
a <- all_mut@data %>%
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>%
  as.data.frame() #%>%
#  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))

#提取基因
gene <- as.character(unique(a$Hugo_Symbol))
#提取样本
sample <- as.character(unique(a$Tumor_Sample_Barcode))

#创建data.frame
mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
#将信息填入
for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}

#创建data.frame
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))
#将信息填入
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

#所有样本突变情况汇总/排序
gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))

colnames(gene_count)[1] = "Gene"
colnames(gene_count)[2] = "Num"
write.table(gene_count,'geneMut.txt', sep="\t", quote=F, row.names = F)
write.mafSummary(maf = all_mut,basename = "input")

#绘制瀑布图oncoplot
pdf(file="maf.pdf", width=6, height=6)
oncoplot(maf = all_mut,
         top = 30, #显示前30个的突变基因信息
         fontSize = 0.6, #设置字体大小
         showTumorSampleBarcodes = F) #不显示病人信息
dev.off()

#计算tmb值
tmb_table = tmb(maf = all_mut,logScale = F)
#保留需要tmb值信息
tmb_table = tmb_table[,c(1,3)]
tmb_table = as.data.frame(tmb_table)
tmb_table[,1]=substr(tmb_table[,1],1,12)
colnames(tmb_table)
tmb_table <- aggregate( . ~ Tumor_Sample_Barcode,data=tmb_table, max)
colnames(tmb_table)[1] = "id"
colnames(tmb_table)[2] = "TMB"
write.table(tmb_table,'TMB.txt', sep="\t", quote=F, row.names = F)

#读取
mo.data <- readRDS("E:/99other/31HNSCC/1data/TCGA/results/mo.data4methy.rds")

mo.data$mut.status = mat_0_1

#保存
saveRDS(mo.data, file = "E:/99other/31HNSCC/1data/TCGA/results/mo.data5mutation.rds")

####1.6.TCGA-clinical####

library(tidyverse)
library(data.table)

#读取生存数据文件
tcga.cli=read.table("E:/99other/31HNSCC/1data/TCGA/clinical/clinical.txt", header=T, sep="\t", check.names=F, row.names=1)

#读入
mo.data <- readRDS("E:/99other/31HNSCC/1data/TCGA/results/mo.data5mutation.rds")

#保留共有的样本
#保留第1个第14-16位为"01A"的样本
fix_sample_names <- function(x){
  x <- x[,substr(colnames(x),14,16)=='01A']
  colnames(x) <- substr(colnames(x),1,12)
  x <- x[,!duplicated(colnames(x))]
  return(x)
}

mo.data <- mo.data %>% lapply(FUN = fix_sample_names )
select_sample <- append(lapply(mo.data,colnames),list(tcga.cli=rownames(tcga.cli))) %>% Reduce(intersect,.)
mo.data <- mo.data %>% lapply(function(x) x[,select_sample]) %>% lapply(function(x) x[rowSums(x)>0,])
tcga.cli <- tcga.cli[select_sample,]

#保存
write.table(data.frame(ID=rownames(tcga.cli),tcga.cli),file="E:/99other/31HNSCC/1data/TCGA/results/tcga.cli.txt", sep="\t", quote=F, row.names = F,col.names = T)
#tcga.cli %>% write.table("E:/99other/31HNSCC/1data/TCGA/results/tcga.cli.txt",sep="\t",quote = F)
mo.data$mRNA.TPM %>% write.table(gzfile("E:/99other/31HNSCC/1data/TCGA/results/tcga.mRNA.TPM.txt.gz"),sep="\t",quote = F)                                       
mo.data$mRNA.exp <- mo.data$mRNA.TPM %>% {log2(1+.)}
mo.data$mRNA.exp %>% write.table(gzfile("E:/99other/31HNSCC/1data/TCGA/results/tcga.mRNA.exp.txt.gz"),sep="\t",quote = F)
mo.data$mRNA.TPM <- NULL                                                        
saveRDS(mo.data,'E:/99other/31HNSCC/1data/TCGA/results/mo.data.rds',compress = F)



####1.7.GSE27020-GPL96####

setwd("E:/99other/31HNSCC/1data/GEO/GSE27020-GPL96")
#下载矩阵
gset <- getGEO("GSE27020",destdir = "E:/99other/31HNSCC/1data/GEO/GSE27020-GPL96",AnnotGPL = F,getGPL = F) 
gset
#提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a=gset[[1]]
#提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
dat=exprs(a)
#展示表达矩阵的前6行，看下数据是否经过log转换，一般数据在20以内，是经过log转换的，若有成百上千的数据，表示这个数据还需要进一步log处理。
head(dat)
#也可以输入上面这段代码，如果数据已经经log后，显示log2 transform not needed；如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NA
dat <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
#查看a的临床信息，为后面选择用于分组的变量做准备
pd=pData(a)
#查看pd
#View(pd)
#输出临床信息
write.csv(pd,'clinicalGSE27020.csv',row.names = TRUE)

#另外一种获取GPL文件方法
#如果没有下载过这个包，就下载一下，代码如下
#BiocManager::install("hgu133plus2.db")
#提取信息
library(hgu133a.db)
ls("package:hgu133a.db") #大致查看R包中的信息，寻找我们需要的SYMBOL
ids <- toTable(hgu133aSYMBOL)  #用toTable（）函数提取
head(ids) #查看提取内容
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")
#exp <- mutate(ex,probe_id=rownames(probe_exp))#将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
#exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id") # 合并数据
#exp2=exp2[!duplicated(probe_exp$symbol),]# 按照symbol列去重
#rownames(exp2)=exp2$symbol # 数据框probe_exp的行名变成symbol
#gene_exp_matrix <- na.omit(exp2 ) # 去空值
#输出文件
#write.table(gene_exp_matrix2,file = "Gastric.cancer.geneid.exprs.txt",sep = "\t",row.names=T,col.names = T)

#另外一种方式
#下载注释文件
#gpl <- getGEO('GPL10558', destdir="")
#我们将平台文件转为list格式，并赋值给gpl1，将gpl1保存为R文件，方便以后调用。
#gpl1<-Table(gpl)
#查看平台文件的列名，我们看到有ID和gene symbol,记住gene symbol在第几列，我们这里在第11列。
#colnames(Table(gpl))
#再了解一下平台文件的数据，当然这里大家可以直接选择gene symbol字段也行，主要了解一下symbol中基因symbol值是否只有一个。
#View(gpl1)
#我们这里的gene symbol字段中的symbol，有的基因就不止一个名称，///后面有重名，我们需要第一个名字，所以需要用字符分割包再处理一下,原理同上面处理title一样。
#提取平台文件gpl1中的ID和gene symbol，并赋值给probe2symbol_df
#gpl1$`Gene Symbol`=str_split(gpl1[,11],'///',simplify = T)[,1]
#probe2symbol_df<-gpl1[,c(1,11)]

probe2symbol_df <- ids
#将列名改为probe_id和symbol
#这两步是因为我懒，不想再调整代码，为了方便后面代码运行，我改的，大家也可以不改。
colnames(probe2symbol_df)=c('probe_id','symbol')
#查看symbol为唯一值的个数
length(unique(probe2symbol_df$symbol))
#查看每个基因出现n次的个数,我们可以看到，symbol出现一次的基因有7050个，出现2次有1493个。。。
table(sort(table(probe2symbol_df$symbol)))
#去掉没有注释symbol的探针（其实这里没有注释的探针数量即为上面出现次数最多的基因440，也就是说有440个探针没有symbol）
ids=probe2symbol_df[probe2symbol_df$symbol != '',]
#%in%用于判断是否匹配，
#注意这里dat是gset表达矩阵的数据，这一步就是把平台文件的ID和矩阵中的ID匹配。
ids=probe2symbol_df[probe2symbol_df$probe_id %in%rownames(dat),]
#取表达矩阵中可以与探针名匹配的那些，去掉无法匹配的表达数据
dat=dat[ids$probe_id,]
#ids新建mean这一列，列名为mean，同时对dat这个矩阵按行操作，取每一行（即每个样本在这个探针上的）的均值，将结果给到mean这一列的每一行，这里也可以改为中位值，median.
ids$mean=apply(dat,1,mean)
#即先按symbol排序，相同的symbol再按照均值从大到小排列，方便后续保留第一个值。
ids=ids[order(ids$symbol,ids$mean,decreasing = T),]
#将symbol这一列取出重复项，'!'为否，即取出不重复的项，去除重复的gene
#取median的话，这样就保留每个基因最大表达量结果.最后得到n个基因。
ids=ids[!duplicated(ids$symbol),]
#新的ids取出probe_id这一列，将dat按照取出的这一列，用他们的每一行组成一个新的dat
dat=dat[ids$probe_id,]
#把ids的symbol这一列中的每一行给dat作为dat的行名
rownames(dat)=ids$symbol
View(dat)
#但是我们看到矩阵最后一行，是没有symbol名字的，我们把他去掉，数字自己更改。
dim(dat)
#dat<-dat[-13437,]
#rownames(dat)[20911] = "A1BG"
#最后我们把结果保存。
#write.csv(dat,file="2GEO/GSE68465/GSE68465.csv")
#write.table(dat,file="GEO.txt", sep="\t", quote=F, row.names = TRUE)
write.table(data.frame(ID=rownames(dat),dat),file="GSE27020.txt", sep="\t", quote=F, row.names = F,col.names = T)



####1.8.GSE41613-GPL570####

setwd("E:/99other/31HNSCC/1data/GEO/GSE41613-GPL570")
#下载矩阵
gset <- getGEO("GSE41613",destdir = "E:/99other/31HNSCC/1data/GEO/GSE41613-GPL570",AnnotGPL = F,getGPL = F) 
gset
#提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a=gset[[1]]
#提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
dat=exprs(a)
#展示表达矩阵的前6行，看下数据是否经过log转换，一般数据在20以内，是经过log转换的，若有成百上千的数据，表示这个数据还需要进一步log处理。
head(dat)
#也可以输入上面这段代码，如果数据已经经log后，显示log2 transform not needed；如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NA
dat <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
#查看a的临床信息，为后面选择用于分组的变量做准备
pd=pData(a)
#查看pd
#View(pd)
#输出临床信息
write.csv(pd,'clinicalGSE41613.csv',row.names = TRUE)

#另外一种获取GPL文件方法
#如果没有下载过这个包，就下载一下，代码如下
#BiocManager::install("hgu133plus2.db")
#提取信息
library(hgu133plus2.db)
ls("package:hgu133plus2.db") #大致查看R包中的信息，寻找我们需要的SYMBOL
ids <- toTable(hgu133plus2SYMBOL)  #用toTable（）函数提取
head(ids) #查看提取内容
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")
#exp <- mutate(ex,probe_id=rownames(probe_exp))#将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
#exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id") # 合并数据
#exp2=exp2[!duplicated(probe_exp$symbol),]# 按照symbol列去重
#rownames(exp2)=exp2$symbol # 数据框probe_exp的行名变成symbol
#gene_exp_matrix <- na.omit(exp2 ) # 去空值
#输出文件
#write.table(gene_exp_matrix2,file = "Gastric.cancer.geneid.exprs.txt",sep = "\t",row.names=T,col.names = T)

#另外一种方式
#下载注释文件
#gpl <- getGEO('GPL10558', destdir="")
#我们将平台文件转为list格式，并赋值给gpl1，将gpl1保存为R文件，方便以后调用。
#gpl1<-Table(gpl)
#查看平台文件的列名，我们看到有ID和gene symbol,记住gene symbol在第几列，我们这里在第11列。
#colnames(Table(gpl))
#再了解一下平台文件的数据，当然这里大家可以直接选择gene symbol字段也行，主要了解一下symbol中基因symbol值是否只有一个。
#View(gpl1)
#我们这里的gene symbol字段中的symbol，有的基因就不止一个名称，///后面有重名，我们需要第一个名字，所以需要用字符分割包再处理一下,原理同上面处理title一样。
#提取平台文件gpl1中的ID和gene symbol，并赋值给probe2symbol_df
#gpl1$`Gene Symbol`=str_split(gpl1[,11],'///',simplify = T)[,1]
#probe2symbol_df<-gpl1[,c(1,11)]

probe2symbol_df <- ids
#将列名改为probe_id和symbol
#这两步是因为我懒，不想再调整代码，为了方便后面代码运行，我改的，大家也可以不改。
colnames(probe2symbol_df)=c('probe_id','symbol')
#查看symbol为唯一值的个数
length(unique(probe2symbol_df$symbol))
#查看每个基因出现n次的个数,我们可以看到，symbol出现一次的基因有7050个，出现2次有1493个。。。
table(sort(table(probe2symbol_df$symbol)))
#去掉没有注释symbol的探针（其实这里没有注释的探针数量即为上面出现次数最多的基因440，也就是说有440个探针没有symbol）
ids=probe2symbol_df[probe2symbol_df$symbol != '',]
#%in%用于判断是否匹配，
#注意这里dat是gset表达矩阵的数据，这一步就是把平台文件的ID和矩阵中的ID匹配。
ids=probe2symbol_df[probe2symbol_df$probe_id %in%rownames(dat),]
#取表达矩阵中可以与探针名匹配的那些，去掉无法匹配的表达数据
dat=dat[ids$probe_id,]
#ids新建mean这一列，列名为mean，同时对dat这个矩阵按行操作，取每一行（即每个样本在这个探针上的）的均值，将结果给到mean这一列的每一行，这里也可以改为中位值，median.
ids$mean=apply(dat,1,mean)
#即先按symbol排序，相同的symbol再按照均值从大到小排列，方便后续保留第一个值。
ids=ids[order(ids$symbol,ids$mean,decreasing = T),]
#将symbol这一列取出重复项，'!'为否，即取出不重复的项，去除重复的gene
#取median的话，这样就保留每个基因最大表达量结果.最后得到n个基因。
ids=ids[!duplicated(ids$symbol),]
#新的ids取出probe_id这一列，将dat按照取出的这一列，用他们的每一行组成一个新的dat
dat=dat[ids$probe_id,]
#把ids的symbol这一列中的每一行给dat作为dat的行名
rownames(dat)=ids$symbol
View(dat)
#但是我们看到矩阵最后一行，是没有symbol名字的，我们把他去掉，数字自己更改。
dim(dat)
#dat<-dat[-13437,]
#rownames(dat)[20911] = "A1BG"
#最后我们把结果保存。
#write.csv(dat,file="2GEO/GSE68465/GSE68465.csv")
#write.table(dat,file="GEO.txt", sep="\t", quote=F, row.names = TRUE)
write.table(data.frame(ID=rownames(dat),dat),file="GSE41613.txt", sep="\t", quote=F, row.names = F,col.names = T)


####1.9.GSE65858-GPL10558####

setwd("E:/99other/31HNSCC/1data/GEO/GSE65858-GPL10558")
#下载矩阵
gset <- getGEO("GSE65858",destdir = "E:/99other/31HNSCC/1data/GEO/GSE65858-GPL10558",AnnotGPL = F,getGPL = F) 
gset
#提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a=gset[[1]]
#提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
dat=exprs(a)
#展示表达矩阵的前6行，看下数据是否经过log转换，一般数据在20以内，是经过log转换的，若有成百上千的数据，表示这个数据还需要进一步log处理。
head(dat)
#也可以输入上面这段代码，如果数据已经经log后，显示log2 transform not needed；如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NA
dat <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
#查看a的临床信息，为后面选择用于分组的变量做准备
pd=pData(a)
#查看pd
#View(pd)
#输出临床信息
write.csv(pd,'clinicalGSE65858.csv',row.names = TRUE)

#另外一种获取GPL文件方法
#如果没有下载过这个包，就下载一下，代码如下
#BiocManager::install("hgu133plus2.db")
#提取信息
library(illuminaHumanv4.db)
ls("package:illuminaHumanv4.db") #大致查看R包中的信息，寻找我们需要的SYMBOL
ids <- toTable(illuminaHumanv4SYMBOL)  #用toTable（）函数提取
head(ids) #查看提取内容
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")
#exp <- mutate(ex,probe_id=rownames(probe_exp))#将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
#exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id") # 合并数据
#exp2=exp2[!duplicated(probe_exp$symbol),]# 按照symbol列去重
#rownames(exp2)=exp2$symbol # 数据框probe_exp的行名变成symbol
#gene_exp_matrix <- na.omit(exp2 ) # 去空值
#输出文件
#write.table(gene_exp_matrix2,file = "Gastric.cancer.geneid.exprs.txt",sep = "\t",row.names=T,col.names = T)

#另外一种方式
#下载注释文件
#gpl <- getGEO('GPL10558', destdir="")
#我们将平台文件转为list格式，并赋值给gpl1，将gpl1保存为R文件，方便以后调用。
#gpl1<-Table(gpl)
#查看平台文件的列名，我们看到有ID和gene symbol,记住gene symbol在第几列，我们这里在第11列。
#colnames(Table(gpl))
#再了解一下平台文件的数据，当然这里大家可以直接选择gene symbol字段也行，主要了解一下symbol中基因symbol值是否只有一个。
#View(gpl1)
#我们这里的gene symbol字段中的symbol，有的基因就不止一个名称，///后面有重名，我们需要第一个名字，所以需要用字符分割包再处理一下,原理同上面处理title一样。
#提取平台文件gpl1中的ID和gene symbol，并赋值给probe2symbol_df
#gpl1$`Gene Symbol`=str_split(gpl1[,11],'///',simplify = T)[,1]
#probe2symbol_df<-gpl1[,c(1,11)]

probe2symbol_df <- ids
#将列名改为probe_id和symbol
#这两步是因为我懒，不想再调整代码，为了方便后面代码运行，我改的，大家也可以不改。
colnames(probe2symbol_df)=c('probe_id','symbol')
#查看symbol为唯一值的个数
length(unique(probe2symbol_df$symbol))
#查看每个基因出现n次的个数,我们可以看到，symbol出现一次的基因有7050个，出现2次有1493个。。。
table(sort(table(probe2symbol_df$symbol)))
#去掉没有注释symbol的探针（其实这里没有注释的探针数量即为上面出现次数最多的基因440，也就是说有440个探针没有symbol）
ids=probe2symbol_df[probe2symbol_df$symbol != '',]
#%in%用于判断是否匹配，
#注意这里dat是gset表达矩阵的数据，这一步就是把平台文件的ID和矩阵中的ID匹配。
ids=probe2symbol_df[probe2symbol_df$probe_id %in%rownames(dat),]
#取表达矩阵中可以与探针名匹配的那些，去掉无法匹配的表达数据
dat=dat[ids$probe_id,]
#ids新建mean这一列，列名为mean，同时对dat这个矩阵按行操作，取每一行（即每个样本在这个探针上的）的均值，将结果给到mean这一列的每一行，这里也可以改为中位值，median.
ids$mean=apply(dat,1,mean)
#即先按symbol排序，相同的symbol再按照均值从大到小排列，方便后续保留第一个值。
ids=ids[order(ids$symbol,ids$mean,decreasing = T),]
#将symbol这一列取出重复项，'!'为否，即取出不重复的项，去除重复的gene
#取median的话，这样就保留每个基因最大表达量结果.最后得到n个基因。
ids=ids[!duplicated(ids$symbol),]
#新的ids取出probe_id这一列，将dat按照取出的这一列，用他们的每一行组成一个新的dat
dat=dat[ids$probe_id,]
#把ids的symbol这一列中的每一行给dat作为dat的行名
rownames(dat)=ids$symbol
View(dat)
#但是我们看到矩阵最后一行，是没有symbol名字的，我们把他去掉，数字自己更改。
dim(dat)
#dat<-dat[-13437,]
#rownames(dat)[20911] = "A1BG"
#最后我们把结果保存。
#write.csv(dat,file="2GEO/GSE68465/GSE68465.csv")
#write.table(dat,file="GEO.txt", sep="\t", quote=F, row.names = TRUE)
write.table(data.frame(ID=rownames(dat),dat),file="GSE65858.txt", sep="\t", quote=F, row.names = F,col.names = T)


####1.10.meta####

setwd("E:/99other/31HNSCC/1data/GEO/META")
dir.create('results')
#
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(factoextra)
  library(sva)
})

#手动处理以下内容
cli1  <- fread("E:/99other/31HNSCC/1data/GEO/GSE27020-GPL96/clinical.GSE27020.txt") %>%
 select(Samples,OS,OS.time) %>%
 filter(OS.time > 0) %>%
 mutate(batch="GSE27020")

cli2  <- fread("E:/99other/31HNSCC/1data/GEO/GSE41613-GPL570/clinical.GSE41613.txt") %>%
 select(Samples,OS,OS.time) %>%
 filter(OS.time > 0) %>%
 mutate(batch="GSE41613")

cli3  <- fread("E:/99other/31HNSCC/1data/GEO/GSE65858-GPL10558/clinical.GSE65858.txt") %>%
  select(Samples,OS,OS.time) %>%
  filter(OS.time > 0) %>%
  mutate(batch="GSE65858")

cli <- rbind(cli1,cli2,cli3) %>% column_to_rownames("Samples")

data1 <- read.table('E:/99other/31HNSCC/1data/GEO/GSE27020-GPL96/GSE27020.txt',sep="\t",header=T,check.names=F,row.names = 1)
data2 <- read.table('E:/99other/31HNSCC/1data/GEO/GSE41613-GPL570/GSE41613.txt',sep="\t",header=T,check.names=F,row.names = 1)
data3 <- read.table('E:/99other/31HNSCC/1data/GEO/GSE65858-GPL10558/GSE65858.txt',sep="\t",header=T,check.names=F,row.names = 1)

genes<-Reduce(intersect,list(rownames(data1),rownames(data2), rownames(data3)))

data1 = data1[genes,]
data2 = data2[genes,]
data3 = data3[genes,]

data.orig <- cbind(data1,data2,data3) %>% .[,rownames(cli)]
data.combat <- ComBat(dat = data.orig, batch = cli$batch)

pca.orig <- prcomp(t(data.orig), scale = TRUE)
pca.combat <- prcomp(t(data.combat), scale = TRUE)

p.orig <- fviz_pca_ind(pca.orig,label=F, col.ind = cli$batch)
p.combat <- fviz_pca_ind(pca.combat,label=F,col.ind = cli$batch)

#save
ggsave(plot = p.orig,filename = "results/Fig.S1A.pdf",height = 5,width = 6)
ggsave(plot = p.combat,filename = "results/Fig.S1B.pdf",height = 5,width = 6)

write.table(data.combat,gzfile("results/meta.combat.txt.gz"),sep='\t',quote = F)
write.table(cli,"results/meta.cli.txt",sep='\t',quote = F)

####2.1.MOVICS####

setwd("E:/99other/31HNSCC/2MOVICS")
dir.create('results')
#devtools::install_github("xlucpu/MOVICS")
suppressPackageStartupMessages({
  library(tidyverse)
  library(MOVICS)
  library(data.table)
})
##### MOVICS 数据准备
mo.data   <- readRDS('E:/99other/31HNSCC/1data/TCGA/results/mo.data.rds')
tcga.cli  <- read.table('E:/99other/31HNSCC/1data/TCGA/results/tcga.cli.txt', header=T, 
                        sep="\t", check.names=F,row.names = 1) %>% select(fustat=OS,futime=OS.time)

#mad, sd, cox, freq和na.action
mo.data$mRMA <- mo.data$mRNA.exp %>% getElites(method= "mad",na.action = "rm",elite.num = 1500) %>% pluck('elite.dat') %>% 
  getElites(method= "cox",surv.info = tcga.cli,p.cutoff = 0.05) %>% pluck('elite.dat')

mo.data$lncRNA <- mo.data$lncRNA.exp %>% getElites(method= "mad",na.action = "rm",elite.num = 1500) %>% pluck('elite.dat') %>% 
  getElites(method= "cox",surv.info = tcga.cli,p.cutoff = 0.05) %>% pluck('elite.dat')

mo.data$miRMA <- mo.data$miRMA.exp %>% getElites(method= "mad",na.action = "rm",elite.num = 1500) %>% pluck('elite.dat') %>% 
  getElites(method= "cox",surv.info = tcga.cli,p.cutoff = 0.05) %>% pluck('elite.dat')

mo.data$meth <- mo.data$meth.beta %>% getElites(method= "mad",na.action = "rm",elite.num = 1500) %>% pluck('elite.dat') %>% 
  getElites(method= "cox",surv.info = tcga.cli,p.cutoff = 0.05) %>% pluck('elite.dat')

mo.data$mut <- mo.data$mut.status %>% as.matrix %>% getElites(method= "freq",na.action = "rm",elite.pct = 0.05) %>% pluck('elite.dat')

mo.data <- mo.data[c("mRMA","lncRNA","miRMA","meth","mut")]
######Fig.S2-CLUSTER NUMBER
optk <- getClustNum(data = mo.data,
                    is.binary = c(F,F,F,F,T), 
                    try.N.clust = 2:8, 
                    fig.path = 'results',
                    fig.name = "Fig.S2")
#### Fig.2C-CONSENSUS HEATMAP
moic.res.list <- getMOIC(data = mo.data, 
                         N.clust = optk$N.clust,
                         type = c("gaussian", "gaussian", "gaussian", "gaussian" ,"binomial"))

cmoic <- getConsensusMOIC(moic.res.list = moic.res.list,
                          distance = "euclidean", 
                          linkage  = "average",
                          fig.path = 'results',
                          fig.name = "Fig.2C")
##savehttp://127.0.0.1:17095/graphics/plot_zoom_png?width=1149&height=778
saveRDS(moic.res.list,"E:/99other/31HNSCC/2MOVICS/results/moic.res.list.rds",compress =F)
saveRDS(cmoic,"E:/99other/31HNSCC/2MOVICS/results/cmoic.rds",compress =F)
moic.res.list<-read_rds('E:/99other/31HNSCC/2MOVICS/results/moic.res.list.rds')
cmoic<-read_rds('E:/99other/31HNSCC/2MOVICS/results/cmoic.rds')
####### Fig.S3-SILHOUETTE
getSilhouette(sil = cmoic$sil,
              fig.path = 'results',
              fig.name = "Fig.S3")
##### Fig.2A-COMPREHENSIVE HEATMAP OF CONSENSUSMOIC
# convert beta value to M value for stronger signal
indata <- mo.data
indata$meth <- log2(indata$meth / (1 - indata$meth))
# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,2,NA),
                     centerFlag = c(T,T,T,T,F),
                     scaleFlag  = c(T,T,T,T,F))

feat   <- moic.res.list$iClusterBayes$feat.res
feat1  <- feat[which(feat$dataset == "mRMA"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "miRMA"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "meth"),][1:10,"feature"]
feat5  <- feat[which(feat$dataset == "mut"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4,feat5)

mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
miRMA.col  <- c("yellow", "black"  , "blue")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, miRMA.col ,meth.col, mut.col)

getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","miRMA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T),
             legend.name   = c("mRNA","lncRNA","miRMA","M value","Mutated"),
             clust.res     = cmoic$clust.res,
             annRow        = annRow,
             color         = col.list,
             width         = 10,
             height        = 5,
             fig.path      = 'results',
             fig.name      = "Fig.2A")
##### Fig.2B-10 CLUSTERING METHODS
cmoic$clust.res %>% select(clust) %>% write.table("E:/99other/31HNSCC/2MOVICS/results/tcga.clust.txt",sep="\t",quote = F)

clust.res <- cmoic$clust.res$clust

for (moic.res in moic.res.list){
  clust.res <- cbind(clust.res,moic.res$clust.res$clust)
}
colnames(clust.res) <- c('Subtype',names(moic.res.list))

col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627","#023E8A", "#9D4EDD")
names(col) <- 1:length(col)

clust.res = as.data.frame(clust.res) %>% arrange(Subtype)
anno = ComplexHeatmap::HeatmapAnnotation(df=clust.res,show_legend = F,col=list(Subtype= col))
library(ComplexHeatmap)
pdf(file = 'E:/99other/31HNSCC/2MOVICS/results/Fig.2B.pdf',width = 5,height = 3)
Heatmap(matrix(nrow = 0, ncol = length(clust.res$Subtype)),top_annotation =anno)
dev.off()
#### Fig.2D-KAPLAN-MEIER CURVE OF CONSENSUSMOIC
surv <- compSurv(moic.res         = cmoic,
                 surv.info        = tcga.cli,
                 convt.time       = "m",
                 surv.median.line = "h", 
                 xyrs.est         = c(5,10),
                 fig.path      = 'results',
                 fig.name         = "Fig.2D")

####3.1.landscape and validation####


setwd("E:/99other/31HNSCC/3LandscapeValidation")
dir.create('results')
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(GSVA)
  library(RTN)
  library(MOVICS)
  library(IOBR)
  library(limma)
  library(GSEABase)
  library(GSVA)
  library(pheatmap)
  library(reshape2)
  library(ggpubr)
  library(readxl)
  library(stringr)
})

tcga.data <- read.table("E:/99other/31HNSCC/1data/TCGA/results/tcga.mRNA.exp.txt.gz",check.names = F)
tcga.clust <- read.table("E:/99other/31HNSCC/2MOVICS/results/tcga.clust.txt")  %>% arrange(clust)
sam.order <- tcga.clust %>% rownames()

#读入
mo.data <- readRDS('E:/99other/31HNSCC/1data/TCGA/results/mo.data.rds')
cmoic <- readRDS('E:/99other/31HNSCC/2MOVICS/results/cmoic.rds')

# Fig.3D-NTP HEATMAP
meta.data <- read.table('E:/99other/31HNSCC/1data/GEO/META/results/meta.combat.txt.gz')
meta.cli <- read.table('E:/99other/31HNSCC/1data/GEO/META/results/meta.cli.txt') %>% dplyr::select(fustat=OS,futime=OS.time)

path.temp <- 'results/TEMP'
dir.create(path.temp)

runDEA(dea.method = "limma",
       res.path   = path.temp,
       expr       = tcga.data,
       moic.res   = cmoic,
       prefix     = "TCGA")

marker.up <- runMarker(moic.res      = cmoic,
                       dea.method    = "limma",
                       prefix        = "TCGA",
                       dat.path      = path.temp,
                       res.path      = path.temp, 
                       p.cutoff      = 0.05, 
                       p.adj.cutoff  = 0.05, 
                       dirct         = "up", 
                       n.marker      = 1000,  #根据分群的数量，调整选择合适的总数
                       doplot        = F)
marker.up$templates %>% write.table('results/marker.up.txt',row.names=F,sep="\t",quote = F)

meta.ntp.pred <- runNTP(expr = meta.data,
                        templates  = marker.up$templates,
                        scaleFlag  = TRUE,
                        centerFlag = TRUE,
                        doPlot     = TRUE,
                        fig.path   = "results",
                        fig.name   = "Fig.3D") 
unlink(path.temp,recursive = T)

# Fig.3E-KAPLAN-MEIER CURVE
compSurv(moic.res         = meta.ntp.pred,
         surv.info        = meta.cli,
         convt.time       = "m",
         surv.median.line = "hv", 
         fig.path         = "results",
         fig.name         = "Fig.3E")

# Fig.3F-H CONSISTENCY HEATMAP
meta.pam.pred <- runPAM(train.expr  = tcga.data,
                        moic.res    = cmoic,
                        test.expr   = meta.data)

tcga.ntp.pred <- runNTP(expr      = tcga.data,
                        templates = marker.up$templates,
                        doPlot    = FALSE)

tcga.pam.pred <- runPAM(train.expr  = tcga.data,
                        moic.res    = cmoic,
                        test.expr   = tcga.data)

runKappa(subt1     = cmoic$clust.res$clust,
         subt2     = tcga.ntp.pred$clust.res$clust,
         subt1.lab = "CSs",
         subt2.lab = "NTP",
         fig.path  = "results",
         fig.name  = "Fig.3F")

runKappa(subt1     = cmoic$clust.res$clust,
         subt2     = tcga.pam.pred$clust.res$clust,
         subt1.lab = "CSs",
         subt2.lab = "PAM",
         fig.path  = "results",
         fig.name  = "Fig.3G")

runKappa(subt1     = meta.ntp.pred$clust.res$clust,
         subt2     = meta.pam.pred$clust.res$clust,
         subt1.lab = "NTP",
         subt2.lab = "PAM",
         fig.path  = "results",
         fig.name  = "Fig.3H")


####3.2.prepare####

setwd("E:/99other/31HNSCC/3LandscapeValidation")
markers <- read.table("E:/99other/31HNSCC/3LandscapeValidation/results/marker.up.txt",header = T)

tcga.data <- read.table('E:/99other/31HNSCC/1data/TCGA/results/tcga.mRNA.exp.txt.gz',check.names = F)
tcga.cli <- read.table('E:/99other/31HNSCC/1data/TCGA/results/tcga.cli.txt', 
                       header=T, sep="\t", check.names=F,row.names = 1)

rt = tcga.data[markers[,1],]


samples = intersect(colnames(rt),rownames(tcga.cli))

rt = rt[,samples]
tcga.cli = tcga.cli[samples,]


uniSigExpTime = cbind(tcga.cli[,c(2,1)],t(rt))

colnames(uniSigExpTime)[1]="futime"
colnames(uniSigExpTime)[2]="fustat"

write.table(data.frame(ID=rownames(uniSigExpTime),uniSigExpTime), file="uniSigExpTime.txt", sep="\t", quote=F, col.names=T,row.names = F)


tcga.cluster <- read.table('tcga.clust.txt',check.names = F)

colnames(tcga.cluster) = "CS"

write.table(data.frame(ID=rownames(tcga.cluster),tcga.cluster), file="tcga.clust.new.txt", sep="\t", quote=F, col.names=T,row.names = F)


####3.3.PCA####

#引用包
library(limma)
library(Rtsne)
library(umap)
library(ggplot2)

expFile="uniSigExpTime.txt"           #表达数据文件
clusterFile="tcga.clust.new.txt"      #分型的结果文件
setwd("E:/99other/31HNSCC/3LandscapeValidation")    #设置工作目录

#读取输入文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)

rt = rt[,-c(2,3)]

rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=t(data)

############PCA分析
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)

#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
cluster=as.vector(cluster[,1])

#设置分型的颜色
bioCol=c("#2EC4B6","#E71D36","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
clusterCol=bioCol[1:length(levels(factor(cluster)))]

#绘制图形
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], cluster=cluster)
PCA.mean=aggregate(PCA[,1:2], list(cluster=PCA$cluster), mean)
PCA$cluster = factor(PCA$cluster)
pdf(file="PCA.pdf", width=5.5, height=4.75)
p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color=cluster, shape=cluster)) +
  scale_colour_manual(name="cluster", values =clusterCol)+
  theme_bw()+
  labs(title ="PCA")+
  theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$cluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()


############UMAP分析
umapOut=umap(data)
umap=data.frame(UMAP1=umapOut$layout[,1], UMAP2=umapOut$layout[,2], cluster=cluster)
umap.mean=aggregate(umap[,1:2], list(cluster=umap$cluster), mean)	
umap$cluster = factor(umap$cluster)
#绘制分型的UAMP图
pdf(file="UMAP.pdf", width=5.5, height=4.75)       #保存输入出文件
p=ggplot(data=umap, aes(UMAP1, UMAP2)) + geom_point(aes(color=cluster, shape=cluster)) +
  scale_colour_manual(name="cluster",  values =clusterCol)+
  theme_bw()+
  labs(title ="UAMP")+
  theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))+
  annotate("text",x=umap.mean$UMAP1, y=umap.mean$UMAP2, label=umap.mean$cluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()


############tSNE分析
tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
tsne=data.frame(tSNE1=tsneOut$Y[,1], tSNE2=tsneOut$Y[,2], cluster=cluster)
tSNE.mean=aggregate(tsne[,1:2], list(cluster=tsne$cluster), mean)	
tsne$cluster = factor(tsne$cluster)
#绘制分型的tSNE图
pdf(file="tSNE.pdf", width=5.5, height=4.75)       #保存输入出文件
p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color=cluster, shape=cluster)) +
  scale_colour_manual(name="cluster",  values =clusterCol)+
  theme_bw()+
  labs(title ="tSNE")+
  theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))+
  annotate("text",x=tSNE.mean$tSNE1, y=tSNE.mean$tSNE2, label=tSNE.mean$cluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()




####3.4.clusterDiff####

#引用包
library(limma)
library(reshape2)
library(ggpubr)

expFile="uniSigExpTime.txt"           #表达数据文件
clusterFile="tcga.clust.new.txt"      #分型的结果文件
setwd("E:/99other/31HNSCC/3LandscapeValidation")    #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
exp = exp[,-c(1,2)]
#colnames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(exp))
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#data=avereps(data)
#data=t(data)

genes = fread("module.gene.txt",header = F)

data=data[,as.vector(genes[,1])[[1]]]


#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(data), row.names(cluster))
expClu=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
expClu$CS = factor(expClu$CS)

#提取差异显著的基因
sigGene=c()
for(i in colnames(expClu)[1:(ncol(expClu)-1)]){
  if(sd(expClu[,i])<0.001){next}
  if(length(levels(factor(expClu[,"CS"])))>2){
    test=kruskal.test(expClu[,i] ~ expClu[,"CS"])
  }else{
    test=wilcox.test(expClu[,i] ~ expClu[,"CS"])
  }
  pvalue=test$p.value
  if(pvalue<0.05){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "CS")
expClu=expClu[,sigGene]

#把数据转换成ggplot2输入文件
data=melt(expClu, id.vars=c("CS"))
colnames(data)=c("CS", "Gene", "Expression")

#设置图形颜色
bioCol=c("#2EC4B6","#E71D36","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"CS"])))]

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color="CS",
            xlab="",
            ylab="Gene expression",
            legend.title="CS",
            palette = bioCol,
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=CS),
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出箱线图
pdf(file="geneboxplot.diff.pdf", width=10, height=4.5)
print(p1)
dev.off()






####3.5.heatmap####

library(pheatmap)       #引用包
expFile="uniSigExpTime.txt"           #表达数据文件
clusterFile="tcga.clust.new.txt"      #分型的结果文件
cliFile="clinical5.txt"           #临床数据文件
setwd("E:/99other/31HNSCC/3LandscapeValidation")    #设置工作目录

#读取表达数据文件
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp = exp[,-c(1,2)]
genes = fread("module.gene.txt",header = F)
exp=exp[,as.vector(genes[,1])[[1]]]


#colnames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(exp))
#exp=t(exp)
#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#合并表达和分型数据
sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample, , drop=F]
cluster=cluster[sameSample, , drop=F]
expCluster=cbind(exp, cluster)
#Project=gsub("(.*?)\\_.*", "\\1", rownames(expCluster))
rownames(expCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(expCluster))
#expCluster=cbind(expCluster, Project)

#合并临床数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))
sameSample=intersect(row.names(expCluster), row.names(cli))
expCluster=expCluster[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expCluster, cli)

#提取热图的数据

data$CS = factor(data$CS)
data=data[order(data$CS),]
Type=data[,((ncol(exp)+1):ncol(data))]
data=t(data[,1:ncol(exp)])

#定义热图注释的颜色
bioCol=c("#2EC4B6","#E71D36","#7CC767","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
prgCluCol=bioCol[1:length(levels(factor(Type$CS)))]
names(prgCluCol)=levels(factor(Type$CS))
ann_colors[["CS"]]=prgCluCol

#热图可视化
pdf("heatmap.pdf", width=8, height=5.5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("#44B1C9",5), "white", rep("#DA4A35",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         show_colnames=F,
         #scale="row",
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()


####3.6.ssGSEA####

#引用包
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="tcga.mRNA.TPM.txt"               #表达数据文件
clusterFile="tcga.clust.new.txt"      #分型的结果文件
gmtFile="immune.gmt"              #免疫基因集文件
setwd("E:/99other/31HNSCC/3LandscapeValidation")    #设置工作目录

#读取表达数据文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
#rownames(rt)=rt[,1]
#exp=rt[,2:ncol(rt)]
#dimnames=list(rownames(exp),colnames(exp))
#data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#data=avereps(data)
data = rt
#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA分析
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#对ssGSEA打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#输出ssGSEA打分结果
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并(将分型的结果与免疫细胞的打分进行合并)
ssgseaScore=t(ssgseaScore)
#row.names(ssgseaScore)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(ssgseaScore))
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

#把数据转换成ggplot2输入文件
data=melt(scoreCluster, id.vars=c("CS"))
colnames(data)=c("CS", "Immune", "Fraction")
data$CS=factor(data$CS)
#绘制箱线图
bioCol=c("#2EC4B6","#E71D36","#7CC767","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"CS"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="CS",
            xlab="",
            ylab="Immune infiltration",
            legend.title="CS",
            palette=bioCol)
p=p+rotate_x_text(50)

#保存图形
pdf(file="boxplot.pdf", width=11, height=5)
p+stat_compare_means(aes(group=CS),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()


####3.7.GSVA####

#引用包
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)

expFile="tcga.mRNA.TPM.txt"               #表达数据文件
clusterFile="tcga.clust.new.txt"      #分型的结果文件
gmtFile="c2.cp.kegg.symbols.gmt"    #基因集文件
setwd("E:/99other/31HNSCC/3LandscapeValidation")    #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
# rownames(rt)=rt[,1]
# exp=rt[,2:ncol(rt)]
# dimnames=list(rownames(exp), colnames(exp))
# data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
# data=avereps(data)
data = rt
#GSVA分析
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)

#读取分型的结果
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并(将分型和GSVA的结果进行合并)
gsvaResult=t(gsvaResult)
#row.names(gsvaResult)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(gsvaResult))
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)

#Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
#gsvaCluster=cbind(gsvaCluster, Project)

#设置比较组

gsvaCluster$CS=gsub(1,"A",gsvaCluster$CS)
gsvaCluster$CS=gsub(2,"B",gsvaCluster$CS)

adj.P.Val.Filter=0.05
allType=as.vector(gsvaCluster$CS)
comp=combn(levels(factor(allType)), 2)

#对比较组进行循环, 观察哪些通路在不同的分型之间是具有差异的
for(i in 1:ncol(comp)){
  #样品分组
  treat=gsvaCluster[gsvaCluster$CS==comp[2,i],]
  con=gsvaCluster[gsvaCluster$CS==comp[1,i],]
  data=rbind(con, treat)
  #对通路进行差异分析
  Type=as.vector(data$CS)
  ann=data[,c(ncol(data), (ncol(data)-1))]
  data=t(data[,-c((ncol(data)-1), ncol(data))])
  design=model.matrix(~0+factor(Type))
  colnames(design)=levels(factor(Type))
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  #输出所有通路的差异情况
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  #输出差异显著的通路
  diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  
  #设置热图注释的颜色
  bioCol=c("#2EC4B6","#E71D36","#7CC767","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  ann_colors=list()
  m6aCluCol=bioCol[1:length(levels(factor(allType)))]
  names(m6aCluCol)=levels(factor(allType))
  ann_colors[["CS"]]=m6aCluCol[c(comp[1,i], comp[2,i])]
  
  #绘制差异通路热图
  termNum=15     #设置显示通路的数目
  term1 = subset(diffSig, logFC>0)
  diffterm1=as.vector(rownames(term1))
  diffterm1Length=length(diffterm1)
  if(diffterm1Length<termNum){termNum=diffterm1Length}
  hmGene1=diffterm1[1:termNum]
  
  termNum=15     #设置显示通路的数目
  term2 = subset(diffSig, logFC<0)
  diffterm2=as.vector(rownames(term2))
  diffterm2Length=length(diffterm2)
  if( diffterm2Length<termNum){termNum=diffterm2Length}
  hmGene2=diffterm2[1:termNum]
  
  hmGene = c(hmGene1,hmGene2)
  hmExp=data[hmGene,]
  
  
  ann$CS = gsub("A","1",ann$CS)
  ann$CS = gsub("B","2",ann$CS)
  ann$CS = factor(ann$CS)
  
  ann_colors=list()
  m6aCluCol=bioCol[1:length(levels(factor(allType)))]
  names(m6aCluCol)=levels(factor(c("1","2")))
  ann_colors[["CS"]]=m6aCluCol[c("1", "2")]
  
  Type = gsub("A","1",Type)
  Type = gsub("B","2",Type)
  
  pdf(file=paste0(contrast,".heatmap.pdf"), width=10, height=4.5)
  pheatmap(hmExp, 
           annotation=ann,
           annotation_colors = ann_colors,
           color = colorRampPalette(c(rep("#1E4A93",2), "white", rep("#D21E1F",2)))(50),
           cluster_cols =F,
           show_colnames = F,cluster_rows = T,
           gaps_col=as.vector(cumsum(table(Type))),
           #scale="row",
           fontsize = 8,
           fontsize_row=6,
           fontsize_col=8)
  dev.off()
}




####4.1.Development####


setwd("E:/99other/31HNSCC/4Development.CLMS")
dir.create('results')
#
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(survival)
  library(survminer)
})
## cox_batch define
cox_batch<-function(dat,time,event){
  coxRun<-function(dat){
    library(survival)
    colnames(dat)=c('time','status','AS')  
    dat=dat[which(!is.na(dat[,1])&!is.na(dat[,3])&!is.na(dat[,2])),]
    #print(nrow(dat))
    if(nrow(dat)<10){
      print(paste0('Sample Num is small:',nrow(dat)))
      return(c(NA,NA,NA,NA))
    }
    #if(quantile(dat[,3])['25%']==quantile(dat[,3])['50%']) return(c(NA,NA,NA,NA))
    fmla <- as.formula("Surv(time, status) ~AS")
    if(table(dat[,2])[1]>1&table(dat[,2])[2]>1){
      cox <- survival::coxph(fmla, data = dat)
      re=c(summary(cox)[[7]][5],summary(cox)[[7]][2],summary(cox)[[8]][3],summary(cox)[[8]][4])
      return(re)
    }else{
      return(c(NA,NA,NA,NA))
    }
  }
  t.inds=which(!is.na(time)&!is.na(event))
  dat1=dat[,t.inds]
  os=time[t.inds]
  ev=event[t.inds]
  
  ct=sum(ev%in%c(0,1))
  if(ct!=length(ev)){
    print('event must be 0(alive) or 1(dead)')
    return(NULL)
  }
  
  res=t(apply(dat1, 1, function(x){
    ep=as.numeric(as.character(x))
    ind2=which(!is.na(ep))
    # print(length(ind2))
    if(length(ind2)>1){
      os1=os[ind2]
      ev1=ev[ind2]
      ep1=ep[ind2]
      return(coxRun(data.frame(os1,ev1,ep1)))
    }else{
      return(c(NA,NA,NA,NA))
    }
  }))
  colnames(res)=c('pvalue','HR','HR.95L','HR.95H')
  row.names(res)=row.names(dat1)
  return(as.data.frame(res))
}
## markers + cox 筛选基因
markers <- read.table("E:/99other/31HNSCC/3LandscapeValidation/results/marker.up.txt",header = T)

tcga.data <- read.table('E:/99other/31HNSCC/1data/TCGA/results/tcga.mRNA.exp.txt.gz',check.names = F)
tcga.cli <- read.table('E:/99other/31HNSCC/1data/TCGA/results/tcga.cli.txt', 
                       header=T, sep="\t", check.names=F,row.names = 1)

meta.data <- read.table('E:/99other/31HNSCC/1data/GEO/META/results/meta.combat.txt.gz',check.names = F)
meta.cli <- read.table('E:/99other/31HNSCC/1data/GEO/META/results/meta.cli.txt')

#meta.data <- read.table('E:/99other/31HNSCC/1data/GEO/GSE157010-GPL570/GSE157010.exp.txt.gz')
#meta.cli <- read.table('E:/99other/31HNSCC/1data/GEO/GSE157010-GPL570/cli.GSE157010.txt') %>% dplyr::select(fustat=OS,futime=OS.time)

#meta.data <- read.table('E:/99other/31HNSCC/1data/GEO/GSE74777-GPL17586/GSE74777.exp.txt.gz')
#meta.cli <- read.table('E:/99other/31HNSCC/1data/GEO/GSE74777-GPL17586/cli.GSE74777.txt') %>% dplyr::select(fustat=OS,futime=OS.time)

#meta.data <- read.table('E:/99other/31HNSCC/1data/GEO/GSE73403-GPL6480/GSE73403.exp.txt.gz')
#meta.cli <- read.table('E:/99other/31HNSCC/1data/GEO/GSE73403-GPL6480/cli.GSE73403.txt') %>% dplyr::select(fustat=OS,futime=OS.time)

#colnames(meta.cli)= c("OS","OS.time")

selectGenes <- markers$probe %>% intersect(row.names(tcga.data)) %>% intersect(row.names(meta.data))

tcga.cox <- cox_batch(dat = tcga.data[selectGenes,rownames(tcga.cli)],
                      time = tcga.cli$OS.time,
                      event = tcga.cli$OS)
meta.cox <- cox_batch(dat = meta.data[selectGenes,rownames(meta.cli)],
                      time = meta.cli$OS.time,
                      event = meta.cli$OS)
#这里以p《0.001筛选预后相关基因
tcga.filtedGenes <- tcga.cox %>% filter(pvalue<0.05) %>% row.names()
meta.filtedGenes <- meta.cox %>% filter(pvalue<0.05) %>% row.names()
com.gene <- intersect(tcga.filtedGenes,meta.filtedGenes)
length(com.gene)
write.table(com.gene,"results/com.gene.txt",quote = F,sep='\t')

tcga.cox[com.gene,] %>% write.table("results/tcga.cox.txt",quote = F,sep='\t')
meta.cox[com.gene,] %>% write.table("results/meta.cox.txt",quote = F,sep='\t')
#
dir.create("results/run.input")

tcga.data[com.gene,] %>% write.table(gzfile('results/run.input/tcga_exp.txt.gz'),quote = F,sep='\t')
meta.data[com.gene,] %>% write.table(gzfile('results/run.input/meta_exp.txt.gz'),quote = F,sep='\t')

tcga.cli %>% mutate(Cohort="TCGA") %>% 
  rownames_to_column("sample") %>% 
  dplyr::select(sample,Cohort,OS,OS.time) %>% 
  write.table("results/run.input/tcga.cli.txt",row.names = F,quote = F,sep='\t')

meta.cli %>% mutate(Cohort="META") %>% 
  rownames_to_column("sample") %>% 
  dplyr::select(sample,Cohort,OS,OS.time) %>% 
  write.table("results/run.input/meta.cli.txt",row.names = F,quote = F,sep='\t')
#
dir.create("results/run.results")

### 将results/run.input 里的文件输入到 http://www.sxdyc.com/prognosticmodel
### 将运行结果解压下载后解压至到 result/run.results

#或以下
# 加载需要使用的R包
library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)

# 设置工作路径
work.path <- "E:/99other/31HNSCC/4Development.CLMS/1model"; setwd(work.path) 

# 设置其他路径
code.path <- file.path(work.path, "Codes") # 存放脚本
data.path <- file.path(work.path, "InputData") # 存在输入数据（需用户修改）
res.path <- file.path(work.path, "Results") # 存放输出结果
fig.path <- file.path(work.path, "Figures") # 存放输出图片

# 如不存在这些路径则创建路径
if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

# 加载模型训练以及模型评估的脚本
source(file.path(code.path, "ML.R"))

# 选择最后生成的模型类型：panML代表生成由不同算法构建的模型； multiCox表示抽取其他模型所用到的变量并建立多变量cox模型
FinalModel <- c("panML", "multiCox")[2]

## Training Cohort
# 训练集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与测试集保持相同类型，表达谱需有一定变异性，以免建模过程报错）
Train_expr <- read.table("E:/99other/31HNSCC/4Development.CLMS/results/run.input/tcga_exp.txt.gz", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 训练集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Train_surv <- read.table("E:/99other/31HNSCC/4Development.CLMS/results/run.input/tcga.cli.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]
Train_surv = Train_surv[,-1]

## Validation Cohort
# 测试集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与训练集保持相同类型）
Test_expr <- read.table("E:/99other/31HNSCC/4Development.CLMS/results/run.input/meta_exp.txt.gz", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 测试集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Test_surv <- read.table("E:/99other/31HNSCC/4Development.CLMS/results/run.input/meta.cli.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

# 提取相同基因
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))

#comgene = comgene[comgene != c("MAMDC2","ITGA3")]
#"MAMDC2","ITGA3"

Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因

# 按队列对数据分别进行标准化（根据情况调整centerFlags和scaleFlags）
## data: 需要表达谱数据（行为样本，列为基因） 
## cohort：样本所属队列，为向量，不输入值时默认全表达矩阵来自同一队列
## centerFlag/scaleFlags：是否将基因均值/标准差标准化为1；
##        默认参数为NULL，表示不进行标准化；
##        为T/F时，表示对所有队列都进行/不进行标准化
##        输入由T/F组成的向量时，按顺序对队列进行处理，向量长度应与队列数一样
##        如centerFlags = c(F, F, F, T, T)，表示对第4、5个队列进行标准化，此时flag顺序应当与队列顺序一致
##        如centerFlags = c("A" = F, "C" = T, "B" = F)，表示对队列C进行标准化，此时不要求flag顺序与data一致
Train_set = scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
names(x = split(as.data.frame(Test_expr), f = Test_surv$Cohort)) # 注意测试集标准化顺序与此一致
Test_set = scaleData(data = Test_expr, cohort = Test_surv$Cohort, centerFlags = T, scaleFlags = T)
# summary(apply(Train_set, 2, var))
# summary(apply(Test_set, 2, var))
# lapply(split(as.data.frame(Test_set), Test_surv$Cohort), function(x) summary(apply(x, 2, var))) # 测试scale结果

# Model training and validation

## method list
# 此处记录需要运行的模型，格式为：算法1名称[算法参数]+算法2名称[算法参数]
# 目前仅有StepCox和RunEnet支持输入算法参数
methods <- read.xlsx(file.path(code.path, "41467_2022_28421_MOESM4_ESM.xlsx"), startRow = 2)$Model
methods <- gsub("-| ", "", methods)

## Train the model
min.selected.var <- 3 # 筛选变量数目的最小阈值
timeVar = "OS.time"; statusVar = "OS" # 定义需要考虑的结局事件，必须出现在Train_surv以及Test_surv中

## Pre-training 
Variable = colnames(Train_expr)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))
preTrain.method

set.seed(seed = 111) # 设置建模种子，使得结果可重复
preTrain.var <- list()
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, # 机器学习方法
                                 Train_expr = Train_set, # 训练集有潜在预测价值的变量
                                 Train_surv = Train_surv, # 训练集生存数据
                                 mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                 classVar = classVar) # 用于训练的生存变量，必须出现在Train_surv中
}
preTrain.var[["simple"]] <- colnames(Train_expr)

model <- list() # 初始化模型结果列表
set.seed(seed =111) # 设置建模种子，使得结果可重复
for (method in methods){ # 循环每一种方法组合
  # method <- "CoxBoost+survivalSVM" # [举例]若遇到报错，请勿直接重头运行，可给method赋值为当前报错的算法来debug
  cat(match(method, methods), ":", method, "\n") # 输出当前方法
  method_name = method # 本轮算法名称
  method <- strsplit(method, "\\+")[[1]] # 各步骤算法名称
  
  if (length(method) == 1) method <- c("simple", method)
  
  selected.var = preTrain.var[[method[1]]]
  # 如果筛选出的变量小于阈值，则该算法组合无意义，置空（尤其针对以RSF筛选变量的情况，需在ML脚本中尝试调参）
  if (length(selected.var) <= min.selected.var) {
    model[[method_name]] <- NULL
  } else {
    
    #selected.var1 = selected.var
    #selected.var = selected.var1[c(2,3,5,1,6)]
    model[[method_name]] <- RunML(method = method[2], # 用于构建最终模型的机器学习方法
                                  Train_expr = Train_expr[, selected.var], # 训练集有潜在预测价值的变量
                                  Train_surv = Train_surv, # 训练集生存数据
                                  mode = "Model",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                  classVar = classVar)  # 用于训练的生存变量，必须出现在Train_surv中
  }
  
  # 如果最终筛选出的变量小于阈值，则该算法组合也无意义，置空
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}

#EXPPP = Train_expr[, selected.var]
#write.table(data.frame(ID=rownames(cor(t(EXPPP))),cor(t(EXPPP))),file="cor.txt", sep="\t", quote=F, row.names = F,col.names = T)

saveRDS(model, file.path(res.path, "model.rds")) # 保存所有模型输出

# 当要求最终模型为多变量cox时，对模型进行更新
if (FinalModel == "multiCox"){
  coxmodel <- lapply(model, function(fit){ # 根据各算法最终获得的变量，构建多变量cox模型，从而以cox回归系数和特征表达计算单样本风险得分
    tmp <- coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                 data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
    return(tmp)
  })
}
saveRDS(coxmodel, file.path(res.path, "coxmodel.rds")) # 保存最终以多变量cox拟合所筛选变量的模型

## Evaluate the model

# 读取已保存的模型列表（请根据需要调整）
model <- readRDS(file.path(res.path, "model.rds")) # 若希望使用各自模型的线性组合函数计算得分，请运行此行
# model <- readRDS(file.path(res.path, "coxmodel.rds")) # 若希望使用多变量cox模型计算得分，请运行此行

methodsValid <- names(model) # 取出有效的模型（变量数目小于阈值的模型视为无效）

# 根据给定表达量计算样本风险评分,测试集
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]], 
                                    new_data = Test_set,
                                    type = "lp") # 同原文，使用linear Predictor计算得分
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) # 输出风险评分文件

# 根据给定表达量计算样本风险评分,训练集
RS_list2 <- list()
for (method in methodsValid){
  RS_list2[[method]] <- CalRiskScore(fit = model[[method]], 
                                     new_data = Train_set,
                                     type = "lp") # 同原文，使用linear Predictor计算得分
}
RS_mat2 <- as.data.frame(t(do.call(rbind, RS_list2)))

rownames(RS_mat2) = substr(rownames(RS_mat2),1,12)

RS_mat3 = rbind(RS_mat2,RS_mat)

write.table(RS_mat3, file.path(res.path, "riskscore_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) # 输出风险评分文件


# 提取所筛选的变量（列表格式）
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]]) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
}

# 提取所筛选的变量（数据框格式）
fea_df <- lapply(model, function(fit){ data.frame(ExtractVar(fit)) }) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"  # 数据框有两列，包含算法以及算法所筛选出的变量
write.table(fea_df, file.path(res.path, "fea_df.txt"),sep = "\t", row.names = F, col.names = T, quote = F)

# 对各模型计算C-index
Cindexlist <- list()
for (method in methodsValid){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], # 预后模型
                                  Test_expr = Test_set, # 测试集预后变量，应当包含训练集中所有的变量，否则会报错
                                  Test_surv = Test_surv, # 训练集生存数据，应当包含训练集中所有的变量，否则会报错
                                  Train_expr = Train_set, # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
                                  Train_surv = Train_surv, # 若需要同时评估训练集，则给出训练集生存数据，否则置NULL
                                  Train_name = "TCGA", # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
                                  #Train_expr = NULL,
                                  #Train_surv = NULL, 
                                  cohortVar = "Cohort", # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
                                  timeVar = timeVar, # 用于评估的生存时间，必须出现在Test_surv中；这里是OS.time
                                  statusVar = statusVar) # 用于评估的生存状态，必须出现在Test_surv中；这里是OS
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = T, col.names = T, quote = F)

# Plot

Cindex_mat <- read.table(file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T) # 计算每种算法在所有队列中平均C-index，并降序排列
Cindex_mat <- Cindex_mat[names(avg_Cindex), ] # 对C-index矩阵排序
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数
fea_sel <- fea_list[[rownames(Cindex_mat)[1]]] # 最优模型（即测试集[或者训练集+测试集]C指数均值最大）所筛选的特征

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") # 设置绘图时的队列颜色
names(CohortCol) <- colnames(Cindex_mat)
CohortCol = CohortCol[c(1,2)]
# 调用简易绘图函数
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, # 主矩阵
                    avg_Cindex = avg_Cindex, # 侧边柱状图
                    CohortCol = CohortCol, # 列标签颜色
                    barCol = "steelblue", # 右侧柱状图颜色
                    col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 热图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类

pdf(file.path(fig.path, "heatmap of cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 6, height = cellheight * nrow(Cindex_mat) * 0.55)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") # 热图注释均放在右侧
invisible(dev.off())


##################################### KM曲线函数定义
mycolor <- ggsci::pal_npg(palette = c("nrc"), alpha =1)(8)
ggsurvplotKM<-function(dat, title = 'Groups', 
                       lables = c(), col = mycolor,
                       risk.table = TRUE,
                       tables.height = 0.25) {
  # dat：数据框，行为样本，列为时间、状态以及分组
  # the color palette to be used. Allowed values include "hue" for the default hue color scale; "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...; or custom color palette e.g. c("blue", "red"); and scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty". See details section for more information. Can be also a numeric vector of length(groups); in this case a basic color palette is created using the function palette.
  library(ggplot2)
  library(survminer)
  library(survival)
  library(ggthemes)
  library(ggplotify)
  library(ggsci)
  
  colnames(dat) <- c("OS.time", "OS", "Groups")
  # 将时间转换成年份
  if (max(dat[, 1]) > 365) {
    dat[, 1] <- dat[, 1] / 365
  } else if (max(dat[, 1]) > 24) {
    dat[, 1] <- dat[, 1] / 12
  } else {
    dat <- dat
  }
  fit <- survfit(Surv(OS.time, OS) ~ Groups,data=dat)
  surv.fit <- ggsurvplot(fit, data = dat, palette = col,
                         pval = TRUE, 
                         pval.method = T,
                         pval.method.size = 4,
                         pval.method.coord = c(0, 0.15),
                         surv.median.line='hv',
                         linetype = 1, 
                         pval.coord=c(0, 0.05), 
                         pval.size = 4,
                         risk.table = risk.table,
                         risk.table.y.text = FALSE,
                         legend.title = title,
                         legend.labs = lables,
                         xlab = 'Time(years)',
                         ggtheme=theme_bw(),
                         tables.height = tables.height)
  # 将图形转换为 ggplot 对象
  if (risk.table) {
    surv.fit1 <- surv.fit$plot + 
      theme(legend.position=c(1,1), 
            legend.justification=c(1,1),
            plot.margin=unit(c(0.1, 0.15, 0, 0.15), "inches"),
            legend.background = element_rect(fill = NA, colour = NA)
            # ,
            # axis.text.x=element_blank(),
            # axis.title.x=element_blank()
      )
    
    surv.fit2 <- surv.fit$table + 
      theme(plot.title=element_blank(),
            plot.margin=unit(c(0, 0.15, 0, 0.15), "inches")) +
      ylab('')
    surv.fit <- ggpubr::ggarrange(surv.fit1,
                                  surv.fit2, 
                                  ncol = 1, 
                                  nrow = 2,
                                  heights = c(1 - tables.height, 
                                              tables.height),
                                  align = "hv")
  } else {
    surv.fit <- surv.fit$plot + 
      theme(legend.position=c(1,1), 
            legend.justification=c(1,1),
            # plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
            legend.background = element_rect(fill = NA, colour = NA))
  }
  return(surv.fit)
}

#
module.name='Enet[alpha=0.1]'
module.risk <- read.delim('E:/99other/31HNSCC/4Development.CLMS/1model/Results/riskscore_mat.txt',sep='\t',header = T,check.names = F,row.names = 1)
module.gene <- read.delim('E:/99other/31HNSCC/4Development.CLMS/1model/Results/fea_df.txt',sep='\t',header = T,check.names = F) %>% 
  filter(algorithm==module.name) %>% pull(features)
module.gene %>% write.table("results/module.gene.txt",quote = F,sep='\t',row.names = F,col.names = F)
##  Fig.4B-module.gene
#多因素cox分析重新计算
calcoef <- function(cli,data,module.gene){
  selectGenes <- intersect(module.gene,row.names(data))
  data <- data[selectGenes,rownames(cli)] 
  cli  <- cbind(cli,t(data))
  fmla <- paste0("Surv(OS.time, OS)~" ,paste0(selectGenes,collapse = '+')) %>% as.formula()
  cox  <- coxph(fmla, data = cli)
  coef <- coef(cox)
  return(coef)
}

df <- calcoef(tcga.cli,tcga.data,module.gene) %>% 
  data.frame(coefficients=.) %>% 
  rownames_to_column("module.gene") %>% 
  arrange(coefficients)
df$module.gene <- factor(df$module.gene,levels = df$module.gene)

write.table(df,file="Results/calcoef.txt", sep="\t", quote=F, row.names = F)



p <- df %>% ggplot() + theme_bw() + 
  geom_segment(color="grey",linewidth=1.5,
               aes(x = 0, y = module.gene,xend = coefficients, yend = module.gene)
  ) +
  geom_point(aes(x=coefficients,y=module.gene,color=module.gene),size=6) +
  theme(
    legend.position = 'none',
    panel.grid.major = element_line(linetype = "dashed"),
    panel.grid.minor = element_blank()
  ) + labs(x="",y="") #+ ggsci::scale_color_simpsons()
ggsave(plot=p,filename = "results/Fig.4B.pdf",width = 8,height = 5)

##################################### Fig.4C-univariate Cox regression analysis
meta.cox <- read.table("results/meta.cox.txt") %>% rownames_to_column("gene") %>% mutate(cohorts="META")
tcga.cox <- read.table("results/tcga.cox.txt") %>% rownames_to_column("gene") %>% mutate(cohorts="TCGA")
df <- rbind(tcga.cox,meta.cox) %>% filter(gene %in% module.gene)

p <- df %>% ggplot() + theme_void() +  
  geom_segment(aes(y=gene,yend=gene,x=HR.95L,xend=HR.95H)) +
  geom_vline(xintercept =1, linetype = "dashed") + 
  geom_point(aes(x=HR,y=gene,fill=pvalue,color=pvalue),shape=23,size=2) + 
  facet_wrap(~cohorts,strip.position = "bottom") +
  scale_color_continuous(type = "viridis") +
  scale_fill_continuous(type = "viridis")  +
  theme(axis.text.y = element_text(),
        strip.text = element_text(size = 15))
ggsave(plot=p,filename = "results/Fig.4C.pdf",height = 5,width = 6)
##################################### Fig.4D-G Survival analysis with high CMLS and low CMLS
tcga.cli <- read.table("E:/99other/31HNSCC/1data/TCGA/results/tcga.cli.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
meta.cli <- read.table("E:/99other/31HNSCC/1data/GEO/META/results/meta.cli.txt")

tcga.cli$RS <- module.risk[row.names(tcga.cli),module.name]
meta.cli$RS <- module.risk[row.names(meta.cli),module.name]

tcga.RS.cut <- surv_cutpoint(tcga.cli, time = "OS.time", event = "OS",variables = "RS")$cutpoint[1,1]
meta.RS.cut <- surv_cutpoint(meta.cli, time = "OS.time", event = "OS",variables = "RS")$cutpoint[1,1]

tcga.cli <- tcga.cli  %>% mutate(binRS = ifelse(RS > tcga.RS.cut,"High","Low"))
meta.cli <- meta.cli  %>% mutate(binRS = ifelse(RS > meta.RS.cut,"High","Low"))

tcga.cli %>% write.table("results/tcga.cli.txt",quote = F,sep='\t')
meta.cli %>% write.table("results/meta.cli.txt",quote = F,sep='\t')

p <- ggsurvplotKM(tcga.cli[, c("OS.time", "OS", "binRS")],
                  lables = c('High', 'Low'),col = c('#AA1701','#0048A1'),
                  title = 'TCGA set')
ggsave(plot=p,filename = "results/Fig.4D.pdf",height = 5,width = 6)

p <- ggsurvplotKM(meta.cli[, c("OS.time", "OS", "binRS")],
                  lables = c('High', 'Low'),col = c('#AA1701','#0048A1'),
                  title = 'META set')
ggsave(plot=p,filename = "results/Fig.4E.pdf",height = 5,width = 6)
##################################### Fig.S4-6 Survival analysis with high and low CMLS genes
dir.create('results/Fig.S4-6')

setwd("E:/99other/31HNSCC/4Development.CLMS")
tcga.cli.genes <- read.table("results/run.input/tcga_exp.txt.gz",check.names = F) %>% t() %>% scale()
tcga.cli.genes <- ifelse(tcga.cli.genes>0,"High","Low")
meta.cli.genes <- read.table("results/run.input/meta_exp.txt.gz",check.names = F) %>% t() %>% scale()
meta.cli.genes <- ifelse(meta.cli.genes>0,"High","Low")

tcga.cli <- read.table("E:/99other/31HNSCC/1data/TCGA/results/tcga.cli.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
meta.cli <- read.table("E:/99other/31HNSCC/1data/GEO/META/results/meta.cli.txt")
tcga.cli <- tcga.cli[rownames(tcga.cli.genes),] %>% cbind(tcga.cli.genes)
meta.cli <- meta.cli[rownames(meta.cli.genes),] %>% cbind(meta.cli.genes)

dir.create("DSS/results")
##DSS PFI
library(tidyverse)
library(data.table)

fread("DSS/PMC6066282-TCGA-CDR-clinical.txt") %>% 
  filter(type=="HNSC")  %>% 
  dplyr::select(bcr_patient_barcode,OS,OS.time,DSS,DSS.time,PFI,PFI.time) %>% 
  data.frame(row.names = 1) %>% 
  write.table("DSS/results/tcga.cli.others.txt",sep="\t",quote = F)


tcga.cli.others = read.table("DSS/results/tcga.cli.others.txt")
tcga.com.samples = intersect(row.names(tcga.cli),row.names(tcga.cli.others))
tcga.com.cli <- tcga.cli[tcga.com.samples,-c(1,2)] %>% cbind(tcga.cli.others[tcga.com.samples,])

# META OS module.gene
for (varialbe in module.gene){
  title = sprintf("%s-%s-%s","META","OS",varialbe)
  dat = meta.cli[,c("OS.time","OS",varialbe)]
  p <- ggsurvplotKM(dat,lables = c('High', 'Low'),col = c('#AA1701','#0048A1'),title = title)
  ggsave(plot=p,filename = sprintf("results/Fig.S4-6/%s.pdf",title),height = 7,width = 7)
}

# TCGA OS module.gene
for (varialbe in module.gene){
  title = sprintf("%s-%s-%s","TCGA","OS",varialbe)
  dat = tcga.com.cli[,c("OS.time","OS",varialbe)]
  p <- ggsurvplotKM(dat,lables = c('High', 'Low'),col = c('#AA1701','#0048A1'),title = title)
  ggsave(plot=p,filename = sprintf("results/Fig.S4-6/%s.pdf",title),height = 7,width = 7)
}

# TCGA DSS module.gene
for (varialbe in module.gene){
  title = sprintf("%s-%s-%s","TCGA","DSS",varialbe)
  dat = tcga.com.cli[,c("DSS.time","DSS",varialbe)]
  p <- ggsurvplotKM(dat,lables = c('High', 'Low'),col = c('#AA1701','#0048A1'),title = title)
  ggsave(plot=p,filename = sprintf("results/Fig.S4-6/%s.pdf",title),height = 7,width = 7)
}

# TCGA PFI module.gene
for (varialbe in module.gene){
  title = sprintf("%s-%s-%s","TCGA","PFI",varialbe)
  dat = tcga.com.cli[,c("PFI.time","PFI",varialbe)]
  p <- ggsurvplotKM(dat,lables = c('High', 'Low'),col = c('#AA1701','#0048A1'),title = title)
  ggsave(plot=p,filename = sprintf("results/Fig.S4-6/%s.pdf",title),height = 7,width = 7)
}

####5.1.Comparison.prognostic.signatures####

setwd("E:/99other/31HNSCC/5Comparison.prognostic.signatures")
dir.create('results')
# 
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(randomForestSRC)
  library(glmnet)
  library(plsRcox)
  library(superpc)
  library(gbm)
  library(CoxBoost)
  library(survivalsvm)
  library(tibble)
  library(BART)
  library(stringr)
  library(ComplexHeatmap)
})
#
tcga.data <- read.table("E:/99other/31HNSCC/1data/TCGA/results/tcga.mRNA.exp.txt.gz",check.names = F)
tcga.cli <- read.table("E:/99other/31HNSCC/4Development.CLMS/results/run.input/tcga.cli.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

meta.data = read.table("E:/99other/31HNSCC/1data/GEO/META/results/meta.combat.txt.gz",check.names = F)
meta.cli = read.table("E:/99other/31HNSCC/1data/GEO/META/results/meta.cli.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

signatures = data.table::fread('signatures.txt')

#
sig2df = function(sig){
  rs = list()
  
  sig1 = sig %>% filter(Gene %in% rownames(tcga.data))
  RS1 = colSums(tcga.data[sig1$Gene,] * sig1$Coef) %>% data.frame(RS=.)
  rs$TCGA <- tcga.cli %>% dplyr::select(OS.time,OS) %>% cbind(RS1)
  
  sig2 = sig %>% filter(Gene %in% rownames(meta.data))
  RS2 = colSums(meta.data[sig2$Gene,] * sig2$Coef) %>% data.frame(RS=.)
  rs$META <- meta.cli %>% dplyr::select(OS.time,OS) %>% cbind(RS2)
  rs2df <- function(rs,modelName){
    cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])})) %>% 
      rownames_to_column('ID')
    cc$Model  <- modelName
    return(cc)
  }
  cc = rs2df(rs,sig$Author[1])
  return(cc)
}
calcoef <- function(cli,data,module.gene){
  selectGenes <- intersect(module.gene,row.names(data))
  selectSamples <- intersect(rownames(cli),colnames(data))
  data <- data[selectGenes,selectSamples]
  cli <- cli[selectSamples,] %>% cbind(t(data))
  fmla <- paste0("Surv(OS.time, OS)~" ,paste0(selectGenes,collapse = '+')) %>% as.formula()
  cox  <- coxph(fmla, data = cli)
  coef <- as.data.frame(coef(cox))
  coef=data.frame(Gene=rownames(coef),Coef=coef[,1])
  return(coef)
}
#
#如果没有系数，就需要重新用多因素cox分析计算系数
head(signatures)
signatures1=signatures[,c(1,2,3,4)]
signatures_new = split(signatures1,signatures1$Author) %>%  
  lapply(function(x){calcoef(cli=tcga.cli,data=tcga.data,module.gene=x$Gene)}) %>% 
  bind_rows(.id = "Author") %>% left_join(signatures1,by = join_by(Author,Gene)) 
signatures_new=signatures_new[,c("ID","PMID","Author","Gene","Coef")]
signatures=signatures_new
signatures_new %>% write.table(.,'results/signatures.new.txt',sep='\t',quote =F,row.names = F)
#
module.name='Enet[alpha=0.1]'
cindex.CMLS <- read.table("E:/99other/31HNSCC/4Development.CLMS/1model/Results/cindex_mat.txt")[module.name,] %>% 
  pivot_longer(names_to = "ID",values_to = "Cindex",cols = 1:2) %>% 
  mutate(Model="Prognostic Model")
cindex.22signatures = split(signatures,signatures$Author) %>% lapply(sig2df) %>% Reduce(rbind,.)
cindex = rbind(cindex.CMLS,cindex.22signatures)
#
cindex.select = cindex %>% filter(ID=="TCGA") %>% arrange(Cindex)
cindex.select$Model <- factor(cindex.select$Model,levels = cindex.select$Model)
p <- cindex.select %>% mutate(current=Model=="Prognostic Model") %>%  ggplot() + theme_classic() +  
  geom_vline(xintercept = 0.6,linewidth=1.5,linetype = "dashed") +
  geom_segment(linewidth=1.5,linetype = "dashed",
               aes(x = 0, y = Model,xend = Cindex, yend = Model,color=current)
  ) + geom_point(aes(x=Cindex,y=Model,color=current),size=6) +
  labs(x="TCGA Cindex",y='') + scale_color_manual(values=c("grey", "red")) +
  theme(legend.position = 'none')
ggsave("results/Fig.5A.pdf",p,height = 4.5,width = 5)
# Fig.5B
cindex.select = cindex %>% filter(ID=="META") %>% arrange(Cindex)
cindex.select$Model <- factor(cindex.select$Model,levels = cindex.select$Model)
p <- cindex.select %>% mutate(current=Model=="Prognostic Model") %>%  ggplot() + theme_classic() +  
  geom_vline(xintercept = 0.6,linewidth=1.5,linetype = "dashed") +
  geom_segment(linewidth=1.5,linetype = "dashed",
               aes(x = 0, y = Model,xend = Cindex, yend = Model,color=current)
  ) + geom_point(aes(x=Cindex,y=Model,color=current),size=6) +
  labs(x="META Cindex",y='') + scale_color_manual(values=c("grey", "red")) +
  theme(legend.position = 'none')
ggsave("results/Fig.5B.pdf",p,height = 4.5,width = 5)

####6.1.prepare####

setwd("E:/99other/31HNSCC/6Application")


tcga.data <- read.table('E:/99other/31HNSCC/1data/TCGA/results/tcga.mRNA.exp.txt.gz',check.names = F)
tcga.riskscore = read.table('E:/99other/31HNSCC/4Development.CLMS/1model/Results/tcga.cli.txt',header=T, sep="\t", check.names=F,row.names = 1)

meta.data <- read.table('E:/99other/31HNSCC/1data/GEO/META/results/meta.combat.txt.gz')
meta.riskscore = read.table('E:/99other/31HNSCC/4Development.CLMS/1model/Results/meta.cli.txt',header=T, sep="\t", check.names=F,row.names = 1)

module.gene <- read.table('module.gene.txt')

tcga.data = tcga.data[module.gene[,1],]
meta.data = meta.data[module.gene[,1],]

colnames(tcga.riskscore)[c(1,2,10,11)] = c("fustat","futime","riskScore","risk")
colnames(meta.riskscore)[c(1,2,4,5)] = c("fustat","futime","riskScore","risk")

tcga.riskscore$futime[tcga.riskscore$futime<=0]=0.003
meta.riskscore$futime[meta.riskscore$futime<=0]=0.003
tcga.riskscore$futime = tcga.riskscore$futime/365
meta.riskscore$futime = meta.riskscore$futime/365

tcga.riskscore$risk = gsub("High","high",tcga.riskscore$risk)
tcga.riskscore$risk = gsub("Low","low",tcga.riskscore$risk)
meta.riskscore$risk = gsub("High","high",meta.riskscore$risk)
meta.riskscore$risk = gsub("Low","low",meta.riskscore$risk)

tcga = cbind(tcga.riskscore[c("futime","fustat")],t(tcga.data),tcga.riskscore[c("riskScore","risk")])
meta = cbind(meta.riskscore[c("futime","fustat")],t(meta.data),meta.riskscore[c("riskScore","risk")])


write.table(data.frame(ID=rownames(tcga),tcga), file="risk.tcga.txt", sep="\t", quote=F, col.names=T,row.names = F)
write.table(data.frame(ID=rownames(meta),meta), file="risk.meta.txt", sep="\t", quote=F, col.names=T,row.names = F)


####6.2.indep####

library(survival)      #引用包
setwd("E:/99other/31HNSCC/6Application")     #设置工作目录

############绘制森林图函数
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width=6.6, height=4)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #绘制右边的森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
  dev.off()
}
############绘制森林图函数

############独立预后分析函数
indep=function(riskFile=null, cliFile=null, project=null){
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #读取临床文件
  
  #数据合并
  sameSample=intersect(row.names(cli),row.names(risk))
  risk=risk[sameSample,]
  cli=cli[sameSample,]
  rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
  
  #单因素独立预后分析
  uniCoxFile=paste0(project,".uniCox.txt")
  uniCoxPdf=paste0(project,".uniCox.pdf")
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniCoxFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniCoxFile, forestFile=uniCoxPdf, forestCol="#0048A1")
  
  
  #多因素独立预后分析
  multiCoxFile=paste0(project,".multiCox.txt")
  multiCoxPdf=paste0(project,".multiCox.pdf")
  uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
  rt1=rt[,c("futime","fustat",as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab, file=multiCoxFile, sep="\t", row.names=F, quote=F)
  bioForest(coxFile=multiCoxFile, forestFile=multiCoxPdf, forestCol="#AA1701")
}
############独立预后分析函数

#独立预后分析
indep(riskFile="risk.tcga.txt", cliFile="clinical1.txt", project="all")

####6.3.C-index####

#引用包
library(dplyr)
library(survival)
library(rms)
library(pec)

riskFile="risk.tcga.txt"     #风险文件
cliFile="clinical2.txt"      #临床数据文件
setwd("E:/99other/31HNSCC/6Application")   #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#定义颜色
bioCol=c("#AA1701","#0048A1","MediumSeaGreen","NavyBlue","Firebrick3","#8B668B","#FF4500")
colnames(rt)
#C-index值计算
riskScore=cph(Surv(futime,fustat)~riskScore, data=rt, surv=TRUE)
Age=cph(Surv(futime,fustat)~Age, data=rt, surv=TRUE)
Gender=cph(Surv(futime,fustat)~Gender, data=rt, surv=TRUE)
Stage=cph(Surv(futime,fustat)~Stage, data=rt, surv=TRUE)
Grade=cph(Surv(futime,fustat)~Grade, data=rt, surv=TRUE)
c_index  <- cindex(list("Risk score"=riskScore, 
                        "Age"=Age,
                        "Gender"=Gender,
                        "Stage"=Stage,
                        "Grade"=Grade),
                   formula=Surv(futime,fustat)~ .,
                   data=rt,
                   eval.times=seq(0,20,1),
                   splitMethod="bootcv",
                   B=1000
)
#输出图形
pdf(file="C-index.pdf", width=6, height=6)
plot(c_index, xlim=c(0,16), ylim=c(0.4,0.8), col=bioCol, legend.x=6, legend.y=0.82, legend.cex=1)
dev.off()

####6.4.Nomo####


#引用包
library(survival)
library(regplot)
library(rms)

riskFile="risk.tcga.txt"      #风险输入文件
cliFile="clinical4.txt"       #临床数据文件
setwd("E:/99other/31HNSCC/6Application")   #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("futime", "fustat", "risk")], cli)

#绘制列线图
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
             plots = c("density", "boxes"),
             clickable=F,
             title="",
             points=TRUE,
             droplines=TRUE,
             observation=rt[50,],
             rank="sd",
             failtime = c(1,3,5),
             prfail = F)

#列线图风险打分
nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

#校准曲线
pdf(file="calibration.pdf", width=6, height=6)
#1年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()


####6.5.riskDiff####


#引用包
library(limma)
expFile="tcga.mRNA.TPM.txt"  #表达数据文件
riskFile="risk.tcga.txt"       #风险文件
logFCfilter=1 #logFC过滤条件
fdrFilter=0.05#fdr过滤条件
setwd("E:/99other/31HNSCC/6Application") #设置工作目录

#读取表达数据文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
# rt=as.matrix(rt)
# rownames(rt)=rt[,1]
# exp=rt[,2:ncol(rt)]
# dimnames=list(rownames(exp), colnames(exp))
# data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
# data=avereps(data)
data = as.matrix(rt)
##去除正常样品
# group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2", "1", group)
# data=data[,group==0]
# data=t(data)
# rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
# data=avereps(data)
# data=t(data)

#读取risk文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(colnames(data), row.names(risk))
data=data[,sameSample]
risk=risk[sameSample,]

#提取low risk和high risk样品
riskLow=risk[risk$risk=="low",]
riskHigh=risk[risk$risk=="high",]
dataLow=data[,row.names(riskLow)]
dataHigh=data[,row.names(riskHigh)]
data=cbind(dataLow,dataHigh)
data=data[rowMeans(data)>0.5,]
conNum=ncol(dataLow)
treatNum=ncol(dataHigh)
Type=c(rep(1,conNum), rep(2,treatNum))

#差异分析
outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff, file="riskDiff.txt", sep="\t", row.names=F, quote=F)



####6.6.GO####


#引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05       #p值过滤条件
qvalueFilter=1       #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")

setwd("E:/99other/31HNSCC/6Application")      #设置工作目录
rt=read.table("riskDiff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#基因名字转换为基因id

rt = subset(rt,logFC>0)

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]#去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存富集结果
write.table(GO, file="GO.up.txt", sep="\t", quote=F, row.names = F)

#定义显示GO的数目
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#柱状图
pdf(file="GObarplot.up.pdf", width=10, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="GObubble.up.pdf", width=10, height=7)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

###########绘制GO圈图
data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

#整理圈图数据
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#绘制圈图主体部分
pdf("GO.circlize.up.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
#circos.par(track.min=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#绘制圈图中间的图例
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#绘制GO分类的图例
main.legend = Legend(
  labels = c("Biological Process", "Molecular Function","Cellular Component"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#绘制pvalue的图例
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()


#引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05       #p值过滤条件
qvalueFilter=1       #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")

setwd("E:/99other/31HNSCC/6Application")      #设置工作目录
rt=read.table("riskDiff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#基因名字转换为基因id

rt = subset(rt,logFC<0)

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]#去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存富集结果
write.table(GO, file="GO.down.txt", sep="\t", quote=F, row.names = F)

#定义显示GO的数目
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#柱状图
pdf(file="GObarplot.down.pdf", width=10, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="GObubble.down.pdf", width=10, height=7)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

###########绘制GO圈图
data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

#整理圈图数据
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#绘制圈图主体部分
pdf("GO.circlize.down.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
#circos.par(track.min=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#绘制圈图中间的图例
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#绘制GO分类的图例
main.legend = Legend(
  labels = c("Biological Process", "Molecular Function","Cellular Component"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#绘制pvalue的图例
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()


####6.7.KEGG####

#引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05      #p值过滤条件
qvalueFilter=1      #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

setwd("E:/99other/31HNSCC/6Application")  #设置工作目录
rt=read.table("riskDiff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#提取差异基因的名称,将基因名字转换为基因id

rt = subset(rt,logFC>0)

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]#去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#KEGG富集分析
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(KEGG, file="KEGG.up.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="KEGGbarplot.up.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

#气泡图
pdf(file="KEGGbubble.up.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()

###########绘制KEGG圈图
Pathway.col=c("#90EE90", "#E7B800", "#00AFBB")
showNum=18
data=KEGG[order(KEGG$p.adjust),]
if(nrow(KEGG)>showNum){
  data=data[1:showNum,]
}
data$Pathway="KEGG"
main.col = Pathway.col[as.numeric(as.factor(data$Pathway))]

#整理圈图数据
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(KEGG=data$ID,start=1,end=max(BgGene))
rownames(df) = df$KEGG
bed2 = data.frame(KEGG=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(KEGG=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(KEGG=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#绘制圈图主体部分
pdf(file="KEGG.circlize.up.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
#circos.par(track.min=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#绘制圈图中间的图例
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',Pathway.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#绘制KEGG分类的图例
main.legend = Legend(
  labels = c("KEGG"),  type="points",pch=15,
  legend_gp = gpar(col=Pathway.col), title_position = "topcenter",
  title = "Pathway", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#绘制pvalue的图例
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()



#引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05      #p值过滤条件
qvalueFilter=1      #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

setwd("E:/99other/31HNSCC/6Application")  #设置工作目录
rt=read.table("riskDiff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#提取差异基因的名称,将基因名字转换为基因id

rt = subset(rt,logFC<0)

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]#去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#KEGG富集分析
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(KEGG, file="KEGG.down.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="KEGGbarplot.down.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

#气泡图
pdf(file="KEGGbubble.down.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()

###########绘制KEGG圈图
Pathway.col=c("#90EE90", "#E7B800", "#00AFBB")
showNum=18
data=KEGG[order(KEGG$p.adjust),]
if(nrow(KEGG)>showNum){
  data=data[1:showNum,]
}
data$Pathway="KEGG"
main.col = Pathway.col[as.numeric(as.factor(data$Pathway))]

#整理圈图数据
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(KEGG=data$ID,start=1,end=max(BgGene))
rownames(df) = df$KEGG
bed2 = data.frame(KEGG=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(KEGG=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(KEGG=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#绘制圈图主体部分
pdf(file="KEGG.circlize.down.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
#circos.par(track.min=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#绘制圈图中间的图例
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',Pathway.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#绘制KEGG分类的图例
main.legend = Legend(
  labels = c("KEGG"),  type="points",pch=15,
  legend_gp = gpar(col=Pathway.col), title_position = "topcenter",
  title = "Pathway", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#绘制pvalue的图例
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()



####7.1.mutation####

#加载
library(maftools)       #引用包

#目录
setwd("E:/99other/31HNSCC/7Immune")
#读入
all.maf = read.maf("input_maftools.maf")

#读取风险文件
risk=read.table("risk.tcga.txt", header=T, sep="\t", check.names=F)
outTab=risk[,c(1, ncol(risk))]
colnames(outTab)=c("Tumor_Sample_Barcode", "Risk")
write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)

#分别获取高低的maf
low.tab = subset(x = outTab,Risk == "low")


low.maf = subsetMaf(maf = all.maf, tsb = low.tab[,1])
write.mafSummary(maf = low.maf,basename = "low.maf")

high.tab = subset(x = outTab,Risk == "high")
high.maf = subsetMaf(maf = all.maf, tsb = high.tab[,1])
write.mafSummary(maf = high.maf,basename = "high.maf")

#读取基因突变的文件
geneNum=20
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]

#注释的颜色
ann_colors=list()
col=c("#1E4A93","#D21E1F")
names(col)=c("low", "high")
ann_colors[["Risk"]]=col

#绘制低风险组瀑布图
pdf(file="low.pdf", width=6, height=6)
maf=read.maf(maf="low.maf_maftools.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

#绘制高风险组瀑布图
pdf(file="high.pdf", width=6, height=6)
maf=read.maf(maf="high.maf_maftools.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

####7.2.TIDE####


#准备TIDE输入文件
setwd("E:/99other/31HNSCC/7Immune")
pairFile="tcga.mRNA.TPM.txt"   

#读取表达文件，并对输入文件整理
rt=read.table(pairFile, header=T, sep="\t", check.names=F)
#rt=as.matrix(rt)

#rownames(rt)=rt[,1]
#exp=rt[,2:ncol(rt)]
exp=rt
#exp=log2(exp+1)
#使用apply函数
#exp = as.numeric(unlist(exp))

Expr <- t(apply(exp, 1, function(x)x-(mean(x)))) 
write.table(data.frame(ID=rownames(Expr),Expr,check.names = F),'tcga_normalize.txt', sep="\t", quote=F, row.names = TRUE)


#准备TIDE输入文件
setwd("E:/99other/31HNSCC/7Immune")
pairFile="tcga.mRNA.TPM.txt"  

#读取表达文件，并对输入文件整理
rt=read.table(pairFile, header=T, sep="\t", check.names=F)
#rt=as.matrix(rt)
#rownames(rt)=rt[,1]
#exp=rt[,2:ncol(rt)]
#exp=rt
exp=log2(exp+1)
#使用apply函数
#exp = as.numeric(unlist(exp))

Expr <- t(apply(exp, 1, function(x)x-(mean(x)))) 
write.table(data.frame(ID=rownames(Expr),Expr,check.names = F),'tcga_normalize_log.txt', sep="\t", quote=F, row.names = TRUE)


#引用包
library(limma)
library(ggpubr)
library(data.table)
tideFile="TIDE.log.csv"          #TIDE的打分文件
riskFile="risk.tcga.txt"      #风险文件
setwd("E:/99other/31HNSCC/7Immune")   #设置工作目录


#读取TIDE数据
tide=fread(tideFile)
tide = as.data.frame(tide)
rownames(tide) = tide[,1]
tide = tide[,"TIDE",drop = F]

# tide=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)
# group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2", "1", group)
# tide=tide[group==0,,drop=F]

# tide[,2] = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))

#dim(tide)
# tide1=tide[!duplicated(tide$V2),]
# 
# rownames(tide1) = tide1[,2]
# 
# tide = as.data.frame(tide1[,-2])
# rownames(tide) = rownames(tide1)
# colnames(tide) = "TIDE"

#tide = tide1[,-2]

#row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
# tide=avereps(tide)

#读取风险数据文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(tide, risk)

#设置比较组
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制小提琴图
gg1=ggviolin(data, x="risk", y="TIDE", fill = "risk", 
             xlab="", ylab="TIDE",
             palette=c("#AA1701","#0048A1"),
             legend.title="Risk",
             add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)

#输出图形	
pdf(file="TIDE.pdf", width=5, height=4.5)
print(gg1)
dev.off()


jco <- c("#AA1701","#0048A1")
gg1= ggplot(data = data,aes(x = risk, y = TIDE, fill = risk))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=2, position = position_jitterdodge(), color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("TIDE")) +
  xlab("")  +
  stat_compare_means(comparisons = my_comparisons) + 
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
        #legend.position = "none",
        #axis.title = element_blank()
  )

#输出图形	
pdf(file="TIDE1.pdf", width=5, height=4.5)
print(gg1)
dev.off()


#引用包
library(survival)
library(survminer)
library(data.table)
TIDEFile="TIDE.log.csv"            #肿瘤突变负荷文件
riskFile="risk.tcga.txt"      #风险文件
setwd("E:/99other/31HNSCC/7Immune")   #修改工作目录

#读取输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
#TIDE=read.table(TIDEFile, header=T, sep="\t", check.names=F, row.names=1)      #读取TIDE数据文件
#读取TIDE数据
tide=fread(TIDEFile)
tide = as.data.frame(tide)
rownames(tide) = tide[,1]
TIDE = tide[,"TIDE",drop = F]
# #合并数据
# TIDE[,2] = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(TIDE))
# #dim(TIDE)
# TIDE1=TIDE[!duplicated(TIDE$V2),]
# rownames(TIDE1) = TIDE1[,2]
# TIDE = as.data.frame(TIDE1[,-2])
# rownames(TIDE) = rownames(TIDE1)
# colnames(TIDE) = "TIDE"
#row.names(TIDE)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(TIDE))
sameSample=intersect(row.names(TIDE), row.names(risk))
TIDE=TIDE[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(risk, TIDE)

#获取肿瘤突变负荷最优cutoff
res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TIDE"))
cutoff=as.numeric(res.cut$cutpoint[1])
TIDEType=ifelse(data[,"TIDE"]<=cutoff, "L-TIDE", "H-TIDE")
scoreType=ifelse(data$risk=="low", "low risk", "high risk")
mergeType=paste0(TIDEType, "+", scoreType)

#定义生存分析函数
bioSurvival=function(surData=null, outFile=null){
  diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
  length=length(levels(factor(surData[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
  #print(surv_median(fit))
  
  #绘制生存曲线
  bioCol=c("#AA1701","#0048A1","NavyBlue","MediumSeaGreen","Firebrick3")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=surData,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="",
                     legend.labs=levels(factor(surData[,"group"])),
                     font.legend=10,
                     legend = c(0.8, 0.8),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=F,
                     cumevents=F,
                     risk.table.height=.25)
  #输出图形
  pdf(file=outFile, onefile = FALSE, width=5.5, height=4.8)
  print(surPlot)
  dev.off()
}

#绘制肿瘤突变负荷的生存曲线
data$group=TIDEType
bioSurvival(surData=data, outFile="TIDE.survival.pdf")

#绘制肿瘤突变负荷联合高低风险的生存曲线
data$group=mergeType
bioSurvival(surData=data, outFile="TIDE-risk.survival.pdf")


####7.3.TMBsurv####

#引用包
library(survival)
library(survminer)
tmbFile="TMB.txt"            #肿瘤突变负荷文件
riskFile="risk.tcga.txt"      #风险文件
setwd("E:/99other/31HNSCC/7Immune")    #修改工作目录

#读取输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)      #读取TMB数据文件

#合并数据
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(risk, tmb)

#获取肿瘤突变负荷最优cutoff
res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
scoreType=ifelse(data$risk=="low", "low risk", "high risk")
mergeType=paste0(tmbType, "+", scoreType)

#定义生存分析函数
bioSurvival=function(surData=null, outFile=null){
  diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
  length=length(levels(factor(surData[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
  #print(surv_median(fit))
  
  #绘制生存曲线
  bioCol=c("Firebrick2","DodgerBlue1","NavyBlue","MediumSeaGreen","Firebrick3")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=surData,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="",
                     legend.labs=levels(factor(surData[,"group"])),
                     font.legend=10,
                     legend = c(0.8, 0.8),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=F,
                     cumevents=F,
                     risk.table.height=.25)
  #输出图形
  pdf(file=outFile, onefile = FALSE, width=5.5, height=4.8)
  print(surPlot)
  dev.off()
}

#绘制肿瘤突变负荷的生存曲线
data$group=tmbType
bioSurvival(surData=data, outFile="TMB.survival.pdf")

#绘制肿瘤突变负荷联合高低风险的生存曲线
data$group=mergeType
bioSurvival(surData=data, outFile="TMB-risk.survival.pdf")

####7.4.immuneCor####


#引用包
library(limma)
library(scales)
library(ggplot2)
library(ggtext)
riskFile="risk.tcga.txt"      #风险输入文件
immFile="infiltration_estimation_for_tcga.csv"     #免疫细胞浸润文件
setwd("E:/99other/31HNSCC/7Immune")     #设置工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取免疫细胞浸润文件
immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
immune=avereps(immune)

#对风险文件和免疫细胞浸润文件取交集，得到交集样品
sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, "riskScore"]
immune=immune[sameSample,]

#对风险打分和免疫细胞进行相关性分析
x=as.numeric(risk)
outTab=data.frame()
for(i in colnames(immune)){
  y=as.numeric(immune[,i])
  corT=cor.test(x, y, method="spearman")
  cor=corT$estimate
  pvalue=corT$p.value
  if(pvalue<0.05){
    outTab=rbind(outTab,cbind(immune=i, cor, pvalue))
  }
}
#输出相关性结果
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)

#绘制气泡图
corResult=read.table("corResult.txt", head=T, sep="\t")
corResult$Software=sapply(strsplit(corResult[,1],"_"), '[', 2)
corResult$Software=factor(corResult$Software,level=as.character(unique(corResult$Software[rev(order(as.character(corResult$Software)))])))
b=corResult[order(corResult$Software),]
b$immune=factor(b$immune,levels=rev(as.character(b$immune)))
colslabels=rep(hue_pal()(length(levels(b$Software))),table(b$Software))     #定义颜色
pdf(file="cor.pdf", width=10, height=11)       #保存图片
ggplot(data=b, aes(x=cor, y=immune, color=Software))+
  labs(x="Correlation coefficient",y="Immune cell")+
  geom_point(size=4.1)+
  theme(panel.background=element_rect(fill="white",size=1,color="black"),
        panel.grid=element_line(color="grey75",size=0.5),
        axis.ticks = element_line(size=0.5),
        axis.text.y = ggtext::element_markdown(colour=rev(colslabels)))
dev.off()


####7.5.immFunction####


#引用包
library(limma)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)
expFile="tcga.mRNA.TPM.txt"             #表达输入文件
gmtFile="immune.gmt"             #免疫功能数据集文件
riskFile="risk.tcga.txt"              #风险文件
socreFile="immFunScore.txt"      #免疫功能打分的输出文件
setwd("E:/99other/31HNSCC/7Immune")      #设置工作目录

#读取表达输入文件，并对输入文件处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
# rt=as.matrix(rt)
# rownames(rt)=rt[,1]
# exp=rt[,2:ncol(rt)]
exp=as.matrix(rt)
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]

#读取数据集文件
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssgsea分析
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#定义ssGSEA score矫正函数
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA score进行矫正
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file=socreFile, sep="\t", quote=F, col.names=F)

#去除正常样品
# group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2", "1", group)
# data=t(data[,group==0])
# rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
# data=avereps(data)
data=t(data)
#读取风险文件
risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)

#合并数据
sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,"risk",drop=F]
rt1=cbind(data, risk)

#对免疫相关功能绘制箱线图
data=melt(rt1,id.vars=c("risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("low","high"))
p=ggboxplot(data, x="Type", y="Score", color = "Risk",
            ylab="Score",add = "none",xlab="",palette = c("#0048A1","#AA1701") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Risk),symnum.s=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")

#输出图片文件
pdf(file="immFunction.pdf", width=10, height=5)
print(p)
dev.off()


jco <- c("#AA1701","#0048A1")

boxplot=ggplot(data = data,aes(x = Type, y = Score, fill = Risk))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)+
  geom_point(shape = 21, size=0.5, position = position_jitterdodge(), color="black", alpha=0.05,)+
  theme_classic() +
  ylab(expression("Score")) +
  xlab("")  +
  
  rotate_x_text(50)+
  stat_compare_means(aes(group=Risk),method = "wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")+
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    #legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
        #legend.position = "none",
        #axis.title = element_blank()
  )


#输出图片
pdf(file="immFunction.boxplot1.pdf", width=11, height=5.5)
print(boxplot)
dev.off()


####7.6.checkpoint####


#引用包
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="tcga.mRNA.TPM.txt"      #表达输入文件
riskFile="risk.tcga.txt"       #风险输入文件
geneFile="ICGs.txt"       #免疫检查点的基因文件
setwd("E:/99other/31HNSCC/7Immune")    #设置工作目录

#读取基因表达文件,并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
# rt=as.matrix(rt)
# rownames(rt)=rt[,1]
# exp=rt[,2:ncol(rt)]
exp=as.matrix(rt)
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因文件
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])
data=log2(data+1)

# #删除正常样品
# group=sapply(strsplit(row.names(data),"\\-"),"[",4)
# group=sapply(strsplit(group,""),"[",1)
# group=gsub("2","1",group)
# data=data[group==0,]
# row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
# data=avereps(data)

#合并数据
risk=read.table(riskFile, sep="\t", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,],risk[sameSample,])
rt1=rt1[,c(sameGene,"risk")]

#提取显著差异的基因
sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
  if(sd(rt1[,i])<0.001){next}
  wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"risk"])
  pvalue=wilcoxTest$p.value
  if(wilcoxTest$p.value<0.05){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "risk")
rt1=rt1[,sigGene]

#把数据转换成ggplot2输入文件
rt1=melt(rt1,id.vars=c("risk"))
colnames(rt1)=c("risk","Gene","Expression")

#设置比较组
group=levels(factor(rt1$risk))
rt1$risk=factor(rt1$risk, levels=c("low","high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

#绘制箱线图
boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="risk",
                  xlab="",
                  ylab="Gene expression",
                  legend.title="Risk",
                  width=0.8,
                  palette = c("DodgerBlue1","Firebrick2") )+
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),
                     method="wilcox.test",
                     symnum.s=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")

#输出图片
pdf(file="checkpoint.diff.pdf", width=10, height=4.5)
print(boxplot)
dev.off()



jco <- c("#AA1701","#0048A1")

boxplot=ggplot(data = rt1,aes(x = Gene, y = Expression, fill = risk))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)+
  geom_point(shape = 21, size=0.5, position = position_jitterdodge(), color="black", alpha=0.05,)+
  theme_classic() +
  ylab(expression("Gene expression")) +
  xlab("")  +
  
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    #legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
        #legend.position = "none",
        #axis.title = element_blank()
  )

#输出图片
pdf(file="checkpoint.diff1.pdf", width=10, height=4.5)
print(boxplot)
dev.off()


####7.7.radiotherapy####

#引用包
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
library(reshape2)
library(ggpubr)
library(readxl)

#设置工作目录
setwd("E:/99other/31HNSCC/7Immune")
pFilter=1  
#读取表达输入文件,并对数据进行处理
rt = read.table("tcga.mRNA.TPM.txt", header=T, sep="\t", check.names=F)
# rt=as.matrix(rt)
# rownames(rt)=rt[,1]
# exp=rt[,2:ncol(rt)]
exp=as.matrix(rt)

dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)


#读入gmt文件
geneSets=getGmt("Radiotherapy Predicited Pathways.gmt", geneIdType=SymbolIdentifier())

#分析 method=c("gsva", "ssgsea", "zscore", "plage")
gsvaResult=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

#对打分进行标准化
normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
gsvaResult=normalize(gsvaResult)

#导出
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="ssgseaOut.radiotherapy.txt", sep="\t", quote=F, col.names=F)


#读取风险输入文件
risk=read.table("risk.TCGA.txt", header=T, sep="\t", check.names=F, row.names=1)

gsvaResult = t(gsvaResult)
#数据合并
sameSample=intersect(row.names(risk), row.names(gsvaResult))
risk=risk[sameSample, "risk",drop=F]
gsvaResult=gsvaResult[sameSample,,drop=F]
rt=cbind(risk, gsvaResult)

#设置比较组
rt$risk=factor(rt$risk, levels=c("low", "high"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#对药物进行循环, 绘制箱线图
for(drug in colnames(rt)[2:ncol(rt)]){
  rt1=rt[,c(drug, "risk")]
  colnames(rt1)=c("Radiotherapy", "Risk")
  rt1=na.omit(rt1)
  rt1$Radiotherapy=log2(rt1$Radiotherapy+1)
  #差异分析
  test=wilcox.test(Radiotherapy ~ Risk, data=rt1)
  diffPvalue=test$p.value
  #对满足条件的药物绘制箱线图
  if(diffPvalue<pFilter){
    boxplot=ggboxplot(rt1, x="Risk", y="Radiotherapy", fill="Risk",
                      xlab="Risk",
                      ylab=paste0(drug, ""),
                      legend.title="Risk",
                      palette=c("#1E4A93","#D21E1F")
    )+ 
      stat_compare_means(comparisons=my_comparisons)
    #输出图形
    pdf(file=paste0("Radiotherapy.", drug, ".pdf"), width=5, height=4.5)
    print(boxplot)
    dev.off()
  }
}


####7.8.oncoPredict####


#引用包
library(limma)
library(oncoPredict)
library(parallel)
set.seed(12345)

expFile="tcga.mRNA.TPM.txt"     #表达数据文件
setwd("E:/99other/31HNSCC/7Immune")     #设置工作目录

#读取表达输入文件,并对数据进行处理
rt = read.table(expFile, header=T, sep="\t", check.names=F)
# rt=as.matrix(rt)
# rownames(rt)=rt[,1]
# exp=rt[,2:ncol(rt)]
exp=as.matrix(rt)

dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
#data=data[rowMeans(data)>0.5,]
# colnames(data)=gsub("(.*?)\\_(.*?)", "\\2", colnames(data))

#读取药物敏感性文件
GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res) 

#药物敏感性
calcPhenotype(trainingExprData = GDSC2_Expr,    #train组的表达数据
              trainingPtype = GDSC2_Res,        #train组的药物数据
              testExprData = data,              #test组的表达数据
              batchCorrect = 'eb',              #批次矫正的方法
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,      #去除波动小的基因
              minNumSamples = 10,               #最小的样品数目
              printOutput = TRUE,               #是否输出结果
              removeLowVaringGenesFrom = 'rawData')

####7.9.boxplot####

#引用包
library(limma)
library(ggplot2)
library(ggpubr)

pFilter=0.05                      #pvalue过滤条件
riskFile="risk.TCGA.txt"            #风险文件
drugFile="DrugPredictions2.csv"     #药物敏感性文件
setwd("E:/99other/31HNSCC/7Immune/calcPhenotype_Output")     #设置工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取药物敏感性文件
senstivity=read.csv(drugFile, header=T, sep=",", check.names=F, row.names=1)
colnames(senstivity)=gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity))

#数据合并
sameSample=intersect(row.names(risk), row.names(senstivity))
risk=risk[sameSample, "risk",drop=F]
senstivity=senstivity[sameSample,,drop=F]
rt=cbind(risk, senstivity)

#设置比较组
rt$risk=factor(rt$risk, levels=c("low", "high"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#对药物进行循环, 绘制箱线图
for(drug in colnames(rt)[2:ncol(rt)]){
  rt1=rt[,c(drug, "risk")]
  colnames(rt1)=c("Drug", "Risk")
  rt1=na.omit(rt1)
  rt1$Drug=log2(rt1$Drug+1)
  #差异分析
  test=wilcox.test(Drug ~ Risk, data=rt1)
  diffPvalue=test$p.value
  #对满足条件的药物绘制箱线图
  if(diffPvalue<pFilter){
    boxplot=ggboxplot(rt1, x="Risk", y="Drug", fill="Risk",
                      xlab="Risk",
                      ylab=paste0(drug, " senstivity"),
                      legend.title="Risk",
                      palette=c("#1E4A93","#D21E1F")
    )+ 
      stat_compare_means(comparisons=my_comparisons)
    #输出图形
    pdf(file=paste0("drugSenstivity.", drug, ".pdf"), width=5, height=4.5)
    print(boxplot)
    dev.off()
  }
}







# attached packages:
# [1] sva_3.50.0             BiocParallel_1.36.0    genefilter_1.84.0      mgcv_1.9-1             factoextra_1.0.7      
# [6] IOBR_0.99.8            patchwork_1.2.0        tidyHeatmap_1.8.1      MOVICS_0.99.17         RTN_2.26.0            
# [11] umap_0.2.10.0          Rtsne_0.17             data.table_1.15.4      BART_2.9.9             nlme_3.1-164          
# [16] survivalsvm_0.0.5      CoxBoost_1.5           gbm_2.2.2              superpc_1.12           plsRcox_1.7.7         
# [21] glmnet_4.1-8           Matrix_1.6-5           randomForestSRC_3.3.0  lubridate_1.9.3        forcats_1.0.0         
# [26] stringr_1.5.1          purrr_1.0.2            readr_2.1.5            tidyr_1.3.1            tibble_3.2.1          
# [31] tidyverse_2.0.0        pec_2023.04.12         prodlim_2023.08.28     rms_6.8-0              Hmisc_5.1-3           
# [36] regplot_1.1            ComplexHeatmap_2.18.0  dplyr_1.1.4            RColorBrewer_1.1-3     circlize_0.4.16       
# [41] enrichplot_1.22.0      org.Hs.eg.db_3.18.0    clusterProfiler_4.10.1 maftools_2.18.0        survminer_0.4.9       
# [46] survival_3.5-8         ggtext_0.1.2           scales_1.3.0           readxl_1.4.3           reshape2_1.4.4        
# [51] pheatmap_1.0.12        GSVA_1.48.3            GSEABase_1.64.0        graph_1.80.0           annotate_1.80.0       
# [56] XML_3.99-0.16.1        AnnotationDbi_1.64.1   IRanges_2.36.0         S4Vectors_0.40.2       Biobase_2.62.0        
# [61] BiocGenerics_0.48.1    oncoPredict_1.2        ggpubr_0.6.0           ggplot2_3.5.0          limma_3.58.1   

