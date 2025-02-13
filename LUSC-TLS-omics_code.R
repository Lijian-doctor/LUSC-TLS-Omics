library(MOVICS)
library(dplyr)
#————————————————提取TLS相关信息——————————————————————————————————————————————————————————————————————
#————————mRNA.expr————————————————————————
# 选择您要筛选的行名
rows_to_select <- c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL11","CXCL13","CXCL13","CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1","CD4","CCR5","CXCR3","CSF2","IGSF6","IL2RA","CD38","CD40","CD5","MS4A1","SDC1","GFI1","IL1R1","IL1R2","IL10","CCL20","IRF4","TRAF6","STAT5A","TNFRSF17")

# 筛选特定的28行
fpkm2 <-fpkm
mRNA.expr <- fpkm2[rownames(fpkm2) %in% rows_to_select, ]

#—————————————————提取高变异信息——————————————————————————————————————————————————————————————

########做相关性，筛选出｜cor｜最大的前40个lncRNA，再导入
write.table(lncRNA.expr,"lncRNA.expr.txt",sep="\t")
write.table(mRNA.expr,"mRNA.expr.txt",sep="\t")
name <- read.table("lncRNA.expr_cor_result_top40.txt",sep="\t",header=T,row.names = 1)
# Get matches 
lncRNA.expr2 <- lncRNA.expr[rownames(lncRNA.expr) %in% 
                              rownames(name),]
####重命名
lncRNA.expr <-lncRNA.expr2

#—————————————————提取高变异信息————————————————————————————————————————————————————
meth.beta2 <-meth.beta

meth.beta3 <- meth.beta2[rownames(meth.beta2) %in% 
                           rownames(mRNA.expr),]
meth.beta <-meth.beta3
#—————————————————提取高变异信息———————————————————————————————————————————————————
mut.status2 <- mut.status
# Get matches 
mut.status3 <- mut.status2[rownames(mut.status2) %in% 
                             rownames(mRNA.expr),]
mut.status <-mut.status3

# 检查数据中是否包含缺失值
missing_values <- any(is.na(meth.beta))

# 如果数据中存在缺失值，可以选择删除包含缺失值的行或列，或进行其他处理
if (missing_values) {
  # 删除包含缺失值的行
  meth.beta1 <- meth.beta[complete.cases(meth.beta), ]
}

#——————————————————————————————————————————————————————————————————————————————————————
lusc.tcga <- list(mRNA.expr = mRNA.expr,
                  lncRNA.expr =lncRNA.expr ,
                  meth.beta=meth.beta,
                  mut.status=mut.status ,
                  count=count ,
                  fpkm=fpkm ,
                  maf=maf  ,
                  segment=segment  ,
                  clin.info=clin.info)

save(lusc.tcga, file = "2.LUSC.tcga.final.RData")  ##这部分是已经筛选了TLS相关基因作为分群数据的数据


#二、GET Module(获取模块)———————————————————————————————————————————————————————————————————————————————

load("2.LUSC.tcga.final.RData")#442个样本

#1）获取数据
names(lusc.tcga)

# extract multi-omics data
mo.data   <- lusc.tcga[1:4]

# extract raw count data for downstream analyses
count    <- lusc.tcga$count

# extract fpkm data for downstream analyses
fpkm     <- lusc.tcga$fpkm

# extract maf for downstream analysis
maf      <- lusc.tcga$maf

# extract segmented copy number for downstream analyses
segment   <- lusc.tcga$segment

# extract survival information
clin.info <- lusc.tcga$clin.info

#3）获得最佳聚类数

# identify optimal clustering number (may take a while)
dev.new()
optk.lusc <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,T), #是否是二进制，mo.data共有四个数据，其中第四个是二进制的
                         #（注：第4个数据是体细胞突变，它是一个二进制矩阵） note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.name    = "3.CLUSTER NUMBER OF TCGA-LUSC")

# perform iClusterBayes (may take a while)
load('2.LUSC.tcga.final.RData')
iClusterBayes.res <- 
  getiClusterBayes(data    = mo.data,
                   N.clust     = 4,
                   type        = c("gaussian","gaussian","gaussian","binomial"), #Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
                   #与矩阵列表相对应的数据类型，可以是高斯、二项式或可能性。
                   n.burnin    = 1800,  #一个整数值，用于指示MCMC burn-in次数。
                   #MCMC burn-in是一个概念，最早来自于统计学和机器学习中的马尔可夫链蒙特卡洛（MCMC）采样。
                   #在开始采样之前，通常会运行一段时间的马尔可夫链以达到平稳分布。这段时间被称为“burn-in”期。
                   #在RNN中，“burn-in”过程是指在开始学习之前先通过网络运行一段时间的序列数据，使得RNN的隐藏状态有机会收敛到某种有意义的状态。
                   #这样，当我们开始学习时，RNN的隐藏状态就已经包含了一些有用的历史信息，而不仅仅是初始的零状态。
                   #这种方法尤其对处理长序列的任务有帮助，因为这些任务可能需要网络记住距离当前时刻较远的历史信息。
                   n.draw      = 1200,    #一个整数值，用于指示MCMC绘制的数量。
                   # MCMC draw的数量是指从后验分布中抽取的样本数量。在MCMC采样中，我们使用马尔可夫链来生成后验分布的样本。
                   #MCMC采样的目的是生成足够多的样本，以便我们可以对后验分布进行准确的估计。
                   #在MCMC采样中，我们通常会运行多个独立的马尔可夫链，以确保我们获得的样本是从后验分布中独立抽取的。
                   #因此，MCMC draw的数量是指从每个独立的马尔可夫链中抽取的样本数量.
                   prior.gamma = c(0.5, 0.5, 0.5, 0.5),  #A numerical vector to indicate the prior probability for the indicator variable gamma of each subdataset.
                   #在统计学中，gamma分布是一种连续概率分布，通常用于建模正值的随机变量。
                   #在贝叶斯统计中，gamma分布通常用作参数的先验分布。在这种情况下，
                   #gamma分布的参数可以被视为超参数，因为它们控制了参数的分布。
                   #在某些情况下，我们可能会使用gamma分布作为指示变量的先验分布。
                   #指示变量是一个二元变量，它的值为0或1，通常用于表示某个事件是否发生。
                   #在这种情况下，gamma分布的参数被解释为指示变量为1的先验概率。
                   sdev        = 0.05,         #A numerical value to indicate the standard deviation of random walk proposal for the latent variable.
                   thin        = 3)            #为了减少自相关，使MCMC链变细的数值。

#> clustering done...
#> feature selection done...
#> 
#> 
#5）同时获取多种算法的结果

#如果同时在getMOIC（）中为methodslist参数指定一个算法列表，它将自动逐个执行具有默认参数的每个算法，
#并最终返回从指定算法派生的结果列表。既然iClusterBayes已经完成，让我们同时尝试其他9种默认参数算法。
#这需要一些时间，所以休息一下喝杯咖啡。

# perform multi-omics integrative clustering with the rest of 9 algorithms
pdf("5.moic.res.list.pdf") #这里的N.clust也要改成和上面（亚型数量）一样的
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 4,
                         type        = c("gaussian", "gaussian", "gaussian", "binomial"))
dev.off()

#使用append（）将iClusterBayes.res作为list附加到moic.res.list，已经有9个结果
moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))
# save moic.res.list to local path
save(iClusterBayes.res,moic.res.list, file = "4-5.moic.res.list.rda")

#6) get consensus from different algorithms、
load("4-5.moic.res.list.rda")

cmoic.lusc <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "6.CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")

###保存数据
save(iClusterBayes.res,moic.res.list,cmoic.lusc, file = "4-6.cmoic.lusc.rda")

#※导出10种聚类方法得到的综合聚类的结果
clust <- cmoic.lusc$clust.res

load("4-6.cmoic.lusc.rda")
getSilhouette(sil      = cmoic.lusc$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "7.1.SILHOUETTE",
              height   = 5.5,
              width    = 5)

#7.2) get multi-omics heatmap based on clustering result
load("4-6.cmoic.lusc.rda")

# convert beta value to M value for stronger signal
indata <- mo.data
indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))

# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # F:no center for mutation
                     scaleFlag  = c(T,T,T,F)) #F: no scale for mutation

#必须提取特征：

feat   <- iClusterBayes.res$feat.res
feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA.expr"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)


# set color for each omics data
# if no color list specified all subheatmaps will be unified to green and red color pattern
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, meth.col, mut.col)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = cmoic.lusc$clust.res, # cluster results，如果此处改为：iClusterBayes.res$clust.res，则意为用单个iClusterBayes聚类的结果来作为分群数据
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "7.2COMPREHENSIVE HEATMAP OF CONSENSUSMOIC_无表型")


'''''
#当然，由于moic.res.list中总共存储了10个结果，您也可以选择其中的任何一个来创建热图。
#在这里，我选择了COCA，它也返回如下所示的样本树状图：
# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = moic.res.list$COCA$clust.res, # cluster results
             clust.dend    = moic.res.list$COCA$clust.dend, # show dendrogram for samples
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF COCA")
'''''

#现在回到cmoic.coad的共识结果，它集成了10个算法，这次还提供了样本的注释来生成热图。
#由于getMoHeatmap（）的核心函数基于ComplexHeatmap R包，因此在创建注释时，
#应始终使用circize:：colorRamp2（）函数为连续变量（例如，本例中的年龄）生成颜色映射函数。

# extract PAM50, pathologic stage and age for sample annotation
library(survival)
surv.info <- lusc.tcga$clin.info
write.csv(annCol,'annCol.csv')

annCol    <- surv.info[,c("PATH_M_STAGE", "PATH_N_STAGE","PATH_T_STAGE", "AGE","SEX"), drop = FALSE]

# generate corresponding colors for sample annotation
annColors <- list(AGE    = circlize::colorRamp2(breaks = c(min(annCol$AGE),
                                                           median(annCol$AGE),
                                                           max(annCol$AGE)), 
                                                colors = c("#0000AA", "#555555", "#AAAA00")),
                  PATH_T_STAGE = c("T1"    = "green",
                                   "T1A"    = "green",
                                   "T1B"    = "green",
                             "T2"    = "blue",
                             "T2A"    = "blue",
                             "T2B"    = "blue",
                             "TB"    = "blue",
                             "T3"    = "red",
                             "T4"    = "yellow", 
                  PATH_M_STAGE = c("M0"    = "green",
                              "M1"    = "blue",
                              "M1A"    = "blue",
                              "M1B"    = "blue",
                              "MX"    = "black"),
                  PATH_N_STAGE = c("N0"    = "green",
                              "N1"    = "blue",
                              "N2"    = "red",
                              "N3"    = "yellow",
                              "NX"    = "black"),
                  SEX = c("Male"    = "green",
                          "Female"    = "blue")
)
)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = cmoic.lusc$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F), # show no dendrogram for features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annRow        = annRow, # mark selected features
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "7.2COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")


#在生存分析中，“censor”和“event”是两个重要的概念。其中，“censor”表示未发生事件的研究对象，而“event”则表示发生了事件的研究对象。
#在生存分析中，对于每个研究对象，记录其是否发生了事件（如死亡、复发等），并在最后的统计分析中，将其结果用0或1来表示。其中，0表示该对象未发生事件，1表示该对象发生了事件。

#三、 COMP Module————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

#After identification of cancer subtypes, it is essential to further characterize each subtype by discovering their difference from multiple aspects. To this end, MOVICS provides commonly used downstream analyses in cancer subtyping researches for easily cohesion with results derived from GET Module. Now let us check them out.

# survival comparison
load("4-6.cmoic.lusc.rda")

# 将列转换为数值型格式 
clin.info$futime <- as.numeric(clin.info$futime)

surv.lusc <- compSurv(moic.res         = cmoic.lusc,
                      surv.info        = clin.info,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      fig.name         = "1.KAPLAN-MEIER CURVE OF CONSENSUSMOIC")
#> --a total of 643 samples are identified.
#> --removed missing values.
#> --leaving 642 observations.
#> --cut survival curve up to 10 years.
#> 
#> 
fpkm<-lusc.tcga$fpkm
ydata<-lusc.tcga$clin.info

#————整理数据——————————————————————————————————————————————————————————————————


xdata.NoneScale.raw_t <- read.csv("fpkm_final.csv",header = T,row.names = 1)
xdata.NoneScale.raw <- xdata.NoneScale.raw_t %>% t()%>% as.data.frame()
ydata.raw <- read.csv("surv.info.csv",header = T,row.names = 1)


xdata.Scale.raw <- xdata.NoneScale.raw %>% 
  { (apply(., 2, sd) != 0) } %>% 
  { xdata.NoneScale.raw[, .] } %>% 
  scale %>% as.data.frame()

#完全整理好数据之后运行这行代码，保存数据为.RData
save(xdata.NoneScale.raw, xdata.Scale.raw, ydata.raw, file = "4.lusc.TCGA.mtx_TLS.3.19.RData")


######################5分群分析#####################
#5分群分析

#免疫浸润：这里用的是过滤过的数据，有19680个基因，442个样本
#fpkm <- read.csv("4-免疫浸润&代谢ssgsea.csv",row.names = 1,header = T)

#基因集需要是list为对象。
#基因的表达量需要是矩阵，行为基因，列为样本
#默认情况下，kcdf=”Gaussian”，适用于输入表达式值连续的情况，如对数尺度的微阵列荧光单元、RNA-seq log-CPMs、log-RPKMs或log-TPMs。当输入表达式值是整数计数时，比如那些从RNA-seq实验中得到的值，那么这个参数应该设置为kcdf=”Poisson”

library(GSEABase)

library(GSVA)

library(tidyverse)
library(dplyr)
######
#1. 获取geneSets 基因背景集
#load(file = "C00基因集（三个）.symbols.Rdata")##(处理好的gmt-rdata数据)

#2. 读取基因的表达矩阵
#rm(list=ls())
## 载入数据-如果数据里有NA，需要先删掉
deg_counts <- read.csv("4-免疫浸润&代谢ssgsea.csv", sep=",", header= T)#(自己的差异基因-基因名+所有的组——fpkm)
######导入数据data1之后，把第一列的名字自动加序号并变为行名
row.names(deg_counts)<-make.names(deg_counts[,1],TRUE)
deg_counts<-deg_counts[,-1] #删除第一列

gene_level_expression_t <-deg_counts%>% t()%>%as.data.frame()

#####log
logTPM <- log(gene_level_expression_t)
##logTPM转置并变为matrix
logTPM <-logTPM%>% t()%>%as.data.frame()%>%as.matrix()

##l变为list
l<- immunecell%>% as.list()
#3. 开始进行ssGSEA
ssgsea<- gsva(logTPM, l, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)#可以改为“gsva”, “ssgsea”, “zscore”, “plage”算法

# Min-Max标准化是指对原始数据进行线性变换，将值映射到[0，1]之间
# 这里是将每个样本中不同的免疫细胞比例标准化到0-1之间
ssgsea.1 <- ssgsea
for (i in colnames(ssgsea)) {
  #i <- colnames(ssgsea)[1]
  ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
  
}
apply(ssgsea.1[,1:6], 2, range)
##保存结果
write.csv(ssgsea,file="TCGA_lusc_4_免疫浸润_ssgsea（原始评分）.csv")         #####原始评分结果
write.csv(ssgsea.1,file="TCGA_lusc_4_免疫浸润_ssgsea（标准化）.csv")   ####标准化的结果


#ssgsea-代谢
deg_counts <- read.csv("4-免疫浸润&代谢ssgsea.csv", sep=",", header= T)#(自己的差异基因-基因名+所有的组——fpkm)
######导入数据data1之后，把第一列的名字自动加序号并变为行名
row.names(deg_counts)<-make.names(deg_counts[,1],TRUE)
deg_counts<-deg_counts[,-1] #删除第一列

gene_level_expression_t <-deg_counts%>% t()%>%as.data.frame()

#####log
logTPM <- log(gene_level_expression_t)
##logTPM转置并变为matrix
logTPM <-logTPM%>% t()%>%as.data.frame()%>%as.matrix()

##l变为list
l<- KEGG_metabolism%>% as.list()
#3. 开始进行ssGSEA
ssgsea<- gsva(logTPM, l, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)#可以改为“gsva”, “ssgsea”, “zscore”, “plage”算法

# Min-Max标准化是指对原始数据进行线性变换，将值映射到[0，1]之间
# 这里是将每个样本中不同的免疫细胞比例标准化到0-1之间
ssgsea.1 <- ssgsea
for (i in colnames(ssgsea)) {
  #i <- colnames(ssgsea)[1]
  ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
  
}
apply(ssgsea.1[,1:6], 2, range)
##保存结果
write.csv(ssgsea,file="TCGA_lusc_4_kegg代谢_ssgsea（原始评分）.csv")         #####原始评分结果
write.csv(ssgsea.1,file="TCGA_lusc_4_kegg代谢_ssgsea（标准化）.csv")   ####标准化的结果


library(ggplot2)
library(ggrepel)
#cluster 1
deg<-read.table("TopGenes_differential_analysis_4.0.tsv",header=T,sep="\t")#读入差异分析结果文件

deg$diff<-"no"#增加一列diff并全部用no填充
colnames(deg)
# [1] "gene_id"     "TS32_tpm"    "TS34_tpm"    "fc"          "log2fc"     
# [6] "pvalue"      "padjust"     "significant" "regulate"    "diff"  
deg$diff[(deg$log2FC>1.5) & (deg$padjust<0.01)]<-"up" #将logFC>1且FDR<0.05的行标注为up
deg$diff[(deg$log2FC<(-1.5)) & (deg$padjust<0.01)]<-"down" #将logFC小于-1且FDR<0.05的行标注为up

as.data.frame(table(deg$diff)) #统计上下调基因数目

deg$diff<-factor(deg$diff,
                 levels = c('up','no','down'),
                 labels = c("Up (96)","No (17698)","Down (1874)")) #将顺序设定为'up','no','down'，同时根据上面的结果标注DEGs的数目

#deg$name<-deg$gene_id #新增一列name，并且把gene的内容复制给name

#deg$name[(deg$log2FoldChange>(-4.5))&(deg$log2FoldChange<(5.5))]<- NA #筛选logFC小于(-2.5)和logFC大于(2.7)的基因，将不在该范围内的name值用NA填充




# cmp <- deg %>% select(name)
# write.csv(cmp,"1.csv")

ggplot(deg)+
  geom_point(aes(x=log2FC,y= -1*log10(padjust),color=diff)) +#画点，用diff进行填充颜色
  geom_vline(xintercept = c(1,-1),linetype="dashed",color="grey")+ #X轴辅助线
  geom_hline(yintercept =-log10(0.5),linetype="dashed",color="grey")+ #Y轴辅助线
  scale_color_manual(values=c("#ff8216","#EEE5D9","#018a89"))+ #将颜色改为"red","grey","green"，对应上面的'up','no','down'
  guides(color=guide_legend(title = " cluster 4|Log2FC|>1.5 & padj<0.01",override.aes = list(size=2)))+ #更改图例，将比较组，差异基因筛选标准和数目标注到图中
  scale_y_continuous(expand = c(0,0))+#坐标轴相交于原点
  theme_bw()+#主题
  theme(legend.position="top")+
  xlab("Log2FC")+
  ylab("-Log10(padj)")+
  # geom_text_repel(aes(x= -1*log10(padjust),y=log2fc,label=name),
  #                 max.overlaps = 10000,                    # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之。
  #                 size=3, # 字体大小box.padding=unit(0.5,'lines'),           # 标记的边距
  #                 point.padding=unit(0.1, 'lines'), 
  #                 segment.color='black',                   # 标记线条的颜色
  #                 show.legend=FALSE)+ #使用ggrepel包进行目标基因标注，可防止相近的点重叠到一起
  theme(panel.grid=element_blank())#去掉次要线条

ggsave(" cluster_4_火山图.pdf",width=4,height =5)
write.csv(deg,"cluster4_DEGs_up.down.csv")


##取5群差异基因的并集&GEPIA2中正常和肿瘤的差异基因的交集：共470个-画韦恩图
data<- read.table("venn.txt",header = TRUE,sep="\t")
head(data)

venn.plot <- venn.diagram(
  x = list(
    deg_our=data$cluster_top, normal_tumor=data$normal_tumor
  ),
  filename = NULL,
  col = "black",#边框颜色
  lwd = 2, # 边框线的宽度
  fill = c("#f5cac3","#588b8b"),
  alpha = 0.70,
  label.col = "black",
  cex = 2.0,
  cat.col = c("#f5cac3","#588b8b"),
  cat.cex = 1.5,#类名字体大小
  margin = 0.04, #边际距离
  scaled = F
)
pdf(file="venn.pdf",width = 5,height =5)
grid.draw(venn.plot)
dev.off()

#差异分析-----------------cs1&cs3----------#


#DESeq2进行差异分析
library(DESeq2)#调用包
data<-read.csv("1&3_des2.csv",header = T,row.names = 1)#读入count值表格，一般为6列，数值必须为整数
head(data)

group<-read.csv("1&3_des2_group.csv",header = T,row.names = 1)#读入分组文件
head(group)

dds <- DESeqDataSetFromMatrix(countData=round(data), 
                              colData=group, 
                              design=~group)
#5.1.指定因子水平（此处需把处理组往前放，对照组往后放）
dds$group <- factor(dds$group, levels = c("CS1","CS3"))
#差异分析，必要步骤
dds1<-DESeq(dds)

resultsNames(dds1) #查看结果
#dds<-DESeqDataSetFromMatrix(countData=data,colData = group,design = ~group)#获取dss文件，countData为count表格，colData是分组表格，design是差异比较组
dds
dds1<-DESeq(dds)#数据标准化，必要步骤
resultsNames(dds1)#查看结果
dds1$condition#默认后面比前面
res<-results(dds1)#结果转换，必要步骤
summary(res)#结果概要

write.csv(res,file="CS1&CS3_DESeq.csv")#导出差异分析结果

deg<-read.csv("CS1&CS3_DESeq.csv",header = T,row.names = 1)
deg$diff<-"no"#增加一列diff并全部用no填充
head(deg)

deg$diff[(deg$log2FoldChange>0.28) & (deg$pvalue<0.05)]<-"up" #将logFC>1且FDR<0.05的行标注为up
deg$diff[(deg$log2FoldChange<(-0.28)) & (deg$pvalue<0.05)]<-"down" #将logFC小于-1且FDR<0.05的行标注为up

as.data.frame(table(deg$diff)) #统计上下调基因数目
write.csv(deg,"CS1&CS3_DESeq.csv") 

#———step3:———————————cox单因素回归————————1766个——————————————————

#——————————log2fc.0.28&p<0.05————————————————
cox <- read.csv("c1&c3差异基因（log2fc.0.28&0.05）.csv",row.names = 1,header = T)
head(cox)

##计算
library(survival)
pFilter=0.05 #按p值保留结果：1是全部输出，0.05是0.05以下的输出
outResult=data.frame()
sigGenes=c("time","status")
for(i in colnames(cox[,3:ncol(cox)])){
  mycox <- coxph(Surv(time, status) ~ cox[,i], data = cox)
  mycoxSummary = summary(mycox)
  pvalue=mycoxSummary$coefficients[,"Pr(>|z|)"]
  if(pvalue<pFilter){
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=mycoxSummary$conf.int[,"exp(coef)"],
                          L95CI=mycoxSummary$conf.int[,"lower .95"],
                          H95CI=mycoxSummary$conf.int[,"upper .95"],
                          pvalue=mycoxSummary$coefficients[,"Pr(>|z|)"])
    )
    
  }
}

##写出结果文件
write.csv(outResult, file='c1&c3_outResult_0.05_单因素cox.csv')

#############单因素森林图
library(survival)
library(survminer)
library(coin)
library(ggpmisc)
library(dplyr)
library(forestplot)
outResult<-read.csv("c1&c3_outResult_0.05_单因素cox_森林图.csv", header = T, sep = ",")
#森林图绘制：
photo <- outResult %>%
  forestplot(labeltext = c(id,HR,pvalue), 
             clip = c(0.5, 3),#森林图中置信区间上下限设置
             xlog = F,
             zero=1.07)%>%
  fp_set_style(box = "#ff8216", #权重方块颜色
               line = "#ff8216") %>%#是否对X轴取对数
  fp_set_zebra_style("#f9f1dd") %>% #更改为斑马样式，并修改汇总菱形的填充颜色
  fp_add_header(id = c("","gene"),
                HR =c("","HR"),
                pvalue = c("","pvalue"))
pdf(file="cox_forestPlot.pdf",width = 6,height =5)
print(photo)
dev.off()

#——step4:lasso:建模————————用cox显著-CS1&CS3差异基因————————33个———————————————————————————————————————
load('4.lusc.TCGA.mtx_TLS.3.19.RData') 
small.subset.ll <- c("CPS1",
                     "CLTCL1",
                     "COCH",
                     "CCN6",
                     "DUSP4",
                     "CDO1",
                     "KREMEN2",
                     "CTSV",
                     "CASP5",
                     "AKR1E2",
                     "CA4",
                     "SLC49A3",
                     "PARM1",
                     "FGA",
                     "PRR15",
                     "LRRC55",
                     "NPIPB15",
                     "SULT1C2",
                     "SFMBT2",
                     "DMD",
                     "RNY1P16",
                     "RNU6.529P",
                     "SP5",
                     "IGKV6D.41",
                     "IGHD2.21",
                     "IGHV7.81",
                     "IGKV2D.30",
                     "UPK3B",
                     "IGHV3.37",
                     "MIR3142HG",
                     "IGKV2.40",
                     "H3C8",
                     "LOC102723996")


#————————————————————————————————————————————————————————————————————————————————————————————————
########确定好模型后，分别在Scale和Nonescale数据表中加入risk分群列，然后保存分群后的数据信息
# AKR1E2	CA4	CASP5	CCN6	CDO1	CLTCL1	CPS1	CTSV	DMD	DUSP4	FGA	H3C8	KREMEN2	LOC102723996	LRRC55	MIR3142HG	PARM1	PRR15	SLC49A3	SULT1C2
####Scale数据
##转换为数据框
xdata.tmp <- xdata.Scale.raw %>% as.data.frame()
##用得到的系数计算指数
index <- +xdata.tmp$AKR1E2*0.229198473259477+
  xdata.tmp$CA4*0.148312217553774+
  xdata.tmp$CASP5*0.00637410093834286 +
  xdata.tmp$CCN6 *0.173900788432551 -
  xdata.tmp$CDO1*0.11839886255889 -
  xdata.tmp$CLTCL1 *0.00694130784369265 +
  xdata.tmp$CPS1*0.0489633240878769 -
  xdata.tmp$CTSV  *0.043331201443875 -
  xdata.tmp$DMD *0.0440632233178641 +
  xdata.tmp$DUSP4*0.0616708981957713 +
  xdata.tmp$FGA *0.0614082707785454 +
  xdata.tmp$H3C8*0.07201920954824 +
  xdata.tmp$KREMEN2 *0.0490970135208903 +
  xdata.tmp$LOC102723996*0.0643937778088998 +
  xdata.tmp$LRRC55 *0.169797408085259+
  xdata.tmp$MIR3142HG*0.0203017920106524 -
  xdata.tmp$PARM1 *0.0181819616359712 -
  xdata.tmp$PRR15*0.00137821893635317 +
  xdata.tmp$SLC49A3 *0.0399429253527924 +
  xdata.tmp$SULT1C2*0.11429977464504 


xdata.tmp$index <- index
##根据index排序,并分为2群
xdata.Scale.risk<-xdata.tmp[order(xdata.tmp$index,decreasing = T),]
group <- c(rep("High-risk",219),rep("Low-risk",219))
cluster <- c(rep("1",219),rep("2",219))

xdata.Scale.risk$risk <- group
xdata.Scale.risk$cluster <- cluster

##转换为数据框
xdata.NoneScale.tmp <- xdata.NoneScale.raw %>% as.data.frame()
##把生成的指数这一列加到原来的表格中
xdata.NoneScale.tmp $index <- index
##根据index排序,并分为2群
xdata.NoneScale.risk<-xdata.NoneScale.tmp[order(xdata.NoneScale.tmp$index,decreasing = T),]
group <- c(rep("High-risk",219),rep("Low-risk",219))
cluster <- c(rep("1",219),rep("2",219))

xdata.NoneScale.risk$risk <- group
xdata.NoneScale.risk$cluster <- cluster

####ydata.risk数据
##重新对xdata和ydata的3个表格排序，使其sample 编号顺序一致，方便后续操作
xdata.NoneScale.risk <- xdata.NoneScale.risk[rownames(xdata.NoneScale.risk) %in% 
                                               rownames(ydata.raw),]
xdata.Scale.risk <- xdata.Scale.risk[rownames(xdata.Scale.risk) %in% 
                                       rownames(ydata.raw),]
ydata.raw_order  <- ydata.raw[rownames(xdata.Scale.risk), ]

##把分群的这一列信息加到ydata表格中，方便后续操作。
ydata.risk<-ydata.raw_order%>%
  as.data.frame()
ydata.risk$risk <- xdata.Scale.risk$risk
ydata.risk$index <- xdata.Scale.risk$index
ydata.risk$cluster <-xdata.Scale.risk$cluster



#———列线图———————————————————————————————————————————————————————————————————————————————————————————————
## 添加变量标签
surv.info <- lusc.tcga$clin.info
surv.info$pstage <- factor(surv.info$PATH_T_STAGE, levels = c("T1","T2","T3","T4","T4A","T4B"))
surv.info$SEX <- factor(surv.info$SEX, levels = c("Male","Female"))
surv.info$PATH_M_STAGE <- factor(surv.info$PATH_M_STAGE, levels = c("MX","M0","M1","M1A","M1B"))
surv.info$PATH_N_STAGE <- factor(surv.info$PATH_N_STAGE, levels = c("NX","N0","N1","N2","N3","N4"))
surv.info2<-surv.info
# Get matches 
surv.info2     <- surv.info2 [rownames(ydata.risk), ]

surv.info2$index<-ydata.risk$index
#data(surv.info2)
head(surv.info2)
## 根据nomogram要求处理数据
library(rms)
write.csv(surv.info2,'surv.info2.csv')
surv.info2 <- read.csv('surv.info2.csv',row.names = 1,header = T)
dd=datadist(surv.info2)
options(datadist="dd")


f2 <- psm(Surv(time,status) ~ index
          +AGE+PATH_T_STAGE+PATH_M_STAGE+PATH_N_STAGE, data =  surv.info2, dist='lognormal') 
med <- Quantile(f2) # 计算中位生存时间
surv <- Survival(f2) # 构建生存概率函数
## 绘制COX回归中位生存时间的Nomogram图
nom <- nomogram(f2, fun=function(x) med(lp=x),
                funlabel="Median Survival Time")
plot(nom)


## 绘制COX回归生存概率的Nomogram图
## 注意surv.info2数据的time是以’天‘为单位
nom <- nomogram(f2, fun=list(function(x) surv(365, x),
                             function(x) surv(1825, x),
                             function(x) surv(3650, x)),
                funlabel=c("1-year Survival Probability",
                           "5-year Survival Probability",
                           "10-year Survival Probability"))
plot(nom, xfrac=.6)


#----------------------------------------------------
## 评价COX回归的预测效果
## 计算c-index
rcorrcens(Surv(time,status) ~ predict(f2), data =  surv.info2)

# Somers' Rank Correlation for Censored Data    Response variable:Surv(time, status)
# 
#                 C   Dxy  aDxy    SD    Z P   n
# predict(f2) 0.601 0.202 0.202 0.048 4.22 0 438

#C-index，concordance index，一致性指数，主要用于计算生存分析中的COX模型预测值与真实之间的区分度，常用在评价患者预后模型的预测精度中。
#C-index在0.5-1之间（任意配对随机情况下一致与不一致刚好是0.5的概率）。0.5为完全不一致,说明该模型没有预测作用，1为完全一致，说明该模型预测结果与实际完全一致。一般情况下C-index在0.50-0.70为准确度较低：在0.71-0.90之间为准确度中等；而高于0.90则为高准确度。
#Calibration校准曲线
## 绘制校正曲线
## 重新调整模型函数f2，也即添加x=T, y=T
f2 <- psm(
  Surv(time,status) ~AGE+SEX+PATH_T_STAGE+PATH_M_STAGE+PATH_N_STAGE+index,
  data =  surv.info2, x=T, y=T, dist='lognormal') 

## 构建校正曲线
cal1 <- calibrate(f2, 
                  cmethod='KM', 
                  method="boot", 
                  u=365,# 需要与之前模型中定义好的time.inc一致，即365或730；
                  m=95, #每次抽样的样本量，可以更改抽样样本量来改善曲线
                  B=1000,
                  smoother="x") 
plot(cal1)

cal2 <- calibrate(f2, 
                  cmethod=c('hare', 'KM'),
                  method="boot", 
                  u=365, 
                  m=150,  ##m要根据样本量来确定，由于标准曲线一般将所有样本分为3组（在图中显示3个点）
                  B=40, 
                  bw=FALSE, 
                  rule="aic", 
                  type="residual", sls=0.05, aics=0, force=NULL,
                  estimates=TRUE,
                  pr=FALSE, what="observed-predicted", tol=1e-12, maxdim=5)
plot(cal2)

##绘制校正曲线
plot(cal1,lwd=2,lty=1,
     conf.int=T,# 是否显示置信区间
     errbar.col="#018a89",#直线曲线bar颜色
     col="#ff8216", # 曲线颜色
     xlim=c(0.6,0.95),ylim=c(0.5,1),
     xlab="Nomogram-Predicted Probability of 1-Year DFS",
     ylab="Actual 1-Year DFS (proportion)",
     subtitles = F)#不显示副标题


#——————————生存状态填充柱状图-原模型图——————————————————————

ydata.risk2<-ydata.risk %>%
  ggplot(aes(x=risk,fill=factor(status)))+
  geom_bar(stat="count",position='fill',width = .6) +
  geom_text(stat='count',aes(label=after_stat(count)), color="white", size=7,position=position_fill(0.5))+
  scale_fill_manual(values=c("#ff8216","#018a89"))+
  coord_flip()+
  theme(panel.background=element_rect(fill='transparent'),
        axis.line = element_line(arrow = arrow(length = unit(0.15, 'cm')),
                                 colour = "black"),panel.grid =element_blank(),
        text=element_text(size= 50),axis.text = element_text(colour="black"))
ggsave("生存状态填充柱状图.PDF",ydata.risk2,width=14,height=5.5)
write.csv(ydata.risk, file='生存状态填充柱状图.csv')

#----------------index&表型相关性-----------------#
##-------------------------相关性分析-----------------------------##
#install.packages("ggpmisc")
library(ggpmisc)

library(ggplot2)
cox  <-   read.csv("index&表型相关性.csv", header = T, sep = ",")

ggplot(cox, aes(x=index,y=Ragnum_Hypoxia_Score))+
  geom_point(aes(color="#ff8216",size=6))+
  geom_smooth(method = "lm",se = T)+
  stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(p.value.label),sep = "~~~")),formula = y~x,parse = TRUE, size=4.5)+
  stat_correlation(mapping = use_label("R"),method = "pearson",label.x = "right",size=4.5)+
  labs(x="Ragnum_Hypoxia_Score",y=NULL,size=18,face="bold")+
  scale_color_manual(values = c("#ff8216"))+
  scale_fill_manual(values = c("#ff8216"))+
  scale_y_continuous(expand = c(0, 0))+
  theme(
    axis.text.x = element_text(size = 12, hjust = 1),
    axis.title.y=element_text(size = 12,face="plain"), 
    plot.title = element_text(size=15,face="bold",hjust = 0.3),
    legend.position="none", #不需要图例
    #legend.position="bottom",  #图例在底下
    panel.background=element_rect(fill='transparent',color='black'),
    panel.grid =element_blank(),
    #axis.text = element_text(colour="black"),
    legend.key = element_blank(),
    legend.title=element_blank())#+
ggsave("Ragnum_Hypoxia_Score.pdf",width=4,height= 4)


#——————表型——————cox单因素回归——————————————————————————
cox <- read.csv("INDEX&表型单多因素cox.csv",row.names = 1,header = T)
head(cox)
# ##筛选出相关基因的列
# cox.select <- cox%>% 
#   as.data.frame() %>% 
#   dplyr::select(status,time,
#     AGE,
#     SEX,
#     M_STAGE,
#     N_STAGE,
    # T_STAGE,
    # index,
    # cluster,
    # Aneuploidy_Score,
    # MSI_MANTIS_Score,
    # MSIsensor_Score,
    # TMB_.nonsynonymous.,
    # Buffa_Hypoxia_Score,
#     Ragnum_Hypoxia_Score) %>%
#   as.data.frame() 
# 
# #ydata.raw<-read.csv("ydata.csv",header = T,row.names = 1)
# ##导入time\status
# #xdata.NoneScale.raw.select$time <- ydata.raw$time
# #xdata.NoneScale.raw.select$status <- ydata.raw$status
# 
# # #调整顺序
# cols <- colnames(cox.select)
# new_cols <- c(cols[length(cols)],cols[length(cols)-1], cols[1:(length(cols) - 2)])
# # 
# # # 然后将 dataframe 按照新的列名顺序排列
# cox.select <- cox.select[, new_cols]

##计算
library(survival)
pFilter=0.1 
outResult=data.frame()
sigGenes=c("time","status")
for(i in colnames(cox[,3:ncol(cox)])){
  mycox <- coxph(Surv(time, status) ~ cox[,i], data = cox)
  mycoxSummary = summary(mycox)
  pvalue=mycoxSummary$coefficients[,"Pr(>|z|)"]
  if(pvalue<pFilter){
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=mycoxSummary$conf.int[,"exp(coef)"],
                          L95CI=mycoxSummary$conf.int[,"lower .95"],
                          H95CI=mycoxSummary$conf.int[,"upper .95"],
                          pvalue=mycoxSummary$coefficients[,"Pr(>|z|)"])
    )
    
  }
}

##写出结果文件
write.csv(outResult, file='outResult_表型_单因素cox.csv')

#############单因素森林图
library(survival)
library(survminer)
library(coin)
library(ggpmisc)
library(dplyr)
library(forestplot)
outResult<-read.csv("outResult_表型_单因素cox_森林图.csv", header = T, sep = ",")
#森林图绘制：
photo <- outResult %>%
  forestplot(labeltext = c(id,HR,pvalue), 
             clip = c(0.5, 3),#森林图中置信区间上下限设置
             xlog = F,
             zero=1.07)%>%
  fp_set_style(box = "#ff8216", #权重方块颜色
               line = "#ff8216") %>%#是否对X轴取对数
  fp_set_zebra_style("#f9f1dd") %>% #更改为斑马样式，并修改汇总菱形的填充颜色
  fp_add_header(id = c("","gene"),
                HR =c("","HR"),
                pvalue = c("","pvalue"))
pdf(file="cox_forestPlot.pdf",width = 7,height =5)
print(photo)
dev.off()


##----------------------------cox多因素---------------------------------------#
#做多因素回归（临床表型）=====================================
myjr<-read.csv('INDEX&表型单多因素cox.csv', header = T, row.names = 1, sep=",")


library(coin)
library("survival")
library("survminer")

#载入并查看数据集
lusc<-myjr

str(lusc)#该数据将所有变量都转换为数值型
#cox 回归分析
#res.cox <- coxph(Surv(time, status) ~ B_cells_naive+B_cells_memory+Plasma_cells+T_cells_CD8+T_cells_CD4_naive+T_cells_CD4_memory_resting+T_cells_CD4_memory_activated+T_cells_follicular_helper+Tregs+T_cells_gamma_delta+NK_cells_resting+NK_cells_activated+Monocytes+Macrophages_M0+Macrophages_M1+Macrophages_M2+Dendritic_cells_resting+Dendritic_cells_activated+Mast_cells_resting+Mast_cells_activated+Eosinophils+Neutrophils,data = lung)
#res.cox <- coxph(Surv(time, status) ~ CNR1+P_stage+T_stage+N_stage+M_stage+gender+race+lymphNodes_num+radiation_therapy+radiation_exposure,data = BRCA)
res.cox <- coxph(Surv(time, status) ~ AGE+
                   SEX+
                   M_STAGE+
                   N_STAGE+
                   T_STAGE+
                 index+
                 cluster+
                 Aneuploidy_Score+
                 MSI_MANTIS_Score+
                 MSIsensor_Score+
                 TMB_.nonsynonymous.+
                 Buffa_Hypoxia_Score+
                   Ragnum_Hypoxia_Score,data = lusc)

res.cox
summary(res.cox)   ###看显示框里的数据，整理后再导入


#绘制森林图(ggforest)


#读取文件绘制森林图
outResult<-read.csv("表型-多因素cox.csv", header = T, row.names=1)#读出来的数据第一列有1，2，3

#森林图绘制：
outResult %>%
  forestplot(labeltext = c(id,HR,Pvalue), #森林图中置信区间上下限设置
             clip = c(0.01, 2.8),
             xlog = F,
             zero=1.2)%>%
  fp_set_style(box = "#ff8216", #权重方块颜色
               line = "#ff8216") %>%#是否对X轴取对数
  fp_set_zebra_style("#f9f1dd") %>% #更改为斑马样式，并修改汇总菱形的填充颜色
  fp_add_header(id = c("","id"),
                HR =c("","HR"),
                Pvalue = c("","Pvalue"))

ggsave("多因素cox_森林图.pdf",width=6,height= 3.5)
dev.off()

#------------------high&low_risk--------高低风险组-ssgsea&免疫浸润--------------#
library(GSEABase)
library(GSVA)

library(tidyverse)
library(dplyr)
######
sessionInfo("GSVA")

deg_counts <- read.csv("low_risk_ssgsea.csv", sep=",", header= T)#(自己的差异基因-基因名+所有的组——fpkm)
######导入数据data1之后，把第一列的名字自动加序号并变为行名
row.names(deg_counts)<-make.names(deg_counts[,1],TRUE)
deg_counts<-deg_counts[,-1] #删除第一列

gene_level_expression_t <-deg_counts%>% t()%>%as.data.frame()

#####log
logTPM <- log(gene_level_expression_t)
##logTPM转置并变为matrix
logTPM <-logTPM%>% t()%>%as.data.frame()%>%as.matrix()

##l变为list
l<- immunecell%>% as.list()
#3. 开始进行ssGSEA
ssgsea<- gsva(logTPM, l, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)#可以改为“gsva”, “ssgsea”, “zscore”, “plage”算法

# Min-Max标准化是指对原始数据进行线性变换，将值映射到[0，1]之间
# 这里是将每个样本中不同的免疫细胞比例标准化到0-1之间
ssgsea.1 <- ssgsea
for (i in colnames(ssgsea)) {
  #i <- colnames(ssgsea)[1]
  ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
  
}
apply(ssgsea.1[,1:6], 2, range)
##保存结果
write.csv(ssgsea,file="low_risk_免疫浸润_ssgsea（原始评分）.csv")         #####原始评分结果
write.csv(ssgsea.1,file="low_risk_免疫浸润_ssgsea（标准化）.csv")   ####标准化的结果

#堆积柱状图
library(ggplot2)
#excel转置，样本在第一列，免疫浸润细胞在第一行
df   <- read.csv("low_risk_免疫浸润_ssgsea（标准化）.csv", header = T,sep = ",")
#df样本在第一列，免疫浸润细胞在第一行
#读入的数据有行号

#导入样本分组信息
Groups=c(rep("High-risk",219))#构建样本分组信息#######z,

df$group <-  Groups

# 转换数据框为长格式
df_long <- tidyr::pivot_longer(df, cols = c("Activated.CD8.T.cell"  ,         
                                             "Central.memory.CD8.T.cell" ,      "Effector.memeory.CD8.T.cell" ,   
                                            "Activated.CD4.T.cell"   ,         "Central.memory.CD4.T.cell"  ,    
                                            "Effector.memeory.CD4.T.cell",     "T.follicular.helper.cell"    ,   
                                           "Gamma.delta.T.cell" ,             "Type.1.T.helper.cell"          , 
                                            "Type.17.T.helper.cell" ,          "Type.2.T.helper.cell"          , 
                                            "Regulatory.T.cell"  ,             "Activated.B.cell"               ,
                                           "Immature..B.cell"     ,           "Memory.B.cell"                  ,
                                            "Natural.killer.cell"  ,           "CD56bright.natural.killer.cell" ,
                                            "CD56dim.natural.killer.cell" ,    "Myeloid.derived.suppressor.cell",
                                            "Natural.killer.T.cell" ,          "Activated.dendritic.cell"       ,
                                            "Plasmacytoid.dendritic.cell" ,    "Immature.dendritic.cell"        ,
                                            "Macrophage"  ,                    "Eosinophil"                     ,
                                            "Mast.cell"  ,                     "Monocyte"  ,                     
                                            "Neutrophil"
), names_to = "variable", values_to = "value")



#转换成百分比格式
library(plyr)
library(tidyverse)

data2<-  df_long 
data3 = ddply(data2,'variable',transform,percent_con=value/sum(value)*100)
data3

# 绘制堆积柱状图
mypalette <- colorRampPalette(brewer.pal(12,"Set2"))

ggplot(data2, aes(x = X, y = value, fill = variable)) +
  geom_bar(stat="identity",position="fill", width=0.7,size=0.25)+ #position可以根据数据更改填充方式（有stack和fill）
  geom_col(position = "fill") + #stack ：：同上
  labs(title = "",fill = "variable",x = "",y = "Estiamted Proportion")+
  scale_y_continuous(breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+
  theme_bw() + xlab("")+
  theme( axis.text.x=element_text(colour="black",size=14), #设置x轴刻度标签的字体属性
         axis.text.y=element_text(size=14,face="plain"), #设置x轴刻度标签的字体属性
         axis.title.y=element_text(size = 14,face="plain"), #设置y轴的标题的字体属性
         axis.title.x=element_text(size = 14,face="plain"))+
  scale_fill_manual(values = mypalette(28))


  ggsave("low_risk_免疫浸润堆积图(标准化).pdf",width=16,height=5)

  ###-----------------------------------(箱线图)免疫相关基因在高低组的表达------------------------------#
  data <- read.csv('免疫相关基因exp.csv',header = T,row.names = 1)
  library(ggsignif)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  # 转换数据框为长格式  把要转换的列名写在括号里，用逗号隔开
  df_long <- tidyr::pivot_longer(data, cols = c(CD8A,
                                                GZMH,
                                                GZMK,
                                                GZMM,
                                                TGFB1,
                                                TGFB2,
                                                TGFB3,
                                                TNF,
                                                GZMA,
                                                GZMB,
                                                CD83,
                                                PDCD1,
                                                CD274,
                                                CTLA4,
                                                IL27,
                                                IL27RA,
                                                BCL2,
                                                HLA.DMA,
                                                HLA.DMB,
                                                HLA.DOB,
                                                HLA.DOA,
                                                HLA.DPA1,
                                                HLA.DPB1,
                                                HLA.DPB2,
                                                HLA.DQA1,
                                                HLA.DQB1,
                                                HLA.DQB2,
                                                HLA.DRA,
                                                HLA.DRB1,
                                                HLA.DRB5,
                                                HLA.E),
                                 names_to = "cell", values_to = "value")
  
  ggplot(df_long, aes(x = risk, y = log2(value), fill = risk)) +   #log2(CNR1+0.000001)
    stat_boxplot(geom = "errorbar",width=0.1,color="black")+
    facet_wrap(.~cell,scales="free")+
    geom_boxplot(size=0.5, alpha=0.8) + #删去了,varwidth = TRUE，就不会样本大的宽度就越大
    stat_compare_means(method = "t.test",label = "p.signif")+             #加显著性
    #coord_cartesian(ylim =  c(-7, 4))+
    #geom_jitter(mapping=aes(x=risk,y=log2(GZMB+0.000001)),colour ="#760077", alpha = 0.2,size=1.5)+#散点
    #geom_jitter(position=position_jitter(width=0.1,height=0.2), size=1.5,alpha=0.3,colour="#8b1a1a") + 
    scale_fill_manual(values=c("#ff8216","#018a89")) +
    #ggtitle("ZNF862")+#设置总的标题
    theme_bw()+
    xlab("")+      #不显示横坐标的title
    theme(legend.position="none", #不需要图例
          axis.text.x = element_text(size = 14, angle=45, hjust = 1),## 调整横坐标字体大小和旋转角度
          axis.text.y = element_text(size = 14),   ## 调整纵坐标字体大小
          axis.title.y=element_text(size = 14,face="plain"), #设置y轴的标题的字体属性
          axis.title.x=element_text(size = 14,face="plain"), #设置x轴的标题的字体属性
          plot.title = element_text(size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
          panel.grid.major = element_blank())#不显示网格线
  
  ggsave("免疫相关基因exp.pdf",width=8,height=14)


  #差异分析-----------------high-risk&low-risk----------#
  
  #DESeq2进行差异分析
  library(DESeq2)#调用包
  data<-read.csv("H&L_DEGs_exp.csv",header = T,row.names = 1)#读入count值表格，一般为6列，数值必须为整数
  
  group<-read.csv("H&L_DEGs_exp_group.csv",header = T,row.names = 1)#读入分组文件
  head(group)
  
  dds <- DESeqDataSetFromMatrix(countData=round(data), 
                                colData=group, 
                                design=~group)
  #5.1.指定因子水平（此处需把处理组往前放，对照组往后放）
  dds$group <- factor(dds$group, levels = c("High-risk","Low-risk"))
  #差异分析，必要步骤
  dds1<-DESeq(dds)
  
  resultsNames(dds1) #查看结果
  #dds<-DESeqDataSetFromMatrix(countData=data,colData = group,design = ~group)#获取dss文件，countData为count表格，colData是分组表格，design是差异比较组
  dds
  dds1<-DESeq(dds)#数据标准化，必要步骤
  resultsNames(dds1)#查看结果
  dds1$condition#默认后面比前面
  res<-results(dds1)#结果转换，必要步骤
  summary(res)#结果概要
  
  write.csv(res,file="High-risk&Low-risk_DESeq.csv")#导出差异分析结果
  
  deg<-read.csv("High-risk&Low-risk_DESeq.csv",header = T,row.names = 1)
  deg$diff<-"no"#增加一列diff并全部用no填充
  head(deg)
  
  deg$diff[(deg$log2FoldChange>0.5) & (deg$padj<0.05)]<-"up" #将logFC>1且FDR<0.05的行标注为up
  deg$diff[(deg$log2FoldChange<(-0.5)) & (deg$padj<0.05)]<-"down" #将logFC小于-1且FDR<0.05的行标注为up
  
  as.data.frame(table(deg$diff)) #统计上下调基因数目
  write.csv(deg,"High-risk&Low-risk_DESeq.csv") 
  
  
  #ggplot2画火山图
  #deg<-read.csv("H&L.all.limmaOur.csv",header=T,sep=",")#读入差异分析结果文件
  library(ggplot2)
  library(ggrepel)
  # deg$diff<-"no"#增加一列diff并全部用no填充
  # colnames(deg)
  # # [1] "gene_id"     "TS32_tpm"    "TS34_tpm"    "fc"          "log2fc"     
  # # [6] "pvalue"      "padjust"     "significant" "regulate"    "diff"  
  deg$diff[(deg$logFC>0.5) & (deg$pvalue<0.05)]<-"up" #将logFC>1且FDR<0.05的行标注为up
  deg$diff[(deg$logFC<(-0.5)) & (deg$pvalue<0.05)]<-"down" #将logFC小于-1且FDR<0.05的行标注为up
  # 
   as.data.frame(table(deg$diff)) #统计上下调基因数目
  
  deg$diff<-factor(deg$diff,
                   levels = c('up','no','down'),
                   labels = c("Up (145)","No (33283)","Down (569)")) #将顺序设定为'up','no','down'，同时根据上面的结果标注DEGs的数目
  
  #deg$name<-deg$gene_id #新增一列name，并且把gene的内容复制给name
  
  #deg$name[(deg$log2FoldChange>(-4.5))&(deg$log2FoldChange<(5.5))]<- NA #筛选logFC小于(-2.5)和logFC大于(2.7)的基因，将不在该范围内的name值用NA填充
  
  
  
  
  # cmp <- deg %>% select(name)
  # write.csv(cmp,"1.csv")
  
  ggplot(deg)+
    geom_point(aes(x=log2FoldChange,y= -1*log10(padj),color=diff)) +#画点，用diff进行填充颜色
    geom_vline(xintercept = c(0.5,-0.5),linetype="dashed",color="grey")+ #X轴辅助线
    geom_hline(yintercept =-log10(0.05),linetype="dashed",color="grey")+ #Y轴辅助线
    scale_color_manual(values=c("#ff8216","#EEE5D9","#018a89"))+ #将颜色改为"red","grey","green"，对应上面的'up','no','down'
    guides(color=guide_legend(title = " C1&C3|log2FoldChange|>0.5 & padj<0.05",override.aes = list(size=1)))+ #更改图例，将比较组，差异基因筛选标准和数目标注到图中
    scale_y_continuous(expand = c(0,0))+#坐标轴相交于原点
    theme_bw()+#主题
    theme(legend.position="top")+
    xlab("log2FoldChange")+
    ylab("-Log10(pvalue)")+
    # geom_text_repel(aes(x= -1*log10(padjust),y=log2fc,label=name),
    #                 max.overlaps = 10000,                    # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之。
    #                 size=3, # 字体大小box.padding=unit(0.5,'lines'),           # 标记的边距
    #                 point.padding=unit(0.1, 'lines'), 
    #                 segment.color='black',                   # 标记线条的颜色
    #                 show.legend=FALSE)+ #使用ggrepel包进行目标基因标注，可防止相近的点重叠到一起
    theme(panel.grid=element_blank())#去掉次要线条
  
  ggsave("  High-risk&Low-risk_(log2FoldChange大于0.5).pdf",width=4.5,height =5)
  #write.csv(deg,"DEGs_up.down.csv")
  
  #write.csv(immunecell,"immunecell.csv")


 
  #------------------high&low_risk--------高低风险组-114个代谢通路--------------#
  library(GSEABase)
  #BiocManager::install("HDF5Array")
  library(GSVA)
  
  library(tidyverse)
  library(dplyr)
  ######
  sessionInfo("GSVA")
  #1. 获取geneSets 基因背景集
  #load(file = "C00基因集（三个）.symbols.Rdata")##(处理好的gmt-rdata数据)
  
  #2. 读取基因的表达矩阵
  #rm(list=ls())
  ## 载入数据-如果数据里有NA，需要先删掉
  deg_counts <- read.csv("高低风险组-ssgsea（代谢&hallmark）.csv", sep=",", header= T)#(自己的差异基因-基因名+所有的组——fpkm)
  ######导入数据data1之后，把第一列的名字自动加序号并变为行名
  row.names(deg_counts)<-make.names(deg_counts[,1],TRUE)
  deg_counts<-deg_counts[,-1] #删除第一列
  
  gene_level_expression_t <-deg_counts%>% t()%>%as.data.frame()
  
  #####log
  logTPM <- log(gene_level_expression_t)
  ##logTPM转置并变为matrix
  logTPM <-logTPM%>% t()%>%as.data.frame()%>%as.matrix()
  
  ##l变为list
  #l<- immunecell%>% as.list()
  #3. 开始进行ssGSEA
  ssgsea<- gsva(logTPM, l, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)#可以改为“gsva”, “ssgsea”, “zscore”, “plage”算法
  
  # Min-Max标准化是指对原始数据进行线性变换，将值映射到[0，1]之间
  # 这里是将每个样本中不同的免疫细胞比例标准化到0-1之间
  ssgsea.1 <- ssgsea
  for (i in colnames(ssgsea)) {
    #i <- colnames(ssgsea)[1]
    ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
    
  }
  apply(ssgsea.1[,1:6], 2, range)
  ##保存结果
  write.csv(ssgsea,file="risk_114个代谢通路_ssgsea（原始评分）.csv")         #####原始评分结果
  write.csv(ssgsea.1,file="risk_114个代谢通路_ssgsea（标准化）.csv")   ####标准化的结果
  
  #------------------high&low_risk--------高低风险组-hallmark--------------#
  library(GSEABase)
  #BiocManager::install("HDF5Array")
  library(GSVA)
  
  library(tidyverse)
  library(dplyr)
  ######
  sessionInfo("GSVA")
  #1. 获取geneSets 基因背景集
  #load(file = "C00基因集（三个）.symbols.Rdata")##(处理好的gmt-rdata数据)
  
  #2. 读取基因的表达矩阵
  #rm(list=ls())
  ## 载入数据-如果数据里有NA，需要先删掉
  deg_counts <- read.csv("高低风险组-ssgsea（代谢&hallmark）.csv", sep=",", header= T)#(自己的差异基因-基因名+所有的组——fpkm)
  ######导入数据data1之后，把第一列的名字自动加序号并变为行名
  row.names(deg_counts)<-make.names(deg_counts[,1],TRUE)
  deg_counts<-deg_counts[,-1] #删除第一列
  
  gene_level_expression_t <-deg_counts%>% t()%>%as.data.frame()
  
  #####log
  logTPM <- log(gene_level_expression_t)
  ##logTPM转置并变为matrix
  logTPM <-logTPM%>% t()%>%as.data.frame()%>%as.matrix()
  
  ##l变为list
  #l<- immunecell%>% as.list()
  #3. 开始进行ssGSEA
  ssgsea<- gsva(logTPM, l, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)#可以改为“gsva”, “ssgsea”, “zscore”, “plage”算法
  
  # Min-Max标准化是指对原始数据进行线性变换，将值映射到[0，1]之间
  # 这里是将每个样本中不同的免疫细胞比例标准化到0-1之间
  ssgsea.1 <- ssgsea
  for (i in colnames(ssgsea)) {
    #i <- colnames(ssgsea)[1]
    ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
    
  }
  apply(ssgsea.1[,1:6], 2, range)
  ##保存结果
  write.csv(ssgsea,file="risk_hallmark_ssgsea（原始评分）.csv")         #####原始评分结果
  write.csv(ssgsea.1,file="risk_hallmark_ssgsea（标准化）.csv")   ####标准化的结果

  
  
  ##-------------------------相关性分析---------6个癌症相关基因&index相关性分析.csv--------------------##
  #install.packages("ggpmisc")
  library(ggpmisc)
  
  library(ggplot2)
  cox  <-   read.csv("6个癌症相关基因&index相关性分析.csv", header = T, sep = ",")
  head(cox)
  # # 转换数据框为长格式  把要转换的列名写在括号里，用逗号隔开
  df_long <- tidyr::pivot_longer(cox, cols = c(FGFR1,ESR1,PDGFB,PDGFRB,FLT3,NTRK1),
                                 names_to = "cell", values_to = "value")

  
  ggplot(df_long, aes(x=index,y=value))+
    geom_point(aes(color="#018a89",size=0.7))+
    geom_smooth(method = "lm",se = T)+
    stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(p.value.label),sep = "~~~")),formula = y~x,parse = TRUE, size=5)+
    stat_correlation(mapping = use_label("R"),method = "pearson",label.x = "right",size=5)+
    labs(x="value",y=NULL,size=18,face="bold")+
    scale_color_manual(values = c("#018a89"))+
    scale_fill_manual(values = c("#018a89"))+
    facet_grid(cols = vars(cell),scales="free")+
    scale_y_continuous(expand = c(0, 0))+
    theme(
      axis.text.x = element_text(size = 12, hjust = 1),## 调整横坐标字体大小和旋转角度
      axis.text.y = element_text(size = 12),   ## 调整纵坐标字体大小
      axis.title.x=element_text(size = 18,face="bold"),
      axis.title.y=element_text(size = 12,face="plain"), #设置y轴的标题的字体属性
      plot.title = element_text(size=15,face="bold",hjust = 0.3),
      legend.position="none", #不需要图例
      #legend.position="bottom",  #图例在底下
      panel.background=element_rect(fill='transparent',color='black'),
      panel.grid =element_blank(),
      #axis.text = element_text(colour="black"),
      legend.key = element_blank(),
      legend.title=element_blank())#+
  ggsave("6个癌症相关基因&index相关性分析.pdf",width=20,height= 4)
  
  #————原癌基因抑癌基因—————————————————————————————————————————————————————————————————————————————————————————————
  
  ##筛选出risk及原癌基因抑癌基因相关基因的列
  ###正图
  xdata.NoneScale.risk.select <- xdata.NoneScale.risk%>% 
    as.data.frame() %>% 
    dplyr::select(
      FLT3,ATM,NRAS,BRCA2,PALB2,PDGFRA,ESR1,RET,BRAF,NTRK3,CHEK2,index,
      risk) %>%
    as.data.frame() 
  #####附图
  xdata.NoneScale.risk.select <- xdata.NoneScale.risk%>% 
    as.data.frame() %>% 
    dplyr::select(  ABL1,ALK,BARD1,BRCA1,BRIP1,CDK12,CHEK1,EGFR,ERBB2,
                    EZH2,FANCL,FGFR1,FGFR2,FGFR3,IDH1,IDH2,KIT,KRAS,MET,NF1,
                    NTRK1,NTRK2,PDGFB,PDGFRB,PIK3CA,RAD51C,
                    RAD54L,ROS1,SMARCB1,TSC1,TSC2,index,
                    risk) %>%
    as.data.frame() 
  
  
  ydata.risk<-ydata.risk%>%
    as.data.frame() 
  
  # 转换数据框为长格式
  xdata.NoneScale.risk.select_long2 <- tidyr::pivot_longer(
    xdata.NoneScale.risk.select, 
    cols = c( ABL1, ALK,BARD1,BRCA1,BRIP1,CDK12,CHEK1,EGFR,ERBB2,
              EZH2,FANCL,FGFR1,FGFR2,FGFR3,IDH1,IDH2,KIT,KRAS,MET,NF1,
              NTRK1,NTRK2,PDGFB,PDGFRB,PIK3CA,RAD51C,
              RAD54L,ROS1,SMARCB1,TSC1,TSC2
    ), names_to = "gene", values_to = "value")
  
  ###########基因-box图
  pdf_file <- paste("原癌基因抑癌基因_lusc-BoxPlot.pdf")
  
  pdf(pdf_file,width=9,height= 12)
  ggplot(xdata.NoneScale.risk.select_long2, aes(x =factor(risk), y=log2(value) ,fill =factor(risk))) +
    geom_boxplot(size=0.5,varwidth = F, alpha=0.8) + 
    geom_jitter(position=position_jitter(width=0.15,height=0.2), size=0.8,alpha=0.3,colour="#bfdbf7") + 
    facet_wrap(gene~.,scales="free",ncol = 6)+
    scale_fill_manual(values=c("#018a89","#ff8216")) +
    theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
          axis.text.x = element_text(angle = 30,vjust = 0.9,hjust = 1), panel.grid =element_blank(),axis.text = element_text(colour="black"),
          strip.background = element_rect(colour="black", fill="white"))+
    geom_signif(
      comparisons = list( c("Low-risk", "High-risk")),
      map_signif_level = TRUE, textsize = 6,step_increase = 0.1
    ) 
  dev.off()
  
  #————————————————————————————————————————————————————————————————————————————————————
  ##筛选出risk及原癌基因抑癌基因相关基因的列
  xdata.NoneScale.risk.select <- xdata.NoneScale.risk%>% 
    as.data.frame() %>% 
    dplyr::select(ABL1,ALK,ATM,BARD1,BRAF,BRCA1,BRCA2,BRIP1,CDK12,CHEK1,CHEK2,EGFR,ERBB2,
                  ESR1,EZH2,FANCL,FGFR1,FGFR2,FGFR3,FLT3,IDH1,IDH2,KIT,KRAS,MET,NF1,NRAS,
                  NTRK1,NTRK2,NTRK3,PALB2,PDGFB,PDGFRA,PDGFRB,PIK3CA,RAD51C,
                  RAD54L,RET,ROS1,SMARCB1,TSC1,TSC2) %>%
    as.data.frame() 
  ydata.risk<-ydata.risk%>%
    as.data.frame() 
  xdata.NoneScale.risk.select_t <- xdata.NoneScale.risk.select %>% t()%>%as.data.frame()
  xdata.NoneScale.risk.select_t = as.data.frame(lapply(xdata.NoneScale.risk.select_t,as.numeric))
  rownames(xdata.NoneScale.risk.select_t)=colnames(xdata.NoneScale.risk.select)
  xdata.NoneScale.risk.select_t <- xdata.NoneScale.risk.select_t %>%as.matrix()
  colnames(xdata.NoneScale.risk.select_t ) <- gsub("\\.", "-", colnames(xdata.NoneScale.risk.select_t ))
  
  #可视化-热图
  library(pheatmap)
  library(RColorBrewer)
  #####先加好risk分组信息：
  ##把risk列提出来变为sample_info
  sample_info<-xdata.NoneScale.risk$risk%>%as.data.frame()
  ##将ssgsea.1的列名加到sample_info的行名上
  rownames(sample_info) = colnames(xdata.NoneScale.risk.select_t)
  ##作图
  pdf(file="pheatmap_hallmark.pdf",width = 10,height =10)
  pheatmap(xdata.NoneScale.risk.select_t,scale="row",border=NA,
           color=colorRampPalette(c("#062E14","white","#ff8216"))(500),
           cluster_row=T,cluster_cols=F,gaps_col =219,
           show_rownames=T,show_colnames = F,treeheight_col =10,treeheight_row =15,legend = T,
           annotation_col =sample_info )
  dev.off()
  
  gene <- "genelist"
  pdf_file <- paste(gene ,"index_癌症相关基因_lusc-pointPlot.pdf", sep="_")
  pdf(pdf_file,width=8,height= 16)
  
  ggplot( xdata.NoneScale.risk.select_long2 , aes(x=index,y=value))+
    geom_point(size=0.8,colour ="#ff8216" )+
    geom_smooth(method = "lm",se = T,colour="#ff8216",fill = "#ff8216")+
    stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(p.value.label),sep = "~~~")),formula = y~x,parse = TRUE, size=2,colour ="#219ebc")+
    stat_correlation(mapping = use_label("R"),method = "pearson",label.x="right",npcy = 0.8,size=2.5,colour ="#219ebc")+
    labs(x=" ",y=NULL)+
    facet_wrap(.~gene,scales="free",ncol=4)+
    scale_y_continuous(expand = c(0, 0))+
    theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
          panel.grid =element_blank(),
          axis.text = element_text(colour="black"),legend.key = element_blank(),
          legend.title=element_blank())+
    labs(fill=' ')
  dev.off()  
  
  #————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
  
  
  gene <- "genelist"
  pdf_file <- paste(gene ,"index_原癌基因抑癌基因_lusc-pointPlot_附图.pdf", sep="_")
  pdf(pdf_file,width=9,height= 10)
  
  ggplot( xdata.NoneScale.risk.select_long2 , aes(x=index,y=value))+
    geom_point(size=0.8,colour ="#219ebc" )+
    geom_smooth(method = "lm",se = T,colour="#219ebc",fill = "#219ebc")+
    stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(p.value.label),sep = "~~~")),formula = y~x,parse = TRUE, size=2,colour ="#219ebc")+
    stat_correlation(mapping = use_label("R"),method = "pearson",label.x="right",npcy = 0.8,size=2.5,colour ="#219ebc")+
    labs(x=" ",y=NULL)+
    facet_wrap(.~gene,scales="free",ncol=6)+
    scale_y_continuous(expand = c(0, 0))+
    theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
          panel.grid =element_blank(),
          axis.text = element_text(colour="black"),legend.key = element_blank(),
          legend.title=element_blank())+
    labs(fill=' ')
  dev.off()  
  
  
  
  #——————确认完之后，把分组信息加到原来表格————————————————————————————————————————————————————————————————————————————————————————————————————
  clust <- cmoic.lusc$clust.res
  ydata.risk$cluster<-as.numeric(ydata.risk$cluster)
  # Get matches
  ydata.risk.tmp <- ydata.risk[rownames(ydata.risk) %in% 
                                 rownames(clust),]
  ydata.risk.tmp     <- ydata.risk.tmp [rownames(clust), ]
  #####把分组信息加到cmoic.coad并重命名为cmoic.coad_TLS
  clust$cluster <- ydata.risk.tmp$cluster
  cmoic.lusc_TLS<-cmoic.lusc
  cmoic.lusc_TLS$clust.res$clust<-ydata.risk.tmp$cluster
  
  '''''
sil<-cmoic.coad_TLS$sil %>% as.data.frame()

clust <- cmoic.coad$clust.res
clust <- cmoic.coad_TLS$clust.res

clust$clust <-as.numeric(clust$clust)
# Get matches
clust.tmp <- clust[rownames(clust) %in% 
                               rownames(sil),]
clust.tmp     <-clust.tmp [rownames(sil), ]
#####把分组信息加到cmoic.coad并重命名为cmoic.coad_TLS
clust.tmp$neighbor <- clust.tmp$clust
clust.tmp$cluster <- clust.tmp$clust
cmoic.coad_TLS2 <- cmoic.coad_TLS
cmoic.coad_TLS2$sil[,1] <-clust.tmp$cluster
cmoic.coad_TLS2$sil[,2] <-clust.tmp$neighbor

sil2<-cmoic.coad_TLS2$sil %>% as.data.frame()
#cmoic.coad_TLS$sil[,1]
'''''
  ###保存数据
  save(iClusterBayes.res,moic.res.list,cmoic.lusc,cmoic.lusc_TLS, file = "4-6.cmoic.lusc_TLS.rda")
  
  
  #7.2) get multi-omics heatmap based on clustering result
  
  #全基因组热图被广泛用于以图形方式显示大型基因组数据集中的潜在模式。
  #它们已被用于揭示样本/基因如何聚集在一起的信息，并提供对潜在样本偏差或其他伪影的见解。
  #在此，MOVICS提供了getMoHeatmap，以直观地处理预聚类的多组学数据，并生成一个精致的热图以满足发布需求。
  #在使用getMoHeatmap（）之前，应使用getStdiz（）函数正确处理组学数据，该函数返回存储规范化组学数据的列表。
  #奥密克戎数据，尤其是表达数据（如RNA和蛋白质），应居中（centerFlag=TRUE）或按比例（scaleFlag=TRUE）或z评分（居中和按比例）。
  #通常，DNA甲基化数据（β矩阵范围从0到1）和体细胞突变（0和1二进制矩阵）不应归一化。
  #然而，根据M=log2β1-β的公式将甲基化β值转换为M值是一个很好的选择，因为它在可视化中具有更强的信号，
  #并且M值适合标准化。该函数还为连续组学数据提供了半宽度的自变量；这样的自变量用于截断归一化后的“极值”；
  #特别地，超过半宽度边界的归一化值将被半宽度代替，这有利于热图中的映射颜色。
  
  load("4-6.cmoic.lusc_TLS.rda")
  
  # convert beta value to M value for stronger signal
  indata <- mo.data
  indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))
  
  # data normalization for heatmap
  plotdata <- getStdiz(data       = indata,
                       halfwidth  = c(2,2,2,NA), # no truncation for mutation
                       centerFlag = c(T,T,T,F), # F:no center for mutation
                       scaleFlag  = c(T,T,T,F)) #F: no scale for mutation
  
  #正如我前面提到的，一些算法也提供了特征选择；这些选定的特征显示出与其他组学数据的复杂串扰，
  #可能具有特殊的生物学意义，从而导致癌症的异质性。因此，我在下面展示了如何使用getMoHeatmap（）基于单一
  #算法（例如，iClusterBayes）和选定的特征生成综合热图。但是，首先必须提取特征：
  
  feat   <- iClusterBayes.res$feat.res
  feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
  feat2  <- feat[which(feat$dataset == "lncRNA.expr"),][1:10,"feature"]
  feat3  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
  feat4  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
  annRow <- list(feat1, feat2, feat3, feat4)
  #iClusterBayes.res中包含的feat.res是根据每个组学数据的特征后验概率进行排序的。
  #通过这种方式，为每个组学数据选择前10个特征，并生成一个特征列表，并将其命名为
  #anRow for heatmap原始注释（别担心，我稍后将讨论如何为样本附加列注释）。
  
  # set color for each omics data
  # if no color list specified all subheatmaps will be unified to green and red color pattern
  mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
  lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
  meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
  mut.col    <- c("grey90" , "black")
  col.list   <- list(mRNA.col, lncRNA.col, meth.col, mut.col)
  
  # comprehensive heatmap (may take a while)
  getMoHeatmap(data          = plotdata,
               row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
               is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
               legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
               clust.res     = cmoic.lusc_TLS$clust.res, # cluster results，如果此处改为：iClusterBayes.res$clust.res，则意为用单个iClusterBayes聚类的结果来作为分群数据
               clust.dend    = NULL, # no dendrogram
               show.rownames = c(F,F,F,F), # specify for each omics data
               show.colnames = FALSE, # show no sample names
               annRow        = annRow, # mark selected features
               color         = col.list,
               annCol        = NULL, # no annotation for samples
               annColors     = NULL, # no annotation color
               width         = 10, # width of each subheatmap
               height        = 5, # height of each subheatmap
               fig.name      = "7.2COMPREHENSIVE HEATMAP OF CONSENSUSMOIC_无表型")
  
  '''''
#当然，由于moic.res.list中总共存储了10个结果，您也可以选择其中的任何一个来创建热图。
#在这里，我选择了COCA，它也返回如下所示的样本树状图：
# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = moic.res.list$COCA$clust.res, # cluster results
             clust.dend    = moic.res.list$COCA$clust.dend, # show dendrogram for samples
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF COCA")
'''''
  #现在回到cmoic.coad的共识结果，它集成了10个算法，这次还提供了样本的注释来生成热图。
  #由于getMoHeatmap（）的核心函数基于ComplexHeatmap R包，因此在创建注释时，
  #应始终使用circize:：colorRamp2（）函数为连续变量（例如，本例中的年龄）生成颜色映射函数。
  
  # extract PAM50, pathologic stage and age for sample annotation
  annCol    <- surv.info[,c("PATH_M_STAGE", "PATH_N_STAGE","PATH_T_STAGE", "AGE","SEX"), drop = FALSE]
  
  # generate corresponding colors for sample annotation
  annColors <- list(age    = circlize::colorRamp2(breaks = c(min(annCol$AGE),
                                                             median(annCol$AGE),
                                                             max(annCol$AGE)), 
                                                  colors = c("#0000AA", "#555555", "#AAAA00")),
                    PATH_T_STAGE = c("T1"    = "green",
                                     "T1A"    = "green",
                                     "T1B"    = "green",
                                     "T2"    = "blue",
                                     "T2A"    = "blue",
                                     "T2B"    = "blue",
                                     "TB"    = "blue",
                                     "T3"    = "red",
                                     "T4"    = "yellow", 
                                     PATH_M_STAGE = c("M0"    = "green",
                                                      "M1"    = "blue",
                                                      "M1A"    = "blue",
                                                      "M1B"    = "blue",
                                                      "MX"    = "black"),
                                     PATH_N_STAGE = c("N0"    = "green",
                                                      "N1"    = "blue",
                                                      "N2"    = "red",
                                                      "N3"    = "yellow",
                                                      "NX"    = "black"),
                                     SEX = c("Male"    = "green",
                                             "Female"    = "blue")
                    ))
  
  # comprehensive heatmap (may take a while)
  getMoHeatmap(data          = plotdata,
               row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
               is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
               legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
               clust.res     = cmoic.lusc_TLS$clust.res, # consensusMOIC results
               clust.dend    = NULL, # show no dendrogram for samples
               show.rownames = c(F,F,F,F), # specify for each omics data
               show.colnames = FALSE, # show no sample names
               show.row.dend = c(F,F,F,F), # show no dendrogram for features
               color         = col.list,
               annCol        = annCol, # annotation for samples
               annRow        = annRow, # mark selected features
               annColors     = annColors, # annotation color
               width         = 10, # width of each subheatmap
               height        = 5, # height of each subheatmap
               fig.name      = "7.2COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")
  
  
  # 1) survival comparison
  write.csv(surv.info,"surv.info.csv")
  surv.info <- read.csv("surv.info.csv",row.names = 1,header = T)#需要把生存状态和时间改了，时间要改为天，不是月
  
  
  
  
  
  load("4-6.cmoic.lusc_TLS.rda")
  surv.lusc <- compSurv(moic.res         = cmoic.lusc_TLS,
                        surv.info        = surv.info,
                        convt.time       = "m", # convert day unit to month
                        surv.median.line = "h", # draw horizontal line at median survival
                        xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                        fig.name         = "1.KAPLAN-MEIER CURVE OF CONSENSUSMOIC")
  #> --a total of 643 samples are identified.
  #> --removed missing values.
  
  print(lusc)
  # > print(surv.coad)
  # $fitd
  # Call:
  #   survdiff(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
  #            na.action = na.exclude)
  # 
  # N Observed Expected (O-E)^2/E (O-E)^2/V
  # Subtype=CS1 219      101     98.7    0.0536     0.116
  # Subtype=CS2 219       84     86.3    0.0613     0.116
  # 
  # Chisq= 0.1  on 1 degrees of freedom, p= 0.7 
  # 
  # $fit
  # Call: survfit(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
  #               na.action = na.exclude, error = "greenwood", type = "kaplan-meier", 
  #               conf.type = "plain")
  # 
  # n events median 0.95LCL 0.95UCL
  # CS1 219    101   53.6    37.3    70.2
  # CS2 219     84   53.0    37.2    69.8
  # 
  # $xyrs.est
  # Call: survfit(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res)
  # 
  # Subtype=CS1 
  # time n.risk n.event survival std.err lower 95% CI upper 95% CI
  # 1825     40      84    0.472  0.0436        0.394        0.566
  # 3650      7      13    0.221  0.0557        0.135        0.362
  # 
  # Subtype=CS2 
  # time n.risk n.event survival std.err lower 95% CI upper 95% CI
  # 1825     27      74    0.463  0.0495        0.376        0.571
  # 3650      5      10    0.204  0.0653        0.109        0.382
  # 
  # 
  # $overall.p
  # [1] 0.7336349
  
  #2)比较临床特征
  #然后比较不同亚型之间的临床特征。MOVICS 提供的功能compClinvar()可以总结连续变量和
  #分类变量并执行适当的统计测试。该函数可以给出.docx格式的表格，易于在医学研究论文中使用。
  surv.info2<-surv.info[,1:8]
  clin.coad <- compClinvar(moic.res      = cmoic.lusc_TLS,
                           var2comp      = surv.info2, # data.frame needs to summarize (must has row names of samples)
                           strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                           factorVars    = c("pstage","fustat","SEX","PATH_M_STAGE","PATH_N_STAGE"), # features that are considered categorical variables
                           nonnormalVars = c("futime","AGE"), # feature(s) that are considered using nonparametric test
                           exactVars     = c("pstage","PATH_M_STAGE","PATH_N_STAGE"), # feature(s) that are considered using exact test
                           doWord        = TRUE, # generate .docx file in local path
                           tab.name      = "2.SUMMARIZATION OF CLINICAL FEATURES")
  #> --all samples matched.
  # print(clin.coad$compTab)
  # level                      CS1                      CS2      p    test
  # 1                      n                             219                      219               
  # 2     AGE (median [IQR])            68.00 [62.00, 73.00]     68.00 [61.00, 74.00]  0.689 nonnorm
  # 3                SEX (%) Female               55 (25.1)                62 (28.3)   0.517        
  # 4                          Male              164 (74.9)               157 (71.7)                
  # 5       PATH_M_STAGE (%)     M0              180 (82.2)               179 (81.7)   0.750   exact
  # 6                            M1                1 ( 0.5)                 3 ( 1.4)                
  # 7                           M1A                1 ( 0.5)                 0 ( 0.0)                
  # 8                           M1B                1 ( 0.5)                 0 ( 0.0)                
  # 9                            MX               36 (16.4)                37 (16.9)                
  # 10      PATH_N_STAGE (%)     N0              137 (62.6)               141 (64.4)   0.081   exact
  # 11                           N1               66 (30.1)                48 (21.9)                
  # 12                           N2               11 ( 5.0)                24 (11.0)                
  # 13                           N3                2 ( 0.9)                 3 ( 1.4)                
  # 14                           NX                3 ( 1.4)                 3 ( 1.4)                
  # 15      PATH_T_STAGE (%)     T1               19 ( 8.7)                19 ( 8.7)   0.909        
  # 16                          T1A               12 ( 5.5)                10 ( 4.6)                
  # 17                          T1B               21 ( 9.6)                17 ( 7.8)                
  # 18                           T2               73 (33.3)                73 (33.3)                
  # 19                          T2A               35 (16.0)                42 (19.2)                
  # 20                          T2B               12 ( 5.5)                16 ( 7.3)                
  # 21                           T3               37 (16.9)                30 (13.7)                
  # 22                           T4               10 ( 4.6)                12 ( 5.5)                
  # 23            fustat (%)      0              118 (53.9)               135 (61.6)   0.122        
  # 24                            1              101 (46.1)                84 (38.4)                
  # 25                 time                    33.70 ± 32.21            29.39 ± 29.86  0.147        
  # 26 futime (median [IQR])        674.62 [349.64, 1435.55] 636.16 [241.15, 1096.76]  0.136 nonnorm
  
  
  '''''会报错，估计是cmoic.coad结果里有的结果没有改回2个分群，自己在网页做吧

#3）比较突变频率
#亚型特异性突变可能有望成为治疗靶点。因此，我比较了不同集群之间的突变频率。MOVICS提供compMut（）
#的函数，该函数处理二进制突变数据（0=野生；1=突变）。
#该函数对每个突变应用独立测试（即Fisher精确测试或χ2测试），并生成一个包含统计结果的表
#（如果指定，也可以使用.docx文件），并使用这些差异突变的基因创建OncoPrint。
pdf("3.ONCOPRINT FOR SIGNIFICANT MUTATIONS.pdf")
mut.coad <- compMut(moic.res     = cmoic.coad_TLS2,
                    mut.matrix   = coad.tcga$mut.status, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.005, # keep those genes that mutated in at least 5% of samples
                    p.adj.cutoff = 0.5, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                    innerclust   = F, # perform clustering within each subtype
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 10, 
                    height       = 5,
                    fig.name     = "3.ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "3.INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")
#> --all samples matched.
dev.off() 
'''''
  
  #4）比较总突变负荷
  #不用说，免疫疗法正在成为现代癌症治疗的支柱。最近的分析将肿瘤基因组景观与抗肿瘤免疫联系起来。
  #特别是，一项新的研究表明，肿瘤特异性基因组病变与免疫检查点激活以及患者对免疫疗法的反应程度和持续时间有关。
  #这些病变包含高突变负荷8和非整倍体9。为了量化这些可能影响免疫治疗的基因组改变，
  #MOVICS提供了两种功能来计算总突变负荷（TMB）和基因组改变部分（FGA）
  #。具体来说，TMB是指在肿瘤基因组中发现的突变数量，而FGA是受拷贝数增加或减少影响的基因组百分比。
  #这两个属性对基因研究人员都很有用，因为它们为他们提供了关于肿瘤基因组构成的更深入的信息。
  #让我从compTMB（）开始，向您展示如何使用这两个函数。首先，此函数的输入maf数据必须至少具有以下10列：
  maf <- lusc.tcga$maf
  head(maf)
  #然后运行compTMB（）。默认情况下，此函数仅在计算体细胞突变频率时考虑非同义变体，
  #包括Frame_Shift_Del、Frame_Shift_Ins、Splice_Site、Translation_Start_Site、Nonsense_mutation、
  #Nonstop_mutation、In_Frame_Del、In_FFrame_Ins和Missense_mutation。
  #除了计算TMB，该函数还将单核苷酸变体分为转变和转变（TiTv），并描述TMB和TiTv的分布。
  # compare TMB
  dev.new()
  tmb.lusc <- compTMB(moic.res     = cmoic.lusc_TLS,
                      maf          = maf,
                      rmDup        = TRUE, # remove duplicated variants per sample
                      rmFLAGS      = FALSE, # keep FLAGS mutations
                      exome.size   = 38, # estimated exome size
                      test.method  = "nonparametric", # statistical testing method
                      fig.name     = "4.DISTRIBUTION OF TMB AND TITV")
  # --all samples matched.
  # -Validating
  # --Removed 3010 duplicated variants
  # -Silent variants: 22237 
  # -Summarizing
  # --Possible FLAGS among top ten genes:
  #   TTN
  # MUC16
  # SYNE1
  # USH2A
  # -Processing clinical data
  # --Missing clinical data
  # -Finished in 1.079s elapsed (1.004s cpu) 
  # Wilcoxon rank sum test p value = 8.26e-02
  head(tmb.lusc$TMB.dat)
  TMB.dat <- tmb.lusc$TMB.dat
  write.csv(TMB.dat,"4.Demo of comparison of TMB among 5 identified subtype of breast cancer in TCGA-lusc cohort.csv")
  
  #5) compare fraction genome altered
  #接下来，compFGA（）不仅计算FGA，还计算每个亚型中每个样本的特定增益（FGG）或损失（FGL）。
  #为了实现此功能，应使用完全相同的列名来准备符合条件的分段副本编号输入，如下所示：      
  # change column names of segment data
  segment <- lusc.tcga$segment
  colnames(segment) <- c("sample","chrom","start","end","value")
  #Let’s see how the input looks like:
  
  head(segment)
  # sample chrom     start       end   value
  # 1 TCGA-3L-AA1B-01A     1   3218610 154074261 -0.1220
  # 2 TCGA-3L-AA1B-01A     1 154080441 154081058 -1.3462
  # 3 TCGA-3L-AA1B-01A     1 154089408 247813706 -0.1003
  # 4 TCGA-3L-AA1B-01A     2    484222  73318303  0.1054
  # 5 TCGA-3L-AA1B-01A     2  73321348  73330514 -0.9064
  # 6 TCGA-3L-AA1B-01A     2  73332815 242476062  0.0987
  #   Then this function can be run without any difficulties. Notably, if your CNA calling procedure did not provide a segmented copy number as value column but the original copy number, argument of iscopynumber must be switched to TRUE instead.
  #那么这个功能就可以毫无困难地运行了。值得注意的是，如果您的CNA调用过程没有提供分段副本编号作为值列，
  #而是提供原始副本编号，则iscopynumber的参数必须切换为TRUE。   
  # compare FGA, FGG, and FGL
  fga.lusc <- compFGA(moic.res     = cmoic.lusc_TLS,
                      segment      = segment,
                      iscopynumber = FALSE, # this is a segmented copy number file
                      cnathreshold = 0.2, # threshold to determine CNA gain or loss
                      test.method  = "parametric", # statistical testing method
                      fig.name     = "5.BARPLOT OF FGA")
  #> --2 samples mismatched from current subtypes.
  #> 5% 10% 15% 21% 26% 31% 36% 41% 46% 51% 57% 62% 67% 72% 77% 82% 88% 93% 98%
  head(fga.lusc$summary)
  summary <- fga.lusc$summary
  write.csv(summary,"5.Demo of comparison of fraction genome altered among 5 identified subtype of breast cancer in TCGA-lusc cohort.csv")
  #6) compare drug sensitivity
  #对于治疗指数较窄的药物，例如化疗药物，预测对药物的反应尤其重要，因为反应是高度可变的，副作用可能
  #是致命的。因此，Paul Geeleher等人（2014）10使用来源于细胞系的基线基因表达和体外药物敏感性，
  #结合体内基线肿瘤基因表达，来预测患者对药物的反应。Paul开发了一种R包pRRophetic，用于从肿瘤基因表达水
  #平预测临床化疗反应11，现在该功能已被用于MOVICS，以检查不同亚型之间的药物敏感性差异。
  #在此，我估计了在TCGA队列中，顺铂和紫杉醇对5种已确定的癌症亚型的IC50。
  # drug sensitivity comparison
  fpkm<-read.csv("fpkm_final.csv",header = T,row.names = 1)
  drug.coad <- compDrugsen(moic.res    = cmoic.lusc_TLS,
                           norm.expr   = fpkm[,cmoic.lusc_TLS$clust.res$samID], # double guarantee sample order
                           drugs       = c(  "Erlotinib",        #####引号里面的药物不要有空格
                                             "Rapamycin",
                                             "Sunitinib",
                                             # "PHA-665752",
                                             #"MG-132",
                                             "Paclitaxel",
                                             # "Cyclopamine",
                                             # "AZ628",
                                             "Sorafenib",
                                             # "VX-680",
                                             "Imatinib",
                                             "Cisplatin",
                                             # "TAE684",
                                             "Crizotinib",
                                             # "Saracatinib",
                                             "S-Trityl-L-cysteine",
                                             #                 "Z-LLNle-CHO",
                                             #                 "Dasatinib",
                                             #                 "GNF-2",
                                             #                 "CGP-60474",
                                             #                 "CGP-082996",
                                             #                 "A-770041",
                                             #                 "WH-4-023",
                                             #                 "WZ-1-84",
                                             #                 "BI-2536",
                                             #                 "BMS-536924",
                                             #                 "BMS-509744",
                                             "CMK",
                                             #                 "Pyrimethamine",
                                             #                 "JW-7-52-1",
                                             #                 "A-443654",
                                             #                 "GW843682X",
                                             #                 "MS-275",
                                             #                 "Parthenolide",
                                             #                 "KIN001-135",
                                             #                 "TGX221",
                                             #                 "Bortezomib",
                                             #                 "XMD8-85",
                                             #                 "Roscovitine",
                                             #                 "Salubrinal",
                                             #                 "Lapatinib",
                                             #                 "GSK269962A",
                                             "Doxorubicin",
                                             #                 "Etoposide",
                                             "Gemcitabine",#
                                             # "MitomycinC",
                                             # "Vinorelbine",
                                             # "NSC-87877",
                                             # "Bicalutamide",
                                             # "QS11",
                                             # "CP466722",
                                             "Midostaurin",
                                             # "CHIR-99021",
                                             # "AP-24534",
                                             # "AZD6482",
                                             # "JNK-9L",
                                             # "PF-562271",
                                             # "HG-6-64-1",
                                             # "JQ1",
                                             # "JQ12",
                                             "DMOG",
                                             # "FTI-277",
                                             # "OSU-03012",
                                             "Shikonin",
                                             "Embelin",
                                             # "FH535",
                                             # "PAC-1",
                                             # "IPA-3",
                                             # "GSK-650394",
                                             # "5-Fluorouracil",
                                             "Thapsigargin",
                                             # "BMS-754807",
                                             # "Lisitinib",
                                             "Bexarotene"
                                             # "Bleomycin",
                                             # "LFM-A13",
                                             # "GW-2580",
                                             # "AUY922",
                                             # "Phenformin",
                                             # "Pazopanib",
                                             # "LAQ824",
                                             # "GSK1904529A",
                                             # "BMS345541",
                                             # "Tipifarnib",
                                             # "BMS-708163",
                                             # "Ruxolitinib",
                                             # "AS601245",
                                             # "TL-2-105",
                                             # "AT-7519",
                                             # "TAK-715",
                                             # "BX-912",
                                             # "ZSTK474",
                                             # "AS605240",
                                             # "GSK1070916",
                                             # "KIN001-102",
                                             # "LY317615",
                                             # "GSK429286A",
                                             # "FMK",
                                             # "QL-XII-47",
                                             # "CAL-101",
                                             # "UNC0638",
                                             # "XL-184",
                                             # "WZ3105",
                                             # "XMD14-99",
                                             # "AC220",
                                             # "CP724714",
                                             # "JW-7-24-1",
                                             # "NPK76-II-72-1",
                                             # "STF-62247",
                                             # "NG-25",
                                             # "TL-1-85",
                                             # "VX-11e",
                                             #  "FR-180204",
                                             #  "Zibotentan",
                                             # "YM155",
                                             # "NSC-207895",
                                             # "AR-42",
                                             # "CUDC-101",
                                             # "Belinostat",
                                             # "I-BET-762",
                                             # "CAY10603",
                                             # "BIX02189",
                                             # "CH5424802",
                                             # "EKB-569",
                                             # "GSK2126458",
                                             # "KIN001-236",
                                             #  "KIN001-244",
                                             # "KIN001-055",
                                             # "KIN001-260",
                                             # "KIN001-266",
                                             # "Masitinib",
                                             # "MP470",
                                             # "MPS-1-IN-1",
                                             # "BHG712",
                                             # "OSI-930",
                                             # "OSI-027",
                                             # "CX-5461",
                                             # "PHA-793887",
                                             #  "PI-103",
                                             # "PIK-93",
                                             # "SB52334",
                                             # "TPCA-1",
                                             # "TG101348",
                                             # "Foretinib",
                                             # "Y-39983",
                                             #  "YM201636",
                                             #  "Tivozanib",
                                             # "GSK690693",
                                             # "SNX-2112",
                                             #  "QL-XI-92",
                                             # "XMD13-2",
                                             # "QL-X-138",
                                             # "XMD15-27"
                           ), # a vector of names of drug in GDSC
                           tissueType  = "digestive_system", # choose specific tissue type to construct ridge regression model
                           test.method = "nonparametric", # statistical testing method
                           prefix      = "6.BOXVIOLIN OF ESTIMATED IC50") 
  
  drug.coad <- compDrugsen(moic.res    = cmoic.lusc_TLS,
                           norm.expr   = fpkm[,cmoic.lusc_TLS$clust.res$samID], # double guarantee sample order
                           drugs       = c(  "5-Fluorouracil",
                                             "AC220",
                                             "BMS-536924",
                                             "BMS-708163",
                                             #"BMS-754806",
                                             "CH5424802",
                                             "Cyclopamine",
                                             "EKB-569",
                                             "GSK690693",
                                             "GW843682X",
                                             "LY317615",
                                             "Lapatinib",
                                             "Lisitinib",
                                             "MG-132",
                                             "MPS-1-IN-1",
                                             "Masitinib",
                                             "Phenformin",
                                             "Roscovitine",
                                             "Sorafenib",
                                             "TGX221",
                                             "VX-11e",
                                             "XMD8-85"
                           ), # a vector of names of drug in GDSC
                           tissueType  = "digestive_system", # choose specific tissue type to construct ridge regression model
                           test.method = "nonparametric", # statistical testing method
                           prefix      = "6.BOXVIOLIN OF ESTIMATED IC50_high") 
  
  ########将药物敏感性的结果变为dataframe,其中行名中有药物信息
  lala <- do.call(rbind, lapply(drug.coad, as.data.frame))
  #save(drug.coad,file = "drug.coad.Rdata")
  #write.csv(drugData2016,"drugData2016.csv")
  ######将行名中的药物提取出来
  lala2 <- lala %>% 
    rownames_to_column(var = "rowname") %>% 
    separate(rowname, into = c("drug", "sample"))
  #####先变为宽格式
  lala_wider <- tidyr::pivot_wider(lala2, 
                                   names_from=drug,
                                   values_from =Est.IC50)
  #####加入index信息
  rownames(lala_wider)<-lala_wider$sample
  lala_wider    <- lala_wider[rownames(ydata.risk), ]
  lala_wider$index<-ydata.risk$index
  
  #####变为长格式
  colnames(lala_wider)
  
  lala_wider_long <- 
    tidyr::pivot_longer(lala_wider,
                        cols = c(
                          Erlotinib,
                          Rapamycin,
                          Sunitinib,
                          Paclitaxel,
                          Sorafenib,
                          Imatinib,
                          Cisplatin,
                          Crizotinib,
                          `S-Trityl-L-cysteine`,
                          CMK,
                          Doxorubicin,
                          Gemcitabine,
                          Midostaurin,
                          DMOG,
                          Shikonin,
                          Embelin,
                          Thapsigargin,
                          Bexarotene ), names_to = "drug", values_to = "value")
  ########做散点图
  cell <- "drug"
  pdf_file <- paste(cell ,"index_药物敏感性_lusc-pointPlot_high.pdf", sep="_")
  pdf(pdf_file,width=6.5,height= 12.5)
  
  ggplot(lala_wider_long , aes(x=log2(index),y=value))+
    geom_point(size=0.8,colour ="#E42E50" )+
    geom_smooth(method = "lm",se = T,colour="#E42E50",fill ="#E42E50")+
    stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(p.value.label),sep = "~~~")),formula = y~x,parse = TRUE, size=3,colour ="black")+
    stat_correlation(mapping = use_label("R"),method = "pearson",label.x="left",npcy = 0.8,size=3,colour ="#588b8b")+
    labs(x=" ",y=NULL)+
    facet_wrap(.~drug,scales="free",ncol=3)+
    scale_y_continuous(expand = c(0, 0))+
    theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
          panel.grid =element_blank(),
          axis.text = element_text(colour="black"),legend.key = element_blank(),
          legend.title=element_blank())+
    labs(fill=' ')
  dev.off()
  
 write.csv(lala_wider,"药物敏感性分数.csv") 
  '''''
#—————药物敏感性———单独计算————————————————————————————————————————————————————————————————————————————————————————————

calcPhenotype(trainingExprData = CTRP2_Expr,
              trainingPtype = CTRP2_Res,
              testExprData = as.matrix(fpkm),#需要matrix
              batchCorrect = 'eb',  
              #IC50是对数转换的，所以表达矩阵也用对数转换过的
              powerTransformPhenotype = F,
              minNumSamples = 20,
              printOutput = T,
              removeLowVaryingGenes = 0.2,
              removeLowVaringGenesFrom = "homogenizeData"
)
'''''


#————————————————————————————GSE验证模型——————————————————————————————————
library(survival)
library(survminer)
library(dplyr)
####Scale数据
xdata.tmp<-read.csv("GSE14814-建模基因.csv",header = T,row.names = 1)
ydata.tmp<-read.csv('GSE14814_clin.info.csv',row.names = 1,header = T)
##转换为数据框
#xdata.tmp <- xdata.Scale.raw %>% as.data.frame()
##用得到的系数计算指数
index <- #+xdata.tmp$AKR1E2*0.229198473259477+
  xdata.tmp$CA4*0.148312217553774+
  xdata.tmp$CASP5*0.00637410093834286 +
  #xdata.tmp$CCN6 *0.173900788432551 -
  xdata.tmp$CDO1*0.11839886255889 -
  xdata.tmp$CLTCL1 *0.00694130784369265 +
  xdata.tmp$CPS1*0.0489633240878769 -
  xdata.tmp$CTSV  *0.043331201443875 -
  xdata.tmp$DMD *0.0440632233178641 +
  xdata.tmp$DUSP4*0.0616708981957713 +
  xdata.tmp$FGA *0.0614082707785454 +
  #xdata.tmp$H3C8*0.07201920954824 +
  xdata.tmp$KREMEN2 *0.0490970135208903 +
  #xdata.tmp$LOC102723996*0.0643937778088998 +
  #xdata.tmp$LRRC55 *0.169797408085259+
  #xdata.tmp$MIR3142HG*0.0203017920106524 -
  xdata.tmp$PARM1 *0.0181819616359712 -
  #xdata.tmp$PRR15*0.00137821893635317 +
  #xdata.tmp$SLC49A3 *0.0399429253527924 +
  xdata.tmp$SULT1C2*0.11429977464504 
##把生成的指数这一列加到原来的表格中
xdata.tmp$index <- index
##根据index排序,并分为2群


#根据中位值，把样品分为两组
group=ifelse(xdata.tmp[,13]>median(xdata.tmp[,13]),"High_risk","Low_risk")
diff=survdiff(Surv(time, status) ~group,data = ydata.tmp)
pValue=1-pchisq(diff$chisq,df=1)

pValue=paste0("p=",sprintf("%.02f",pValue))   #保留小数点后几位

#ydata.raw<-read.csv("PACA-AC-clin.info.csv",header = T,row.names = 1)
fit <- survfit(Surv(time, status) ~ group, data = ydata.tmp) #剩余数目

#绘制生存分析图
#dev.new()
hh=ggsurvplot(fit, 
              data= ydata.tmp,
              conf.int=F,  #置信区间
              pval=pValue,
              pval.size=5,
              legend.labs=c("High_risk","Low_risk"),
              legend.title= "GSE14814",
              xlab="Time(Days)",
              break.time.by = 10,
              risk.table.title="",
              risk.table=F,
              risk.table.height=.25,
              palette =c("#ff8216","#018a89"),
              add.all = F)
pdf(file="GSE14814_surivalPlot.pdf",width = 5,height =5)
print(hh)
dev.off()



write.csv(xdata.tmp,'exp_seq.GSE157010_index.tsv')
#————————————————————————————GSE14814验证模型——————————————————————————————————
########确定好模型后，分别在Scale和Nonescale数据表中加入risk分群列，然后保存分群后的数据信息
# ABI3 DOK3 GYPC LAT2 RASGRP3 TBC1D10C
library(survival)
library(survminer)
library(dplyr)
####Scale数据
xdata.tmp<-read.csv("GSE14814-建模基因.csv",header = T,row.names = 1)
ydata.tmp<-read.csv('GSE14814_clin.info.csv',row.names = 1,header = T)
##转换为数据框
#xdata.tmp <- xdata.Scale.raw %>% as.data.frame()
##用得到的系数计算指数
index <- #+xdata.tmp$AKR1E2*0.229198473259477+
  xdata.tmp$CA4*0.148312217553774+
  xdata.tmp$CASP5*0.00637410093834286 +
  #xdata.tmp$CCN6 *0.173900788432551 -
  xdata.tmp$CDO1*0.11839886255889 -
  xdata.tmp$CLTCL1 *0.00694130784369265 +
  xdata.tmp$CPS1*0.0489633240878769 -
  xdata.tmp$CTSV  *0.043331201443875 -
  xdata.tmp$DMD *0.0440632233178641 +
  xdata.tmp$DUSP4*0.0616708981957713 +
  xdata.tmp$FGA *0.0614082707785454 +
  #xdata.tmp$H3C8*0.07201920954824 +
  xdata.tmp$KREMEN2 *0.0490970135208903 +
  #xdata.tmp$LOC102723996*0.0643937778088998 +
  #xdata.tmp$LRRC55 *0.169797408085259+
  #xdata.tmp$MIR3142HG*0.0203017920106524 -
  xdata.tmp$PARM1 *0.0181819616359712 -
  #xdata.tmp$PRR15*0.00137821893635317 +
  #xdata.tmp$SLC49A3 *0.0399429253527924 +
  xdata.tmp$SULT1C2*0.11429977464504 
##把生成的指数这一列加到原来的表格中
xdata.tmp$index <- index
##根据index排序,并分为2群


#根据中位值，把样品分为两组
group=ifelse(xdata.tmp[,3]>median(xdata.tmp[,3]),"High_risk","Low_risk")
diff=survdiff(Surv(time, status) ~group,data = ydata.tmp)
pValue=1-pchisq(diff$chisq,df=1)

pValue=paste0("p=",sprintf("%.02f",pValue))   #保留小数点后几位

#ydata.raw<-read.csv("PACA-AC-clin.info.csv",header = T,row.names = 1)
fit <- survfit(Surv(time, status) ~ group, data = ydata.tmp) #剩余数目

#绘制生存分析图
#dev.new()
hh=ggsurvplot(fit, 
              data= ydata.tmp,
              conf.int=F,  #置信区间
              pval=pValue,
              pval.size=5,
              legend.labs=c("High_risk", "Low_risk"),
              legend.title= "GSE14814",
              xlab="Time(Days)",
              break.time.by = 10,
              risk.table.title="",
              risk.table=F,
              risk.table.height=.25,
              palette =c("#ff8216","#018a89"),
              add.all = F)
pdf(file="GSE14814_surivalPlot.pdf",width = 5,height =5)
print(hh)
dev.off()



write.csv(xdata.tmp,'exp_seq.GSE14814_index.tsv')


#————————————————————————————GSE73403，GSE157010，GSE19188、GSE18842、GSE37745验证模型——————————————————————————————————
########确定好模型后，分别在Scale和Nonescale数据表中加入risk分群列，然后保存分群后的数据信息
# ABI3 DOK3 GYPC LAT2 RASGRP3 TBC1D10C
library(survival)
library(survminer)
library(dplyr)
####Scale数据
xdata.tmp<-read.csv("GSE37745-建模基因.csv",header = T,row.names = 1)
ydata.tmp<-read.csv('GSE37745_clin.info.csv',row.names = 1,header = T)
##转换为数据框
#xdata.tmp <- xdata.Scale.raw %>% as.data.frame()
##用得到的系数计算指数
index <- +xdata.tmp$AKR1E2*0.229198473259477+
  xdata.tmp$CA4*0.148312217553774+
  xdata.tmp$CASP5*0.00637410093834286 +
  #xdata.tmp$CCN6 *0.173900788432551 -
  xdata.tmp$CDO1*0.11839886255889 -
  xdata.tmp$CLTCL1 *0.00694130784369265 +
  xdata.tmp$CPS1*0.0489633240878769 -
  xdata.tmp$CTSV  *0.043331201443875 -
  xdata.tmp$DMD *0.0440632233178641 +
  xdata.tmp$DUSP4*0.0616708981957713 +
  xdata.tmp$FGA *0.0614082707785454 +
  #xdata.tmp$H3C8*0.07201920954824 +
  xdata.tmp$KREMEN2 *0.0490970135208903 +
  #xdata.tmp$LOC102723996*0.0643937778088998 +
  xdata.tmp$LRRC55 *0.169797408085259+
  #xdata.tmp$MIR3142HG*0.0203017920106524 -
  xdata.tmp$PARM1 *0.0181819616359712 -
  xdata.tmp$PRR15*0.00137821893635317 +
  #xdata.tmp$SLC49A3 *0.0399429253527924 +
  xdata.tmp$SULT1C2*0.11429977464504 
##把生成的指数这一列加到原来的表格中
xdata.tmp$index <- index
##根据index排序,并分为2群


#根据中位值，把样品分为两组
group=ifelse(xdata.tmp[,16]>median(xdata.tmp[,16]),"High_risk","Low_risk")
diff=survdiff(Surv(time, status) ~group,data = ydata.tmp)
pValue=1-pchisq(diff$chisq,df=1)

pValue=paste0("p=",sprintf("%.02f",pValue))   #保留小数点后几位

#ydata.raw<-read.csv("PACA-AC-clin.info.csv",header = T,row.names = 1)
fit <- survfit(Surv(time, status) ~ group, data = ydata.tmp) #剩余数目

#绘制生存分析图
#dev.new()
hh=ggsurvplot(fit, 
              data= ydata.tmp,
              conf.int=F,  #置信区间
              pval=pValue,
              pval.size=5,
              #legend.labs=c("High_risk","Low_risk",),
              legend.title= "GSE37745",
              xlab="Time(Days)",
              break.time.by = 10,
              risk.table.title="",
              risk.table=F,
              risk.table.height=.25,
              #palette =c("#ff8216","#018a89"),
              add.all = F)
pdf(file="GSE37745_surivalPlot.pdf",width = 5,height =5)
print(hh)
dev.off()


#write.csv(xdata.tmp,'exp_seq.GSE73403_index.tsv')



#————————————————————————————GSE33532验证模型——————————————————————————————————
########确定好模型后，分别在Scale和Nonescale数据表中加入risk分群列，然后保存分群后的数据信息
# ABI3 DOK3 GYPC LAT2 RASGRP3 TBC1D10C
library(survival)
library(survminer)
library(dplyr)
####Scale数据
xdata.tmp<-read.csv("GSE33532-建模基因.csv",header = T,row.names = 1)
ydata.tmp<-read.csv('GSE33532_clin.info.csv',row.names = 1,header = T)
##转换为数据框
#xdata.tmp <- xdata.Scale.raw %>% as.data.frame()
##用得到的系数计算指数
index <- +xdata.tmp$AKR1E2*0.229198473259477+
  xdata.tmp$CA4*0.148312217553774+
  xdata.tmp$CASP5*0.00637410093834286 +
  #xdata.tmp$CCN6 *0.173900788432551 -
  xdata.tmp$CDO1*0.11839886255889 -
  xdata.tmp$CLTCL1 *0.00694130784369265 +
  xdata.tmp$CPS1*0.0489633240878769 -
  #xdata.tmp$CTSV  *0.043331201443875 -
  xdata.tmp$DMD *0.0440632233178641 +
  xdata.tmp$DUSP4*0.0616708981957713 +
  xdata.tmp$FGA *0.0614082707785454 +
  #xdata.tmp$H3C8*0.07201920954824 +
  xdata.tmp$KREMEN2 *0.0490970135208903 +
  #xdata.tmp$LOC102723996*0.0643937778088998 +
  xdata.tmp$LRRC55 *0.169797408085259+
  #xdata.tmp$MIR3142HG*0.0203017920106524 -
  xdata.tmp$PARM1 *0.0181819616359712 -
  xdata.tmp$PRR15*0.00137821893635317 +
  #xdata.tmp$SLC49A3 *0.0399429253527924 +
  xdata.tmp$SULT1C2*0.11429977464504 
##把生成的指数这一列加到原来的表格中
xdata.tmp$index <- index
##根据index排序,并分为2群


#根据中位值，把样品分为两组
group=ifelse(xdata.tmp[,3]>median(xdata.tmp[,3]),"High_risk","Low_risk")
diff=survdiff(Surv(time, status) ~group,data = ydata.tmp)
pValue=1-pchisq(diff$chisq,df=1)

pValue=paste0("p=",sprintf("%.02f",pValue))   #保留小数点后几位

#ydata.raw<-read.csv("PACA-AC-clin.info.csv",header = T,row.names = 1)
fit <- survfit(Surv(time, status) ~ group, data = ydata.tmp) #剩余数目

#绘制生存分析图
#dev.new()
hh=ggsurvplot(fit, 
              data= ydata.tmp,
              conf.int=F,  #置信区间
              pval=pValue,
              pval.size=5,
              legend.labs=c("High_risk", "Low_risk"),
              legend.title= "GSE33479",
              xlab="Time(Days)",
              break.time.by = 10,
              risk.table.title="",
              risk.table=F,
              risk.table.height=.25,
              palette =c("#ff8216","#018a89"),
              add.all = F)
pdf(file="GSE33479_surivalPlot.pdf",width = 5,height =5)
print(hh)
dev.off()


write.csv(xdata.tmp,'exp_seq.GSE33479_index.tsv')



###-----------------------样本随机分群------------------------------------------##
library(caret)

grouplist <- createDataPartition(y=xdata.Scale.raw[,2] , p = 0.6 , list = F)

train.x.data <- xdata.Scale.raw[grouplist,]
test.x.data <- xdata.Scale.raw[-grouplist,]

train.y.data <- ydata.raw[rownames(ydata.raw) %in% 
                            rownames(train.x.data),]

test.y.data <- ydata.raw[rownames(ydata.raw) %in% 
                           rownames(test.x.data),]

test.x.data<-t(test.x.data)%>%as.data.frame()
train.x.data<-t(train.x.data)%>%as.data.frame()

write.csv(train.y.data,"train.y.data.csv")
write.csv(test.y.data,"test.y.data.csv")
write.csv(train.x.data,"train.x.data.csv")
write.csv(test.x.data,"test.x.data.csv")
#—————————————————————————train———TCGA数据集-验证模型————————————————用了它——————————————————
########确定好模型后，分别在Scale和Nonescale数据表中加入risk分群列，然后保存分群后的数据信息
# ABI3 DOK3 GYPC LAT2 RASGRP3 TBC1D10C
library(survival)
library(survminer)
library(dplyr)
####Scale数据

xdata.tmp<-read.csv("train.xdata.modle.csv",row.names = 1,header = T)
ydata.tmp<-read.csv("train.y.data.csv",row.names = 1,header = T)
##转换为数据框
#xdata.tmp <- xdata.Scale.raw %>% as.data.frame()
##用得到的系数计算指数
index <- +xdata.tmp$AKR1E2*0.229198473259477+
  xdata.tmp$CA4*0.148312217553774+
  xdata.tmp$CASP5*0.00637410093834286 +
  xdata.tmp$CCN6 *0.173900788432551 -
  xdata.tmp$CDO1*0.11839886255889 -
  xdata.tmp$CLTCL1 *0.00694130784369265 +
  xdata.tmp$CPS1*0.0489633240878769 -
  xdata.tmp$CTSV  *0.043331201443875 -
  xdata.tmp$DMD *0.0440632233178641 +
  xdata.tmp$DUSP4*0.0616708981957713 +
  xdata.tmp$FGA *0.0614082707785454 +
  xdata.tmp$H3C8*0.07201920954824 +
  xdata.tmp$KREMEN2 *0.0490970135208903 +
  xdata.tmp$LOC102723996*0.0643937778088998 +
  xdata.tmp$LRRC55 *0.169797408085259+
  xdata.tmp$MIR3142HG*0.0203017920106524 -
  xdata.tmp$PARM1 *0.0181819616359712 -
  xdata.tmp$PRR15*0.00137821893635317 +
  xdata.tmp$SLC49A3 *0.0399429253527924 +
  xdata.tmp$SULT1C2*0.11429977464504 
##把生成的指数这一列加到原来的表格中
xdata.tmp$index <- index
##根据index排序,并分为2群


#根据中位值，把样品分为两组
group=ifelse(xdata.tmp[,21]>median(xdata.tmp[,21]),"High_risk","Low_risk")
diff=survdiff(Surv(time, status) ~group,data = ydata.tmp)
pValue=1-pchisq(diff$chisq,df=1)

pValue=paste0("p=2.25029799383325e-06")   #保留小数点后几位

#ydata.raw<-read.csv("PACA-AC-clin.info.csv",header = T,row.names = 1)
fit <- survfit(Surv(time, status) ~ group, data = ydata.tmp) #剩余数目

#绘制生存分析图
#dev.new()
hh=ggsurvplot(fit, 
              data= ydata.tmp,
              conf.int=F,  #置信区间
              pval=pValue,
              pval.size=5,
              legend.labs=c("High_risk","Low_risk"),
              legend.title= "train",
              xlab="Time(Days)",
              break.time.by = 10,
              risk.table.title="",
              risk.table=F,
              risk.table.height=.25,
              palette =c("#ff8216","#018a89"),
              add.all = F)
pdf(file="train_surivalPlot.pdf",width = 5,height =5)
print(hh)
dev.off()


write.csv(xdata.tmp,'exp_seq.GSE33479_index.tsv')

#—————————————————————————test———TCGA数据集-验证模型——————————————————————————————————
########确定好模型后，分别在Scale和Nonescale数据表中加入risk分群列，然后保存分群后的数据信息
# ABI3 DOK3 GYPC LAT2 RASGRP3 TBC1D10C
library(survival)
library(survminer)
library(dplyr)
####Scale数据

xdata.tmp<-read.csv("test.xdata.modle.csv",row.names = 1,header = T)
ydata.tmp<-read.csv("test.y.data.csv",row.names = 1,header = T)
##转换为数据框
#xdata.tmp <- xdata.Scale.raw %>% as.data.frame()
##用得到的系数计算指数
index <- +xdata.tmp$AKR1E2*0.229198473259477+
  xdata.tmp$CA4*0.148312217553774+
  xdata.tmp$CASP5*0.00637410093834286 +
  xdata.tmp$CCN6 *0.173900788432551 -
  xdata.tmp$CDO1*0.11839886255889 -
  xdata.tmp$CLTCL1 *0.00694130784369265 +
  xdata.tmp$CPS1*0.0489633240878769 -
  xdata.tmp$CTSV  *0.043331201443875 -
  xdata.tmp$DMD *0.0440632233178641 +
  xdata.tmp$DUSP4*0.0616708981957713 +
  xdata.tmp$FGA *0.0614082707785454 +
  xdata.tmp$H3C8*0.07201920954824 +
  xdata.tmp$KREMEN2 *0.0490970135208903 +
  xdata.tmp$LOC102723996*0.0643937778088998 +
  xdata.tmp$LRRC55 *0.169797408085259+
  xdata.tmp$MIR3142HG*0.0203017920106524 -
  xdata.tmp$PARM1 *0.0181819616359712 -
  xdata.tmp$PRR15*0.00137821893635317 +
  xdata.tmp$SLC49A3 *0.0399429253527924 +
  xdata.tmp$SULT1C2*0.11429977464504 
##把生成的指数这一列加到原来的表格中
xdata.tmp$index <- index
##根据index排序,并分为2群


#根据中位值，把样品分为两组
group=ifelse(xdata.tmp[,21]>median(xdata.tmp[,21]),"High_risk","Low_risk")
diff=survdiff(Surv(time, status) ~group,data = ydata.tmp)
pValue=1-pchisq(diff$chisq,df=1)

pValue=paste0("p=0.00385")   #保留小数点后几位

#ydata.raw<-read.csv("PACA-AC-clin.info.csv",header = T,row.names = 1)
fit <- survfit(Surv(time, status) ~ group, data = ydata.tmp) #剩余数目

#绘制生存分析图
#dev.new()
hh=ggsurvplot(fit, 
              data= ydata.tmp,
              conf.int=F,  #置信区间
              pval=pValue,
              pval.size=5,
              legend.labs=c("High_risk","Low_risk"),
              legend.title= "test",
              xlab="Time(Days)",
              break.time.by = 10,
              risk.table.title="",
              risk.table=F,
              risk.table.height=.25,
              palette =c("#ff8216","#018a89"),
              add.all = F)
pdf(file="test_surivalPlot.pdf",width = 5,height =5)
print(hh)
dev.off()



#————————————————————————————TCGA数据集-验证模型——————————————————————————————————
########确定好模型后，分别在Scale和Nonescale数据表中加入risk分群列，然后保存分群后的数据信息
# ABI3 DOK3 GYPC LAT2 RASGRP3 TBC1D10C
library(survival)
library(survminer)
library(dplyr)
####Scale数据
xdata.tmp<-xdata.ll%>%as.data.frame()
xdata.tmp<-train.x.data
ydata.tmp<-ydata.ll
##转换为数据框
#xdata.tmp <- xdata.Scale.raw %>% as.data.frame()
##用得到的系数计算指数
index <- +xdata.tmp$AKR1E2*0.229198473259477+
  xdata.tmp$CA4*0.148312217553774+
  xdata.tmp$CASP5*0.00637410093834286 +
  xdata.tmp$CCN6 *0.173900788432551 -
  xdata.tmp$CDO1*0.11839886255889 -
  xdata.tmp$CLTCL1 *0.00694130784369265 +
  xdata.tmp$CPS1*0.0489633240878769 -
  xdata.tmp$CTSV  *0.043331201443875 -
  xdata.tmp$DMD *0.0440632233178641 +
  xdata.tmp$DUSP4*0.0616708981957713 +
  xdata.tmp$FGA *0.0614082707785454 +
  xdata.tmp$H3C8*0.07201920954824 +
  xdata.tmp$KREMEN2 *0.0490970135208903 +
  xdata.tmp$LOC102723996*0.0643937778088998 +
  xdata.tmp$LRRC55 *0.169797408085259+
  xdata.tmp$MIR3142HG*0.0203017920106524 -
  xdata.tmp$PARM1 *0.0181819616359712 -
  xdata.tmp$PRR15*0.00137821893635317 +
  xdata.tmp$SLC49A3 *0.0399429253527924 +
  xdata.tmp$SULT1C2*0.11429977464504 
##把生成的指数这一列加到原来的表格中
xdata.tmp$index <- index
##根据index排序,并分为2群


#根据中位值，把样品分为两组
group=ifelse(xdata.tmp[,27]>median(xdata.tmp[,27]),"High_risk","Low_risk")
diff=survdiff(Surv(time, status) ~group,data = ydata.tmp)
pValue=1-pchisq(diff$chisq,df=1)

pValue=paste0("p=",sprintf("%.02f",pValue))   #保留小数点后几位

#ydata.raw<-read.csv("PACA-AC-clin.info.csv",header = T,row.names = 1)
fit <- survfit(Surv(time, status) ~ group, data = ydata.tmp) #剩余数目

#绘制生存分析图
#dev.new()
hh=ggsurvplot(fit, 
              data= ydata.tmp,
              conf.int=F,  #置信区间
              pval=pValue,
              pval.size=5,
              legend.labs=c( "High_risk","Low_risk"),
              legend.title= "建模数据",
              xlab="Time(Days)",
              break.time.by = 10,
              risk.table.title="",
              risk.table=F,
              risk.table.height=.25,
              palette =c("#ff8216","#018a89"),
              add.all = F)
pdf(file="建模数据_surivalPlot.pdf",width = 5,height =5)
print(hh)
dev.off()

#7）与其他亚型比较一致性
#目前，许多癌症都有传统的分类，评估新亚型与以前分类的一致性对于反映聚类分析的稳健性和确定潜在但新的亚型至关重要。
#为了测量当前子类型与其他预先存在的分类之间的一致性（相似性），MOVICS提供了compAgree（）
#的功能来计算四个统计数据：Rand指数（RI）12、调整后的相互信息（AMI）13、Jaccard指数（JI）14和Fowlkes Mallows指数（FM）15；
#所有这些测量值都在0到1之间，值越大，两个评估对象就越相似。该功能还可以生成冲积图，以可视化两个评估对象与当前子类型的一致性作为参考。
surv.info<-read.csv("surv.info.csv",header = T,row.names = 1)

# customize the factor level for pstage
surv.info$pstage <- factor(surv.info$PATH_T_STAGE, levels = c("T1","T1A","T1B","T2","T2A","T2B","TB","T3","T4"))
surv.info$SEX <- factor(surv.info$SEX, levels = c("Male","Female"))
surv.info$M_STAGE <- factor(surv.info$PATH_M_STAGE, levels = c("MX","M0","M1","M1X","M1B"))
surv.info$N_STAGE <- factor(surv.info$PATH_N_STAGE, levels = c("NX","N0","N1","N1A","N1B","N1C","N2","N2A","N2B"))

# agreement comparison (support up to 6 classifications include current subtype)
agree.coad <- compAgree(moic.res  = cmoic.lusc,
                        subt2comp = surv.info[,c("SEX","PATH_T_STAGE","PATH_M_STAGE","PATH_N_STAGE")],
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "7.AGREEMENT OF CONSENSUSMOIC WITH PAM50 AND PSTAGE")
#> --all samples matched.

print(agree.coad)
write.csv(agree.coad,"7.Agreement of 5 identified subtypes with PAM50 classification and pathological stage in TCGA-BRCA cohort.csv")

#四、 RUN Module——————————————————————————————————————————————————————————————————————————————————————————————————
#深呼吸，我们成功进入最后一个RUN模块，胜利在望。在本模块中，MOVICS旨在通过识别不同亚型的潜在预测
#生物标志物和功能途径来表征它们。识别和应用分子生物标志物有效预测亚型对疾病管理和治疗尤其重要，
#从而提高临床疗效。在这种情况下，MOVICS从差异表达分析（DEA）开始寻找亚型特异性生物标志物。

#1） 运行差分表达式分析

#MOVICS提供runDEA（）功能，嵌入了三种最先进的DEA方法来识别差异表达基因（DEG），
#包括用于RNA-Seq counts数据的edgeR和DESeq2，以及用于微阵列图谱或标准化表达数据的limma。
#值得注意的是，由于runDEA（）在选择limma算法时会自动检查数据量表，因此建议提供微阵列表达谱或
#标准化表达数据（例如，RSEM、FPKM、TPM），而无需z-score或log2转换。这一步骤也相当耗时，
#尤其是在使用原始计数数据的情况下，所以请将手从键盘上松开，闭上眼睛，放松片刻。

# run DEA with edgeR
runDEA(dea.method = "edger",      ##三选一：edger和deseq2，以及limma。
       expr       = count, # raw count data
       moic.res   = cmoic.lusc,
       prefix     = "TCGA-LUSC") # prefix of figure name
# --all samples matched.
# --you choose edger and please make sure an RNA-Seq count data was provided.
# edger of CS1_vs_Others done...
# edger of CS2_vs_Others done...

# run DEA with DESeq2
runDEA(dea.method = "deseq2",
       expr       = count,
       moic.res   = cmoic.lusc,
       prefix     = "TCGA-LUSC")
#> --all samples matched.
#> --you choose deseq2 and please make sure an RNA-Seq count data was provided.
#> deseq2 of CS1_vs_Others done...
#> deseq2 of CS2_vs_Others done...


# run DEA with limma
runDEA(dea.method = "limma",
       expr       = fpkm, # normalized expression data
       moic.res   =cmoic.lusc,
       prefix     = "TCGA-LUSC")
#> --all samples matched.
#> --you choose limma and please make sure a microarray profile or a normalized expression data [FPKM or TPM without log2 transformation is recommended] was provided.
#> --log2 transformation done for expression data.
#> limma of CS1_vs_Others done...
#> limma of CS2_vs_Others done...

#每个已识别的癌症亚型将与其他（其他）亚型进行比较，并根据res.path的参数存储相应的.txt文件。
#默认情况下，这些文件将保存在当前工作目录下。

#2) run biomarker identification procedure
#在该程序中，选择按log2FoldChange排序的差异表达最多的基因作为每个亚型的生物标志物
#（默认情况下，每个亚型有200个生物标志物）。这些生物标志物应通过显著性阈值
#（例如，标称P值<0.05和调整后P值<0.005），并且不得与其他亚型的任何生物标志物重叠。

# choose edgeR result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.lusc,
                       dea.method    = "limma",   ##三选一：edger和deseq2，以及limma。
                       prefix        = "TCGA-LUSC", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = fpkm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "2.UPREGULATED BIOMARKER HEATMAP_limma")
#> --all samples matched.
#> --log2 transformation done for expression data.
# check the upregulated biomarkers
head(marker.up$templates)
templates <- marker.up$templates
write.csv(templates,"2.Demo of subtype-specific upregulated biomarkers for 5 identified subtypes of breast cancer in TCGA-COAD cohort.csv")


#然后尝试limma的结果来鉴定亚型特异性下调的生物标志物。请随时尝试DESeq2。
# choose limma result to identify subtype-specific down-regulated biomarkers
marker.dn <- runMarker(moic.res      = cmoic.lusc,
                       dea.method    = "limma",        ##三选一：edger和deseq2，以及limma。
                       prefix        = "TCGA-LUSC",
                       dirct         = "down",
                       n.marker      = 100, # switch to 50
                       doplot        = TRUE,
                       norm.expr     = fpkm,
                       annCol        = annCol,
                       annColors     = annColors,
                       fig.name      = "2.DOWNREGULATED BIOMARKER HEATMAP_limma")
#> --all samples matched.
#> --log2 transformation done for expression data.
head(marker.dn$templates)
templates <- marker.dn$templates

write.csv(templates,"2.Demo of subtype-specific downregulated biomarkers for 5 identified subtypes of breast cancer in TCGA-COAD cohort.csv")

#3) run gene set enrichment analysis
#类似地，基于其相应的DEA结果对每个亚型运行GSEA，以识别亚型特异性功能途径。
#为此，我准备了一个基因集背景，该背景包括来自分子特征数据库（MSigDB，https://www.gsea-msigdb.org/gsea/msigdb/index.jsp).你可以下载其他感兴趣的背景资料供自己学习。

# MUST locate ABSOLUTE path of msigdb file
MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)
#"/Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library/MOVICS/extdata/c5.bp.v7.1.symbols.xls"

#同样，这些确定的特定途径应通过显著性阈值（例如，标称P值<0.05和调整后P值<0.025），
#并且不得与其他亚型的任何途径重叠。在具有亚型特异性途径后，检索途径内的基因，
#通过使用GSVA R包计算单个样本富集分数。随后，亚型特异性富集分数将由亚型内的平均值或中值
#（默认为平均值）表示，并将通过对角线热图进一步可视化。

# run GSEA to identify up-regulated GO pathways using results from edgeR
gsea.up <- runGSEA(moic.res     = cmoic.lusc,
                   dea.method   = "limma",   ##三选一：edger和deseq2，以及limma。
                   prefix       = "TCGA-LUSC", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = fpkm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "3.UPREGULATED PATHWAY HEATMAP_limma")
#> --all samples matched.
#> GSEA done...
#> --log2 transformation done for expression data.
#> Estimating GSVA scores for 50 gene sets.
#> Estimating ECDFs with Gaussian kernels


#   |===================================================================== |  98%
# |                                                                            
#   |======================================================================| 100%
#> gsva done...
#> heatmap done...

#Check some columns of GSEA results for the first subtype (CS1).#
print(gsea.up$gsea.list$CS1[1:6,3:6])

#Also check results of subtype-specific enrichment scores.

head(round(gsea.up$grouped.es,3))

#然后尝试来源于DESeq2的结果来鉴定亚型特异性下调途径。请随意尝试limma。

# run GSEA to identify down-regulated GO pathways using results from DESeq2
gsea.dn <- runGSEA(moic.res     = cmoic.lusc,
                   dea.method   = "limma",      ##三选一：edger和deseq2，以及limma。
                   prefix       = "TCGA-LUSC",
                   msigdb.path  = MSIGDB.FILE,
                   norm.expr    = fpkm,
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "gsva", # switch to ssgsea/gsva
                   norm.method  = "mean", # switch to median/mean
                   fig.name     = "3.DOWNREGULATED PATHWAY HEATMAP_limma") 
#> --all samples matched.
#> GSEA done...
#> --log2 transformation done for expression data.
#> Estimating ssGSEA scores for 50 gene sets.

#|======================================================================| 100%
#> ssgsea done...
#> heatmap done...

#4) run gene set variation analysis
#对于所有新定义的分子亚型，描述其通过基因集的不同特征验证的特征是至关重要的。
#MOVICS提供了一个简单的函数，该函数使用基因集变异分析来基于给定的感兴趣基因集列表计算
#每个亚型中每个样本的富集分数。首先，我们必须准备一份感兴趣的基因列表，保存为GMT格式。

# MUST locate ABSOLUTE path of gene set file
GSET.FILE <- 
  system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "c2.cp.kegg.v2022.1.Hs.symbols.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "h.all.v2022.1.Hs.symbols.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "h.all.v2022.1.Hs.symbols_select.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "REACTOME_metabolism.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "28_immune_cell_GeneSignatures.gmt", package = "MOVICS", mustWork = TRUE)
GSET.FILE <- 
  system.file("extdata", "CMS1.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "CMS3.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "CMS4.gmt", package = "MOVICS", mustWork = TRUE)


GSET.FILE <- 
  system.file("extdata", "TF.gmt", package = "MOVICS", mustWork = TRUE)
GSET.FILE <- 
  system.file("extdata", "TF_wx.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "immune_gene.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "chromatin.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "signature_metabolism.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "signature_metabolism_select.gmt", package = "MOVICS", mustWork = TRUE)

'''''
     ########将药物敏感性的结果变为dataframe,其中行名中有药物信息
     lala <- do.call(rbind, lapply(signature_metabolism, as.data.frame))
     #save(drug.coad,file = "drug.coad.Rdata")
     #write.csv(drugData2016,"drugData2016.csv")
     ######将行名中的药物提取出来
     lala2 <- lala %>% 
       rownames_to_column(var = "rowname") %>% 
       separate(rowname, into = c("dx", "sample"), sep = "\\.")
     head(lala2,1)
     #####先变为宽格式
     lala_wider <- tidyr::pivot_wider(lala2, 
                                      names_from=dx,
 
                                      
   lala_wider_t<-t(lala_wider)                                   
        colnames(lala_wider_t)  =lala_wider_t[1,]  
        
        lala_wider_t2 <- lala_wider_t[-1,]                             
                                      
   write.table(lala_wider_t2,"signature_metabolism.txt",col.names = F,sep="\t",quote = F) 
 '''''    

#然后我们可以根据指定的方法和给定的感兴趣的基因集来计算单个样本的富集分数。
# run GSVA to estimate single sample enrichment score based on given gene set of interest
gsva.res <- 
  runGSVA(moic.res      = cmoic.lusc,
          norm.expr     = fpkm,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "ssgsea", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          fig.path      = getwd(),
          fig.name      = "4.GENE SETS OF INTEREST HEATMAP_hallmarker_TLS",
          height        = 7,
          width         = 8)

# check raw enrichment score
print(gsva.res$raw.es[1:3,1:3])

signature_metabolism114 <-gsva.res$raw.es
write.csv(signature_metabolism114 ,"4.GENE SETS OF INTEREST HEATMAP_hallmarker_TLS.csv")

# check z-scored and truncated enrichment score

#print(gsva.res$scaled.es[1:3,1:3])



##整理GEO数据
clin.info<-read.csv("GSE30219_clin.info.csv",row.names = 1,header = T)
mRNA.expr<-read.csv("GSE30219_xdata.csv",header = T,row.names = 1)
lusc.geo <- list(mRNA.expr = mRNA.expr,
                  clin.info=clin.info)

save(lusc.geo, file = "LUSC.geo.验证GEO30219.RData")  ##这部分是已经筛选了TLS相关基因作为分群数据的数据

#四、 RUN Module——————————————geo验证数据————————————————————————————————————————————————————————————————————————————————————

#NTP预测外部数据
# run NTP in Yau cohort by using up-regulated biomarkers
yau.ntp.pred <- runNTP(expr       = lusc.geo$mRNA.expr,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR YAU") 

# --original template has 400 biomarkers and 291 are matched in external expression profile.
# cosine correlation distance
# 278 samples; 4 classes; 33-94 features/class
# serial processing; 1000 permutation(s)...
# predicted samples/class (FDR<0.05)
# 
# CS1  CS2  CS3  CS4 <NA> 
#   37   32   66   66   77 

head(yau.ntp.pred$ntp.res)

ntp.res  <-yau.ntp.pred$ntp.res
write.csv(ntp.res,"5.Demo of predicted subtypes in Yau cohort by NTP using subtype-specific upregulated biomarkers identified from TCGA-BRCA cohort.csv")
#上面准备了一个对象yau.ntp.pred，该对象在结构上与getMOIC（）返回的对象相似，
#但只存储cluster.res，如果有额外的数据可用，这些cluster.res可以传递给COMP模块内的函数。
#例如，我在此首先比较Yau队列中预测的5种癌症亚型的生存结果。

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = lusc.geo$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "h", # switch to both
                     fig.name         = "5.KAPLAN-MEIER CURVE OF NTP FOR YAU",
                     #p.adjust.method='BY'   #c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none' 
) 
#> --a total of 682 samples are identified.
#> --cut survival curve up to 10 years.
print(surv.yau)
# $fitd
# Call:
#   survdiff(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
#            na.action = na.exclude)
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# Subtype=CS1 66       49     38.9    2.6066    3.3433
# Subtype=CS2 56       27     47.0    8.5160   11.5614
# Subtype=CS3 81       52     53.6    0.0491    0.0698
# Subtype=CS4 75       60     48.4    2.7573    3.7742
# 
# Chisq= 14.2  on 3 degrees of freedom, p= 0.003 
# 
# $fit
# Call: survfit(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
#               na.action = na.exclude, error = "greenwood", type = "kaplan-meier", 
#               conf.type = "plain")
# 
# n events median 0.95LCL 0.95UCL
# CS1 66     49   27.5    13.8    58.0
# CS2 56     27  124.9    83.6      NA
# CS3 81     52   62.0    39.3   102.3
# CS4 75     60   44.3    31.5    74.8
# 
# $xyrs.est
# [1] "[Not Available]: argument of xyrs.est was not specified."
# 
# $overall.p
# [1] 0.002678778
# 
# $pairwise.p
# 
# Pairwise comparisons using Log-Rank test 
# 
# data:  mosurv.res and Subtype 
# 
# CS1    CS2    CS3   
# CS2 0.0036 -      -     
#   CS3 0.2285 0.0430 -     
#   CS4 0.9119 0.0025 0.2285
# 
# P value adjustment method: BH 

# predict subtype in discovery cohort using PAM
tcga.pam.pred <- runPAM(train.expr  = fpkm,
                        moic.res    = tcga.ntp.pred,
                        test.expr   = fpkm)
# --all samples matched.
# --a total of 21654 genes shared and used.
# --training expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.
# --testing expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.


# check consistency between current and NTP-predicted subtype in discovery TCGA-BRCA
runKappa(subt1     = tcga.pam.pred$clust.res$clust,
         subt2     = tcga.ntp.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "NTP",
         fig.name  = "7.CONSISTENCY HEATMAP FOR TCGA between CMOIC and NTP")
# check consistency between current and PAM-predicted subtype in discovery TCGA-BRCA
runKappa(subt1     = tcga.ntp.pred$clust.res$clust,
         subt2     = tcga.pam.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "PAM",
         fig.name  = "7.CONSISTENCY HEATMAP FOR TCGA between CMOIC and PAM")

# check consistency between NTP and PAM-predicted subtype in validation Yau-BRCA
runKappa(subt1     = tcga.ntp.pred$clust.res$clust,
         subt2     = tcga.pam.pred$clust.res$clust,
         subt1.lab = "NTP",
         subt2.lab = "PAM",
         fig.name  = "7.CONSISTENCY HEATMAP FOR YAU")