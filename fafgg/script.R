### 1 ####
library("tidyr")

setwd('E:/GuangZhou/Job/WG/Data')

filelist = list.files(path = "./RNA_seq/Science/GSE57299_RAW",pattern = "*.gz")
dir = "./RNA_seq/Science/GSE57299_RAW"
library(dplyr)
library(data.table)
da = data.frame()
attach(frame)
for (file in filelist){ 
  fiedir = paste(dir,file,sep = "/")
  d = read.table(gzfile(fiedir),header = F,sep = "\t")
  to = strsplit(file, "_", fixed= T)
  d$file = paste(to[[1]][3],to[[1]][4],sep = "_")
  


  #counid_d = merge(da,d,by.x = "V4",by.y = "V4",all.x = FALSE,all.y=FALSE)

  counid_d = d[!duplicated(d[,4]),]
  da = rbind(da,counid_d[,c(4:5,7)])
  rm(counid_d,d)
  #########
  #da = rbind(da,counid_d[,c(4:5,7)])
}

##########  web sample 1 
# mydata<-data.frame(
#   name=c("store1","store2","store3","store4"),
#   address=c("普陀区","黄浦区","徐汇区","浦东新区"),
#   sale2014=c(3000,2500,2100,1000),
#   sale2015=c(3020,2800,3900,2000),
#   sale2016=c(5150,3600,2700,2500),
#   sale2017=c(4450,4100,4000,3200)
# )
# mydata1<-melt(
#   mydata,
#   id.vars=c("address","name"),#要保留的主字段
#   variable.name = "Year",#转换后的分类字段名称（维度）
#   value.name = "Sale" #转换后的度量值名称
# )
# 
# kk = spread(
#   data=mydata1,
#   key=Year,
#   value=Sale
# )

################

#长转宽——dcast

#mk = pivot_wider(da,names_from = file)

kl = spread(
  data=da,
  key=file,
  value=V5
)


paper = c("TIPARP","GPR68","AHRR","CYP1A1","CYP1B1","CPSS23","CTTNBP2")
pa = data.frame(pap =paper)

papa = merge(pa,kl,by.x = "pap",by.y = "V4")

 rownames(papa) = papa$pap
 papa=papa[,2:21]

?pheatmap


 mf = papa[,c(2,6,14,18,4,8,16,20,3,7,15,19,1,5,13,17)]
 
 
library(pheatmap)
 
pheatmap(mf,scale = "row",cluster_cols = F,
         color = colorRampPalette(colors = c("#0DB33F","black","#EB1D25"))(50)
         )


### 2 ####

library(tidyverse)

setwd('E:/GuangZhou/Job/WG/Data')

da = data.frame()
for (i in list.files("./RNA_seq",pattern = "*.gtf"))
{
  fiedir = paste("./RNA_seq",i,sep = "/")
  d = read.table(fiedir,header = F,sep = "\t")
  d$fen = sapply(as.character(d$V2), function(x) if(grepl(".", x)){strsplit(x,"[.]")[[1]][1]})
  counid_d = d[!duplicated(d[,10]),]
  counid_d$sample = i
  da = rbind(da,counid_d[,7:11])
  rm(counid_d,d)
}
table(da$sample)

tpm = da[,c(4,3,5)]
kl = spread(
  data=tpm,
  key=sample,
  value=V9
)
kl = kl[grep("Gene",kl$fen,invert = T),]


paper = c("TIPARP","GPR68","AHRR","CYP1A1","CYP1B1","CPSS23","CTTNBP2")
pa = data.frame(pap =paper)

papa = merge(pa,kl,by.x = "pap",by.y = "fen")
papa[,2:5] = apply(papa[,2:5],2,as.numeric)

colnames(papa) = c("gene","control","heme","VM171","SR1")
row.names(papa) = papa$gene


papa = papa[,2:5]


library(pheatmap)


pheatmap(papa,scale = "row",cluster_cols = F,
         color = colorRampPalette(colors = c("#0DB33F","black","#EB1D25"))(50)
)

### for up and down

colnames(kl) = c("gene","control","heme","VM171","SR1")
rownames(kl) = kl$gene
kl[,2:5] = apply(kl[,2:5],2,as.numeric)
kl = kl[,2:5] + 0.001
kl$ratio = log2(kl$heme / kl$control)

DEgene =  kl[which(abs(kl$ratio) > 1),]


library("clusterProfiler")
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)



up = DEgene[which(DEgene$ratio > 1),]
down = DEgene[which(DEgene$ratio < -1),]
#利用bitr函数进行id转换，使用bioconductor的db系列包进行
hg<-bitr(row.names(up[order(up$ratio,decreasing = T),]),fromType="SYMBOL",toType=c("ENTREZID","SYMBOL"),OrgDb="org.Hs.eg.db")

go <- enrichGO(hg$SYMBOL,OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
write.csv(go,file="go_up.txt")
dotplot(go,showCategory=8,split = 'ONTOLOGY', title = "Up gene of heme vs control (BP,CC,MF)",)

## down
hg_down<-bitr(row.names(down[order(down$ratio,decreasing = F),]),fromType="SYMBOL",toType=c("ENTREZID","SYMBOL"),OrgDb="org.Hs.eg.db")

go <- enrichGO(hg_down$SYMBOL,OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
write.csv(go,file="go_down.txt")
dotplot(go,showCategory=8,split = 'ONTOLOGY', title = "Down gene of heme vs control (BP,CC,MF)",)


########### for edger DEG ######

library(edgeR)

setwd('E:/GuangZhou/Job/WG/Data')
sw <- read.table("./RNA_seq/featureCounts_m.txt",header = T,sep = "\t")
sw = sw[,c(1,7,8)]
sw$fen = sapply(as.character(sw$Geneid), function(x) if(grepl(".", x)){strsplit(x,"[.]")[[1]][1]})
counid_d = sw[!duplicated(sw[,4]),]
row.names(counid_d) = counid_d[,4]

group <- 1:2
y <- DGEList(counts = counid_d[,2:3],genes = counid_d[,4],group = group)
y

keep <- rowSums(cpm(y)>0) >= 1
y <- y[keep, , keep.lib.sizes=FALSE]
y

y <- calcNormFactors(y)

y_bcv <- y
bcv <- 0.2

et <- exactTest(y_bcv, dispersion = bcv ^ 2)

topTags(et,n=100)


gene1 <- decideTestsDGE(et, p.value = 0.05, lfc = 0)
summary(gene1)

library(ggplot2)
res <- et$table
ggplot(res, aes(x=logFC, y=-log10(df$PValue)) ) +
  geom_point() +
  ylab("-log10(p value)") +
  theme_bw()


detags <- rownames(y)[as.logical(gene1)]
plotSmear(et, de.tags=detags)
abline(h=c(-4, 4), col="blue")


res$Group <- "Not"
######
res$Group[which((res$PValue < 0.05) & (res$logFC > 1))] = "Up"
res$Group[which((res$PValue < 0.05) & (res$logFC < -1))] = "Down"

table(res$Group)

write.table(res, 'control_treat.txt', sep = '\t', col.names = T, quote = FALSE)


p <- ggplot(data = res, aes(x = logFC, y = -log10(PValue), color = Group)) +
  geom_point(size = 1) +  #绘制散点图
  scale_color_manual(values = c('green', 'gray', 'red'), limits = c('Down', 'Not','Up')) +  #自定义点的颜色
  labs(x = 'log2 Fold Change', y = '-log10 p-value', title = '', color = '') +  #坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线>、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
  geom_hline(yintercept = 1.3, lty = 3, color = 'black')
  #xlim(-12, 12) + ylim(0, 35)
  #pdf(paste(outname, ".volcano.pdf",sep=""))
p


### for go ##
library("clusterProfiler")
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

hg<-bitr(row.names(res[which(res$Group == "Up"),]),fromType="SYMBOL",toType=c("ENTREZID","SYMBOL"),OrgDb="org.Hs.eg.db")

go <- enrichGO(hg$SYMBOL,OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
write.csv(go,file="go_up_edger.txt")
dotplot(go,showCategory=8,split = 'ONTOLOGY', title = "Up gene of heme vs control (BP,CC,MF)",)

## down
hg_down<-bitr(row.names(res[which(res$Group == "Down"),]),fromType="SYMBOL",toType=c("ENTREZID","SYMBOL"),OrgDb="org.Hs.eg.db")

go <- enrichGO(hg_down$SYMBOL,OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.1, qvalueCutoff = 0.1,keyType = 'SYMBOL')
write.csv(go,file="go_down_edger.txt")
dotplot(go,title = "Down gene of heme vs control (BP,CC,MF)")



