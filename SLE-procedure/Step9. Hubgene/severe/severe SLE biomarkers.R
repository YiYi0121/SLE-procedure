rm(list = ls()) #清除缓存
#修改1：设定工作目录
setwd("C:/Users/Administrator/Desktop/SLE/procedure/Step9. Hubgene/severe")

library(tidyverse)
library(dplyr)
#轻型核心药物靶基因和疾病相关基因的交集基因
zx_data1 <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step7. network/severe/ppi/severe_diseaseRelated-DrugTargets-intersectgenes.txt",header = F, sep="\t")

#GEO差异基因
zx_data2 <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO/severe.results/severe.diff.txt",header = T, sep="\t") 

zx_jj <- list()
zx_jj[["DrugTargets & DiseaseRelatedGenes"]]=unique(zx_data1$V1)
zx_jj[["GEO DEGs"]]=unique(zx_data2$id)

#绘制venn图
library(venn)
mycol=c("#5D90BA","#D8D155")
pdf(file="./results/Figure 7F.zx.venn.pdf",width=5,height=5)
venn(zx_jj,col=mycol,zcolor=mycol,box=F)
dev.off()

#保存交集基因
zx_intersect_Genes=zx_data2[which(zx_data2$id %in% intersect(zx_jj[[1]],zx_jj[[2]])),] %>% arrange(desc(logFC))
write.table(file="./results/zx_intersect_Genes.txt", zx_intersect_Genes, sep="\t", quote=F, col.names=T, row.names=F)
write.table(file="./preRcicros/gene.txt", zx_intersect_Genes, sep="\t", quote=F, col.names=T, row.names=F)

table(zx_intersect_Genes$regulate)

severe.GO <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO/severe.results/severe.GO.txt",sep="\t",header=T)

zx_intersect_GO <- severe.GO[grepl(paste(zx_intersect_Genes$id, collapse = "|"), severe.GO$geneID), ]
write.table(file="./results/zx_intersect_GO.csv", zx_intersect_GO, sep=",", quote=F, col.names=T, row.names=F)
nrow(severe.GO)
nrow(zx_intersect_GO)
#染色体上的位置:gene已保存在preRcicros，运行geoCRG08.preRcircos.pl
#得到Rcircos.geneLabel
library("RCircos") 
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t", check.names=F)
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 4
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=0.7
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)

pdf(file="./results/Figure 7H.zx.intersect.genes.RCircos.pdf", width=7, height=7)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Label.Data=read.table("./preRcicros/Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()


#交集基因热图
zx.exp <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO/severe.results/severe.diffGeneExp.txt",header = T, sep="\t") %>%
  filter(id %in% zx_intersect_Genes$id) %>%
  column_to_rownames("id") %>%
  as.data.frame()

group <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO/group.txt",header = T, sep="\t") %>%
  filter(group=="Healthy control" | group=="Severe") %>%
  mutate(group=case_when(group=="Healthy control"~"Healthy control",
                         group=="Severe"~"Severe SLE")) %>%
  dplyr::select(GSM,group) %>%
  column_to_rownames("GSM")

library(pheatmap)
pdf(file="./results/Figure 7G.heatmap.pdf", width=8, height=6)
pheatmap(zx.exp, 
         annotation=group, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =T,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
dev.off()

#细胞死亡基因
celldeath <- read.table("celldeath.csv",sep=",",header=T,check.names=F) 
zx.celldeath <- bind_rows(data.frame(ferr=celldeath[which(celldeath$ferrgene %in% zx_intersect_Genes$id),]$ferrgene),
                          data.frame(pyro=celldeath[which(celldeath$pyrogene %in% zx_intersect_Genes$id),]$pyrogene),
                          data.frame(anoikis=celldeath[which(celldeath$anoikisgene %in% zx_intersect_Genes$id),]$anoikisgene),
                          data.frame(copper=celldeath[which(celldeath$coppergene %in% zx_intersect_Genes$id),]$coppergene) )

zx.celldeath <- as.data.frame(lapply(zx.celldeath,function(x)c(x[!is.na(x)],x[is.na(x)])))
zx.celldeath <- zx.celldeath[rowSums(is.na(zx.celldeath))!=ncol(zx.celldeath),]
zx.celldeath[is.na(zx.celldeath)]<-""

write.table(file="./results/zx.celldeath.txt", zx.celldeath, sep="\t", quote=F, col.names=T, row.names=F)

#与免疫浸润细胞相关性
#引用包
library(limma)
library(reshape2)
library(tidyverse)
library(dplyr)
library(ggplot2)

zx.exp <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO/severe.results/severe.diffGeneExp.txt",header = T, sep="\t") %>%
  column_to_rownames("id")
#去除对照组样品
severe <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO/group.txt",header = T, sep="\t") %>%
  dplyr::select(GSM,group) %>%
  filter(group=="Severe")
data <- zx.exp[which(rownames(zx.exp) %in% zx_intersect_Genes$id),
               which(colnames(zx.exp) %in% severe$GSM)]  
data=t(data)

#读取免疫细胞结果文件，并对数据进行整理
immune.diff =read.table("C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO/severe.results/severe.immuneDiff.xls",sep="\t",header=T,  check.names=F)
immune=read.table("C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO/severe.results/severe.CIBERSORT-Results.txt" , header=T, sep="\t", check.names=F, row.names=1)
immune=immune[,which(colnames(immune) %in% immune.diff$Cell)]
sameSample=intersect(row.names(data), row.names(immune))
data=data[sameSample,,drop=F]
immune=immune[sameSample,,drop=F]

#相关性分析
outTab=data.frame()
for(cell in colnames(immune)){
  if(sd(immune[,cell])==0){next}
  for(gene in colnames(data)){
    x=as.numeric(immune[,cell])
    y=as.numeric(data[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, Immune=cell, cor, text, pvalue))
  }
}

#绘制相关性热图
outTab$cor=as.numeric(outTab$cor)
pdf(file="./results/Figure 7I.zx.immun-cor.pdf", width=7, height=5)
ggplot(outTab, aes(Immune, Gene)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +    #去掉背景
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),   #x轴字体
        axis.text.y = element_text(size = 8, face = "bold")) +       #y轴字体
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","correlation")) +   #设置图例
  scale_x_discrete(position = "bottom")      #X轴名称显示位置
dev.off()

#PPI交集基因
ppi <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step7. network/severe/ppi/score2.out.gene.txt",sep="\t",header=F,  check.names=F)
zx.hub <- intersect(zx_intersect_Genes$id,ppi$V1)
zx.hub
zx.biomarker <- zx_intersect_Genes[which(zx_intersect_Genes$id %in% zx.hub),]
write.table(file="./results/zx.biomarker.txt", zx.biomarker, sep="\t", quote=F, col.names=T, row.names=F)
write.table(file="C:/Users/Administrator/Desktop/SLE/procedure/step10. vina/severe/zx.biomarker.txt", zx.biomarker, sep="\t", quote=F, col.names=T, row.names=F)

length(zx.hub)
severe.GO <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO/severe.results/severe.GO.txt",sep="\t",header=T)
zx_hub_GO <- severe.GO[grepl(paste(zx.hub, collapse = "|"), severe.GO$geneID), ]
nrow(severe.GO)
nrow(zx_hub_GO)


