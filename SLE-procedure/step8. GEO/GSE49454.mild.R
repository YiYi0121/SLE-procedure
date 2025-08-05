#load("./mild.immune.RData")

rm(list = ls()) #清除缓存
geo.data.series.number = "GSE49454"
setwd("C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO")



########################
Sys.setenv(LIB_XML = "$(MINGW_PREFIX)") 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("Biobase", force = TRUE)
#BiocManager::install("pkgname")
#BiocManager::install("GEOquery")
#BiocManager::install("limma")
#install.packages("gplots")
#install.packages("stringr")
#install.packages("Rcpp", dependencies = TRUE) 
#install.packages(c("devtools","usethis"))
#install.packages("backports")
#install.packages("dplyr")


#devtools::install_github("vqv/ggbiplot", force = TRUE)
#目前不能升级#
######################################
# library package
library(Biobase)
library(GEOquery)
library(limma) #标准化
library(stringr)
library(dplyr)
library(tidyverse)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(ComplexHeatmap)
library(affyPLM)


########################################
## For easier use, define arguments first
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 8)
# Differential expression analysis with limma
# load series and platform data from GEO
gset <- getGEO(geo.data.series.number, GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset[[1]]
# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))
# group names for all samples
samples <- colnames(exprs(gset))
geo.platform.number = gset@annotation    #测序平台类型gset@annotation
probes_expr = exprs(gset);dim(probes_expr) #探针表达矩阵
# probes_expr = log2(probe_expr+1)  #如果数值够大，进行log2的处理
# boxplot(probes_expr2, las = 2)


####################列出表型构建phenoType.csv文件进行分组########################
phenoDat = pData(gset) 

## selection = c("GSM724489","GSM724490","GSM724491","GSM724492","GSM724493","GSM724494","GSM724512","GSM724513","GSM724514","GSM724515","GSM724516","GSM724517","GSM724518","GSM724519","GSM724520")
## phenoDat = phenoDat %>% filter(geo_accession %in% selection) #选定特定对象


############################## 


############################## log2 transform
# log2 transform
ex = exprs(gset)
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("impute")
library(impute)
#KNN法计算缺失值,to impute missing expression data, using nearest neighbor averaging
imputed_gene_exp = impute.knn(ex,k=10,rowmax = 0.5,
                              colmax=0.8,maxp =3000, rng.seed=362436069)
ex = imputed_gene_exp$data

ex = as.data.frame(ex)
##ex = subset(ex,select= selection) # 选择特定对象
ex$ID = rownames(ex)

ids = as.data.frame(gset@featureData@data[["ID"]])
ids$Gene.symbol = gset@featureData@data[["Gene.symbol"]]  #探针号对应
colnames(ids) = c("ID","Gene.Symbol")

library(dplyr)
ids = ids[-grep('///',ids$Gene.Symbol),]          # 一个探针对应多个基因，去除
exprSet = inner_join(ids,ex,by = 'ID')          # 探针没有对应基因，去除
exprSet = na.omit(exprSet)
exprSet = exprSet[-which(exprSet$Gene.Symbol==''),]

library(limma)
exprSet= avereps(exprSet[,-c(1,2)],              # 多个探针对应一个基因，取均值
                 ID = exprSet$Gene.Symbol)
exprSet = as.data.frame(exprSet)
exprSet = normalizeBetweenArrays(exprSet)

exprSet = as.data.frame(exprSet)

library(tibble)
group <- read.table("group.txt",sep="\t",header=T,check.names=F) %>% arrange(group)
table(group$group)

exp_SLE <- exprSet[,which(colnames(exprSet) %in% group$GSM)]
# 将排序依据的列转换为向量
order_vector <- as.vector(group$GSM)
# 获取order_vector中每个元素在colnames(data)中的位置
order_index <- match(order_vector, colnames(exp_SLE))
# 按照order_index对data框中的列进行排序
exp_SLE <- exp_SLE[, order_index]
exprSet_out <- rownames_to_column(exp_SLE, var="geneNames")
exprSet_out <- exprSet_out[order(exprSet_out[,1]),]
write.table(exprSet_out, file=paste(geo.data.series.number,".txt",sep=""), row.names=F, col.names=T, quote=F, sep=",")


###############轻型SLE和正常对照差异基因分析##########################
##### 读取分组数据 ####
table(group$group)
mild.group <- group %>% filter(group=="Healthy control" | group=="Mild") %>%
  mutate(group=case_when(group=="Healthy control" ~ "control",
                         group=="Mild" ~ "mild"))
mild.exp <- exp_SLE[,which(colnames(exp_SLE) %in% mild.group$GSM)]
mild.exp.output <- mild.exp %>% rownames_to_column("Gene") 
write.table(mild.exp.output, file="./mild.results/mild.exp.txt", sep="\t", quote=F, col.names=T,row.names =F)

##### 分组矩阵 ####
mild.rt<-as.matrix(mild.exp)
mild.class <- as.character(mild.group$group)
mild.design <- model.matrix(~0+factor(mild.class)) 
colnames(mild.design) <- c("control","mild")
table(mild.group$group)


#绘制归一化结果：经过归一化处理后，在每个样本中，数据的表达平齐了
col <- c(rep("skyblue",table(mild.group$group)[1]),rep("pink",table(mild.group$group)[2]))
pdf("./mild.results/Figure 6A.mild SLE and control data set preprocessing normalization effect.pdf", height = 6, width = 6) 
par(mar=c(6,5,4,3) + 0.2)
boxplot(mild.exp,outline=FALSE,notch=T,col=col,las=2,cex.lab=0.5)
dev.off()


#differential
mild.fit <- lmFit(mild.rt,mild.design)
mild.cont<-makeContrasts(mild-control,levels=mild.design)
mild.fit2 <- contrasts.fit(mild.fit, mild.cont)
mild.fit2 <- eBayes(mild.fit2)
mild.allDiff <- topTable(mild.fit2,adjust='fdr',number=200000)
##输出差异结果
logFoldChange=0.58
adjustP=0.05
mild.diffSig <- mild.allDiff[with(mild.allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
mild.diffSig <- mild.diffSig[order(-mild.diffSig$logFC),] %>%
  mutate(regulate=case_when(logFC>0~"up",logFC<0~"down"))
table(mild.diffSig$regulate)
mild.diffSigOut <- rbind(id=colnames(mild.diffSig),mild.diffSig)
write.table(mild.diffSigOut, file="./mild.results/mild.diff.txt", sep="\t", quote=F, col.names=F)
#输出差异基因表达量
mild.diffGeneExp<- mild.rt[rownames(mild.diffSig),]
mild.diffGeneExpOut<- rbind(id=colnames(mild.diffGeneExp),mild.diffGeneExp)
write.table(mild.diffGeneExpOut,file="./mild.results/mild.diffGeneExp.txt",sep="\t",quote=F,col.names=F)
#绘制火山图
#定义显著性
library(ggplot2)
mild.Sig=ifelse((mild.allDiff$adj.P.Val<adjustP) & (abs(mild.allDiff$logFC)>logFoldChange), ifelse(mild.allDiff$logFC>logFoldChange,"Up","Down"), "Stable")
mild.hs=cbind(mild.allDiff, mild.Sig=mild.Sig)
mild.p <- ggplot(
  # 数据、映射、颜色
  mild.hs, aes(x = logFC, y = -log10(adj.P.Val), colour=mild.Sig)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-logFoldChange,logFoldChange),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(adjustP),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(Fold Change)",
       y="-log10 (p-value)")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
#输出火山图
pdf(file="./mild.results/Figure 6B.mild SLE vs control.volcano.pdf", width=6, height=6)
print(mild.p)
dev.off()

#轻型差异基因GO圈图
rt=mild.diffSig %>% rownames_to_column("id")



pvalueFilter=0.05       #p值过滤条件
qvalueFilter=0.05       #矫正后的p值过滤条件

#基因名字转换为基因id
colnames(rt)[1]="Gene"
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(GO,file="./mild.results/mild.GO.txt",sep="\t",quote=F,row.names = F)

#library(dplyr)
#GO.class <- read.table("tan.total.pathway.csv",sep = ",",header = T)  %>% distinct()
#GO.class <- GO.class[which(GO.class$ID  %in%  GO$ID),]
#GO.class.out <- GO[-which(GO$ID %in% GO.class$ID),] %>% 
#  dplyr::select(ID,Pathway.name=Description) %>%
#  mutate(class="") %>%
#  rbind(GO.class) %>% 
#  distinct()
#write.table(GO.class.out,file="./mild.results/mild.GO.class.csv",sep=",",quote=F,row.names = F)

#整理后重新读取
#GO.class1 <- read.table("./mild.results/mild.GO.class1.csv",sep = ",",header = T)  

#library(ggplot2)
#pdf("./mild.results/Figure 6C.mild.GO.class.pdf",width=8,height=8)
#ggplot(data=GO.class, mapping=aes(x=Class,fill=Class))+
#  geom_bar(stat="count",fill="pink",colour="black")+
  #geom_text(stat='count',aes(label=..count..), vjust=1.6, color="white", size=5)+
#  theme(text = element_text(size=14), axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))+
#  labs(x="Function",y="GO terms")
#dev.off()
###########绘制GO圈图###########
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")
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

#绘制圈图的主体部分
pdf("./mild.results/Figure 6C.mild.GO.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
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
  labels = c("Biological Process","Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#绘制显著性pvalue的图例
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()

#细胞死亡基因
celldeath <- read.table("celldeath.csv",sep=",",header=T,check.names=F) 
qx.celldeath <- bind_rows(data.frame(ferr=celldeath[which(celldeath$ferrgene %in% rownames(mild.diffSig)),]$ferrgene),
                      data.frame(pyro=celldeath[which(celldeath$pyrogene %in% rownames(mild.diffSig)),]$pyrogene),
                      data.frame(anoikis=celldeath[which(celldeath$anoikisgene %in% rownames(mild.diffSig)),]$anoikisgene),
                      data.frame(copper=celldeath[which(celldeath$coppergene %in% rownames(mild.diffSig)),]$coppergene) )

qx.celldeath <- as.data.frame(lapply(qx.celldeath,function(x)c(x[!is.na(x)],x[is.na(x)])))
qx.celldeath <- qx.celldeath[rowSums(is.na(qx.celldeath))!=ncol(qx.celldeath),]
qx.celldeath[is.na(qx.celldeath)]<-""

write.table(file="./mild.results/qx.celldeath.txt", qx.celldeath, sep="\t", quote=F, col.names=T, row.names=F)


#save.image("mild.diff.RData")
#########################mild SLE immunoinfiltration######################################
rm(list = ls()) #清除缓存
inputFile="C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO/mild.results/mild.exp.txt"      #表达数据文件
setwd("C:/Users/Administrator/Desktop/SLE/procedure/step8. GEO")      #设置工作目录
source("geoCRG.CIBERSORT.R")       #引用包

#免疫细胞浸润分析
outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=T)

#对样品进行分组
group <- read.table("Group.txt",sep = "\t",header = T) %>% dplyr::select(GSM,group)
conNum <- as.numeric(table(group$group)[1])
treatNum <- as.numeric(table(group$group)[2])

outTab1 <- outTab %>% as.data.frame() %>% rownames_to_column("GSM") %>%
  left_join(group,by="GSM") %>% column_to_rownames("GSM")
table(outTab1$group) #过滤前样本数

#对免疫浸润结果过滤，并且保存免疫细胞浸润结果
outTab2 <- outTab1[outTab1[,"P-value"]<0.05,]
table(outTab2$group) #过滤后样本数

outTab3=as.matrix(outTab2[,1:(ncol(outTab2)-4)])
outTab3.out=rbind(id=colnames(outTab3),outTab3)
write.table(outTab3.out, file="./mild.results/mild.CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


data <- t(outTab3)
rt <- outTab3


#把数据转换成ggplot2输入文件
library(reshape)
library(ggpubr)
Type <- outTab2$group
data1=cbind(as.data.frame(t(data)), Type)
data1=melt(data1, id.vars=c("Type"))
colnames(data1)=c("Type", "Immune", "Expression")

#输出小提琴图
#install.packages("vioplot")
library(vioplot)      
outTab2=data.frame()
pdf(file="./mild.results/Figure 6D.mild.immune.vioplot.pdf", height=8, width=13)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，正常组用蓝色表示，肿瘤组用红色表示
for(i in 1:ncol(rt)){
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }
  conData=rt[1:conNum,i]
  treatData=rt[(conNum+1):(conNum+treatNum),i]
  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'forestgreen')
  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'orange')
  wilcoxTest=wilcox.test(conData,treatData)
  p=wilcoxTest$p.value
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab2=rbind(outTab2,cellPvalue)
  }
  mx=max(c(conData,treatData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topleft", 
       c("Healthy control", "Mild SLE"),
       lwd=3,bty="n",cex=1,
       col=c("forestgreen","orange"))
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

library(dplyr)

mean = aggregate(x=data1$Expression, by=list(data1$Type,data1$Immune),mean)
sd = aggregate(x=data1$Expression, by=list(data1$Type,data1$Immune),sd)
immuneDiff = left_join(mean,sd,by=c("Group.1","Group.2"))

immuneDiff <- immuneDiff[which(immuneDiff$Group.2 %in% outTab2$Cell),] %>%
  group_by(Group.2) %>%
  dplyr::select(Cell=Group.2, group=Group.1,mean=x.x,sd=x.y) %>%
  left_join(outTab2,by="Cell")


#输出免疫细胞和p值表格文件
write.table(immuneDiff,file="./mild.results/mild.immuneDiff.xls",sep="\t",row.names=F,quote=F)

#绘制柱状图
pdf(file="./mild.results/mild.barplot.pdf", width=14.5, height=8.5)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"normal",cex=1.8)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"patient",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()

#save.image("mild.immune.RData")
