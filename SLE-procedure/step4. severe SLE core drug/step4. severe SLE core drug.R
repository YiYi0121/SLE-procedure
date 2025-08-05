rm(list = ls()) #清除缓存

#修改1：设定工作目录
setwd("C:/Users/Administrator/Desktop/SLE/procedure/step4. severe SLE core drug")

library(tidyverse)
library(dplyr)
library(arules) #分析关联规则
library(arulesViz) #根据关联规则作图
library(RColorBrewer)
library(shinythemes)# 交互分析
library(showtext)
showtext_auto(enable=T)
#修改2：设定方子数量
N=136

#未转换为事务型数据#
transdata0 <- read.transactions("drug.csv", format =c("basket"), 
                                sep = ",", cols = 1, header = TRUE, rm.duplicates = TRUE)

#看看数据#
inspect(transdata0[1:10])


#统计
herbfreq <- as.data.frame(sort(itemFrequency(transdata0, "absolute"), decreasing = TRUE)) 
colnames(herbfreq)<-"freq"
herbfreq$'percent(%)' <- round(herbfreq$freq/N*100,digits = 2) 
herbfreq <- herbfreq %>% rownames_to_column("drug")
#修改3：挑选高频中药
top <- nrow(herbfreq[which(herbfreq$freq>=27),])
hmod <- herbfreq[which(herbfreq$freq>=27),]$drug
#可视化
#高频中药频数分布图
pdf("./results/Figure 4F.severe SLE. high Frequency.pdf", height = 8, width = 14)
itemFrequencyPlot(transdata0, topN = top,
                        type ="absolute",
                        horiz = F, col = "#B0C4DE",
                        xlab="Herb",ylim=c(0,round(herbfreq[1,2],digits = -1)),
                        cex.names=0.8, cex.axis=0.9,
                        ylab="Frequency",
                        #main="A. Distribution chart of high-frequency TCM item counts"
                        
)
dev.off()


####分析关联规则####
##构建模型，设置参数
##support为支持度，confidence为置信度，lift为提升度
myrules <- apriori(transdata0, parameter = list( support = 0.07, confidence = 0.6, target = "rules"))                   
##查看模型摘要
summary(myrules)

#高关联规则分析
strong_rules <- apriori(transdata0, parameter = list( support = 0.1, confidence = 0.8,
                                                       maxlen = 10, minlen = 2, target = "rules"))
strong_rules <- sort(strong_rules, by = "support", decreasing = T ) 
strong_rules <-  subset(strong_rules, lift >= 3)
rules_s <- as(strong_rules,'data.frame') %>% as_tibble() %>%
     separate(rules, c("lhs","rhs"), sep = "=>") 

write.table(rules_s,file="./results/severe.strongrules.txt",sep="\t",quote=F,col.names=T,row.names = T)


##强规则项中包含的药物
smod=unique(
  trimws(
    c(
      unique(unlist(strsplit(gsub("\\{|\\}","",rules_s$lhs),","))),
      unique(unlist(strsplit(gsub("\\{|\\}","",rules_s$rhs),",")))
    )))


smod

##关联规则可视化
pdf("./results/Figure 4G.severe SLE.strongrules TOP50.pdf", height = 8, width = 14)
# 设置字体大小为1.5
par(cex = 1.5)
plot(strong_rules,method = "graph",measure = "confidence",
     control = list(shading="lift"),colors=c("#76EEC6","#FAEBD7"))
dev.off()


#聚类分析数据准备
drug <- read.table("drug.csv",sep=",",header=T,check.names=F, row.names = 1)  %>% t()  %>% as_tibble() %>%
  sapply(table) %>% unlist() 

b <- names(drug) %>% as_tibble() %>% 
  separate(value, c("prescription", "drug"),sep="\\.") %>%
  filter(drug != " " ) %>%
  mutate(count=1) %>%
  spread(drug,count) 

b[is.na(b)]<-0
b<-data.frame(b,stringsAsFactors = F)  %>% column_to_rownames("prescription")

colsum <- colSums(b)

c<- rbind(b,colsum) 
cc <- c[,which((c[nrow(c),] >7))] #一个药至少在5%以上方子里出现过
database <- cc[-nrow(cc),] %>% t() 
rownames(database)=str_replace_all(rownames(database),"\\W+"," ")

#应用elbow method，Silhouette method，Gap statistic三种方法综合确定聚类的最佳簇数
#参考：https://www.programminghunter.com/article/94701468914/
#画完图后根据最佳聚类数改参数xintercept
library(factoextra)
nb1 <- fviz_nbclust(database, hcut, method = "wss")+
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method",title = "")

nb2 <- fviz_nbclust(database, hcut, method = "silhouette")+
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Silhouette method",title = "")


nb3 <- fviz_nbclust(database, hcut, method = "gap_stat")+
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Gap statistic",title = "")

#选定K值等于4
library(cowplot)
pdf("./results/Figure 4H.severe SLE.best K-value.pdf", height = 4, width = 12)
plot_grid(nb1,nb2,nb3,
          #labels = c("A", "B","C"),
          nrow=1)
dev.off()


#层次聚类分析http://127.0.0.1:18735/graphics/plot_zoom_png?width=1745&height=917
hc <- hclust(dist(database,method = "binary"),method = "complete")
#hc[["labels"]] <- str_replace_all(hc[["labels"]],"\\W+"," ")
kmod <- as.data.frame(cutree(hc, k=4)) #根据最佳K值修改K值
colnames(kmod)="cluster"
kmod <- kmod %>% rownames_to_column("drug") %>% arrange(cluster) 
#mutate(drug=str_replace_all(drug,"\\W+"," "))

library(ape)
#最佳K值是4
mypal = c("#556270","#FF6B6B","#4ECDC4","#C44D58")
clus4 = cutree(hc,4) #修改最佳K值
#op = par(bg="#E8DDCB")
pdf("./results/Figure 4I. severe SLE.Cluster Dendrogram.pdf", height = 12, width = 12)  
plot(as.phylo(hc),cex=0.8,type="cladogram",tip.color = mypal[clus4],  col = "red")
dev.off() 

#最大簇聚类找到的核心药物
#核心药物满足：freq>=20、高关联规则（support > 0.1, confidence > 0.8, lift > 3)、最大簇聚类药物
jmod=kmod[which(kmod$cluster==1),]$drug
jmod
smod
hmod
cmod <- intersect(jmod,intersect(smod,hmod)) %>% as.data.frame() #cmod就是核心药物 
colnames(cmod) <- "drug"
cmod_zx<-merge(cmod,herbfreq,by = 'drug',all.x=TRUE,all.y=FALSE) %>% arrange(desc(freq))
cmod_zx
write.table(cmod_zx,file="./results/severe SLE.core drug.txt",sep="\t",quote=F,col.names=T,row.names = T)

#绘制韦恩图
library(ggVennDiagram)
hx <- list(
  A = hmod, 
  B = smod,
  C = jmod
)
pdf("./results/Figure 4J.severe SLE.core drug.pdf", height = 8, width = 8) 
ggVennDiagram(hx, 
                    label_alpha = 0,  #设置背景透明度
                    label = "count", #c("both", "count", "percent", "none") 四选一，默认选第一个
                    edge_size = 0.1, edge_lty = "solid",
                    category.names = c(paste0("Inclusion criteria 1:","\n","Freqency >=27"),paste0("Inclusion criteria 2:","\n","In accord with strong association rules"),paste0("Inclusion criteria 3:","\n","In maximum cluster clustering")),
              set_size = 4)+
  scale_fill_distiller(palette = "Blues", direction = -3) +
  scale_color_brewer(palette = "Paired")+
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  guides(fill=FALSE) #移除图例
dev.off()
#Qualitative
#Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3

#Sequential
#Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd

#绘制venn图---轻型核心药物和重型核心药物的交集
library(ggvenn) 
cmod_qx <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step3. mild SLE core drug/results/mild SLE.core drug.txt",header = TRUE, sep="\t") 

hxj=list()
hxj[["Mild SLE"]]=cmod_qx$drug
hxj[["Severe SLE"]]=cmod_zx$drug

mycol=c("skyblue","pink")

pdf(file="./results/Figure 4K.mild and severe core drug ggvenn.pdf", width=6, height=6)
ggvenn(hxj, show_elements = TRUE,fill_color = c("skyblue", "pink"),
       label_sep = "\n", stroke_size = 0.5,set_name_size = 3,text_color = "black",
       text_size = 2.5)
dev.off()

#保存交集基因
intersect=intersect(cmod_qx$drug,cmod_zx$drug)
intersect
write.table(file="./results/intersect between mild and severe.txt", intersect, sep="\t", quote=F, col.names=F, row.names=F)



