rm(list = ls()) #清除缓存

#修改1：设定工作目录
setwd("C:/Users/Administrator/Desktop/SLE/procedure/step2. analysis of all prescription")

library(tidyverse)
library(dplyr)
library(arules) #分析关联规则
library(arulesViz) #根据关联规则作图
library(RColorBrewer)
library(shinythemes)# 交互分析
library(showtext)
showtext_auto(enable=T)

#修改2：设定方子数量
N=239

#未转换为事务型数据#
transdata0 <- read.transactions("C:/Users/Administrator/Desktop/SLE/procedure/step1. prescriptions preprocessing/drug.csv", format =c("basket"), 
                                sep = ",", cols = 1, header = TRUE, rm.duplicates = TRUE)

#看看数据#
inspect(transdata0[1:10])


#统计
herbfreq <- as.data.frame(sort(itemFrequency(transdata0, "absolute"), decreasing = TRUE)) 
colnames(herbfreq)<-"freq"
herbfreq$'percent(%)' <- round(herbfreq$freq/N*100,digits = 2) 
herbfreq <- herbfreq %>% rownames_to_column("drug")  %>% mutate(drug = trimws(gsub("[?]", " ", drug))) %>% 
  mutate(drug = trimws(gsub("<a1><af>", "'", drug)))

#修改3：挑选高频中药
top=nrow(herbfreq[which(herbfreq$freq>=48),])  #>=n*20%
#可视化
#高频中药频数分布图
pdf("./results/Figure 3A.pdf", height = 8, width = 14)
itemFrequencyPlot(transdata0, topN = top,
                        type ="absolute",
                        horiz = F, col = "#B0C4DE",
                        xlab="Herb",ylim=c(0,round(herbfreq[1,2],digits = -1)),
                        cex.names=0.8, cex.axis=0.9,
                        ylab="Frequency",
                       # main="A. Distribution chart of high-frequency TCM item counts"
                        
)
dev.off()


herbfreq[122,]$drug<-"Fragrant Solomon’S Seal Rhizome"
herbfreq[166,]$drug<-"Duchesneae  Indicae  Herba"

#四气五味归经
info <-  read.table("info.txt",header = TRUE, sep="\t") 
info$drug <- trimws(info$drug)
herbfreq_out <- merge(herbfreq, info, by="drug") %>% arrange(desc(freq))
write.table(herbfreq_out,file="./results/herbfreq.txt",sep="\t",quote=F,col.names=T,row.names = F)
#not <- herbfreq[-which(herbfreq$drug %in% info$drug),]

#Five Flavors#
ww1 <- aggregate(herbfreq_out$freq, by=list(type=herbfreq_out$Five.flavours.1),sum) %>% filter(type != "")
ww1 <- t(data.frame(ww1,row.names=1,stringsAsFactors = F)) %>% as.data.frame() 

ww2 <- aggregate(herbfreq_out$freq, by=list(type=herbfreq_out$Five.flavours.2),sum) %>% filter(type != "")
ww2 <- t(data.frame(ww2,row.names=1,stringsAsFactors = F)) %>% as.data.frame() 

ww3 <- aggregate(herbfreq_out$freq, by=list(type=herbfreq_out$Five.flavours.3),sum) %>% filter(type != "")
ww3 <- t(data.frame(ww3,row.names=1,stringsAsFactors = F)) %>% as.data.frame()

ww <- bind_rows(ww1,ww2,ww3)
ww[is.na(ww)]<- 0
ww["total",] <- colSums(ww) 
colnames(ww) <- trimws(colnames(ww))
ww <- ww %>% rownames_to_column("ww") %>% filter(ww=="total") %>%
  select(ww,Sour,Bitter,Sweet,Pungent,Salty,Bland,Astringent)
  
write.table(ww,file="./results/ww.txt",sep="\t",quote=F,col.names=T,row.names = F)

#绘制五味雷达图
library(ggradar)
library(ggplot2)
P1 <- ggradar(
  ww, 
  values.radar = c("0", "1000", "2000"),
  grid.min = 0, grid.mid = 1000, grid.max = 2000,
  base.size = 10,axis.label.size = 4,grid.label.size = 5,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = "#00AFBB",
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "#00AFBB",
  gridline.min.colour = "grey"
) + ggtitle("Frequency of herbs' flavors") + theme(
  plot.title = element_text(size = 13, face = "bold",family = "times new roman", hjust = 0.5, vjust = 0.5) ,

  )

#Natures#
yx <- aggregate(herbfreq_out$freq, by=list(type=herbfreq_out$Four.qi),sum) %>% filter(type != "")
yx1 <- t(data.frame(yx,row.names=1,stringsAsFactors = F)) %>% as.data.frame() %>% rownames_to_column("yx") %>%
  select(yx,Cold,Heat,Warm,Cool,Neutral)
write.table(yx,file="./results/yx.txt",sep="\t",quote=F,col.names=T,row.names = F)
#绘制药性雷达图
P2 <- ggradar(
  yx1, 
  values.radar = c("0", "1000", "2000"),
  grid.min = 0, grid.mid = 1000, grid.max = 2000,
  base.size = 10,axis.label.size = 4,grid.label.size = 5,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = "#E7B800",
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "#E7B800",
  gridline.min.colour = "grey"
)+ ggtitle("Frequency of herbs' natures") + theme(
  plot.title = element_text(size = 13, face = "bold",family = "times new roman", hjust = 0.5, vjust = 0.5) 
)

#Meridian tropism#
gj1 <- aggregate(herbfreq_out$freq, by=list(type=herbfreq_out$Meridian.affinity.1),sum) %>% filter(type != "")
gj1 <- t(data.frame(gj1,row.names=1,stringsAsFactors = F)) %>% as.data.frame() 

gj2 <- aggregate(herbfreq_out$freq, by=list(type=herbfreq_out$Meridian.affinity.2),sum) %>% filter(type != "")
gj2 <- t(data.frame(gj2,row.names=1,stringsAsFactors = F)) %>% as.data.frame() 

gj3 <- aggregate(herbfreq_out$freq, by=list(type=herbfreq_out$Meridian.affinity.3),sum) %>% filter(type != "")
gj3 <- t(data.frame(gj3,row.names=1,stringsAsFactors = F)) %>% as.data.frame()

gj4 <- aggregate(herbfreq_out$freq, by=list(type=herbfreq_out$Meridian.affinity.4),sum) %>% filter(type != "")
gj4 <- t(data.frame(gj4,row.names=1,stringsAsFactors = F)) %>% as.data.frame()

gj5 <- aggregate(herbfreq_out$freq, by=list(type=herbfreq_out$Meridian.affinity.5),sum) %>% filter(type != "")
gj5 <- t(data.frame(gj5,row.names=1,stringsAsFactors = F)) %>% as.data.frame()

gj6 <- aggregate(herbfreq_out$freq, by=list(type=herbfreq_out$Meridian.affinity.6),sum) %>% filter(type != "") 
colnames(gj6)[2] <- gj6$type
gj6 <- gj6 %>% select(`Large intestine`)

gj <- bind_rows(gj1,gj2,gj3,gj4,gj5,gj6)
gj[is.na(gj)]<- 0
gj["total",] <- colSums(gj) 
gj <- gj %>% rownames_to_column("gj") %>% filter(gj=="total") %>%
  select(gj,Heart,Liver,Spleen,Lung,Kidney,Stomach,`Large intestine`,`Small intestine`,Bladder,Gallbladder,Pericardium,Sanjiao) 
write.table(gj,file="./results/gj.txt",sep="\t",quote=F,col.names=T,row.names = F)
#绘制归经雷达图
P3 <- ggradar(
  gj, 
  values.radar = c("0", "900", "1800"),
  grid.min = 0, grid.mid = 900, grid.max = 1800,
  base.size = 10,axis.label.size = 4,grid.label.size = 5,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = "#FC4E07",
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "#FC4E07",
  gridline.min.colour = "grey"
)+ ggtitle("Frequency of herbs' meridian affinities") + theme(
  plot.title = element_text(size = 13, face = "bold",family = "times new roman", hjust = 0.5, vjust = 0.5) 
)



library(gridExtra)
pdf("./results/Figure 3B.pdf", height = 8, width = 19)
grid.arrange(P1,P2,P3,nrow=1)
dev.off()


####分析关联规则####
##构建模型，设置参数
##support为支持度，confidence为置信度，lift为提升度
myrules <- apriori(transdata0, parameter = list( support = 0.07, confidence = 0.6, target = "rules"))                   
##查看模型摘要
summary(myrules)

## 对得出的关联规则按照置信度进行排序，也可根据需要根据支持度或者提升度进行排序
myrules.a <- sort(myrules, by = "confidence", decreasing = T ) 
## 排序后查看结果
a <- inspect(myrules.a[1:50])
rules_top50 <- as(myrules.a[1:50],'data.frame') %>% as_tibble()
write.table(rules_top50,file="./results/rules_top50.txt",sep="\t",quote=F,col.names=T,row.names = T)

##关联规则可视化
#图1C
pdf("./results/Figure 3C.pdf", height = 8, width = 14) 
plot(myrules.a[1:50],method = "graph",measure = "confidence",
     control = list(shading="lift"),colors=c("#76EEC6","#FAEBD7"))
dev.off()



