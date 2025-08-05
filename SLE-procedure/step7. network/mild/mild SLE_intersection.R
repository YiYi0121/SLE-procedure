rm(list = ls()) #清除缓存

#修改1：设定工作目录
setwd("C:/Users/Administrator/Desktop/SLE/procedure/step7. network/mild")

library(tidyverse)
library(dplyr)
#轻型核心药物靶基因
data1 <-  data.table::fread("C:/Users/Administrator/Desktop/SLE/procedure/Step5. TCMSP&pubchem&swissADME/CoreDrug.TargetPrediction.csv",data.table=F)
  
qx <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step3. mild SLE core drug/results/mild SLE.core drug.txt",header = T, sep="\t")
data1_qx <- data1[which(data1$drug_English %in% qx$drug),] 
table(data1_qx$drug_English)
table(data1_qx$Target_Class)
allTargets.symbol <- data1_qx %>% 
  select(Drug=drug_English, MolId=MOL_ID,MolName=Molecule,Symbol=Common_name)
write.table(file="./network/allTargets.symbol.txt", allTargets.symbol, sep="\t", quote=F, col.names=T, row.names=F)

#疾病相关基因
data2 <-  read.table("C:/Users/Administrator/Desktop/SLE/procedure/step6. disease-related genes/disease-related genes.txt",header = T, sep="\t") %>% 
  select(Total) %>% distinct()
write.table(file="./network/Disease.txt", data2, sep="\t", quote=F, col.names=F, row.names=F)

jj_qx <- list()
jj_qx[["Drug targets"]]=unique(data1_qx$Common_name)
jj_qx[["Disease-related genes"]]=unique(data2$Total)

#绘制venn图
library(venn)
mycol=c("#029149","#E0367A","#5D90BA","#431A3D","#91612D","#FFD121","#D8D155","#223D6C","#D20A13","#088247","#11AA4D","#7A142C","#5D90BA","#64495D","#7CC767")
pdf(file="Figure 5A.mild_diseaseRelated-DrugTargets.venn.pdf",width=5,height=5)
venn(jj_qx,col=mycol[1:length(jj_qx)],zcolor=mycol[1:length(jj_qx)],box=F)
dev.off()

#保存疾病和药物的交集基因
qx_intersectGenes=intersect(jj_qx[["Drug targets"]],jj_qx[["Disease-related genes"]])
write.table(file="./ppi/mild_diseaseRelated-DrugTargets-intersectgenes.txt", qx_intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)

