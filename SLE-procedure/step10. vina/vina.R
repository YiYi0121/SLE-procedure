rm(list = ls()) #清除缓存
#修改1：设定工作目录
setwd("C:/Users/Administrator/Desktop/SLE/procedure/step10. vina")

library(tidyverse)
library(dplyr)
qx_hub <- read.table("./mild/qx.biomarker.txt",header = T, sep="\t")
zx_hub <- read.table("./severe/zx.biomarker.txt",header = T, sep="\t")

qx_drug <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step3. mild SLE core drug/results/mild SLE.core drug.txt",header = T,sep="\t")
zx_drug <- read.table("C:/Users/Administrator/Desktop/SLE/procedure/step4. severe SLE core drug/results/severe SLE.core drug.txt",header = T,sep="\t")

SwissPrediction <-  data.table::fread("C:/Users/Administrator/Desktop/SLE/procedure/Step5. TCMSP&pubchem&swissADME/CoreDrug.TargetPrediction.csv",data.table=F) %>%
  mutate(id=paste0(paste0(MOL_ID,Molecule),Common_name))

TCMSP <- data.table::fread("TCMSP.csv",data.table=F,header = T,drop=1) %>%
  rename(drug=drug_English,drug_py=drug,MolId=`MOL ID`,MolName=MOlecule) %>%
  mutate(id=paste0(paste0(drug,MolId),MolName)) %>%
  select(id,MOlecule_url,InChIKey)

#network
qx_vina <- data.table::fread("./mild/net.network.txt",data.table=F)  %>%
  filter(Node2 %in% qx_hub$id) %>%
  arrange(Node2,Node1) %>%
  mutate(mild=1) %>%
  mutate(id=paste0(paste0(Node1,MolName),Node2)) %>%
  left_join(SwissPrediction,by="id") %>%
  distinct(Node1,Node2,MolName, .keep_all = TRUE) %>%
  select(drug=drug_English,drug_py=drug,MolId=Node1,MolName,Target=Node2,TargetFullName=Target,Uniprot_ID) %>%
  mutate(id=paste0(paste0(drug,MolId),MolName)) %>%
  left_join(TCMSP,by="id") %>%
  select(drug,drug_py,MolId,MolName,MOlecule_url,InChIKey,Target,TargetFullName,Uniprot_ID)  %>%
  filter(drug %in% qx_drug$drug) %>%
  distinct() %>%
  arrange(Target,drug,MolId,MolName)

table(qx_vina$drug_py) %>% as.data.frame()
table(qx_vina$drug) %>% as.data.frame()
table(qx_vina$Target) %>% as.data.frame()

zx_vina <-  data.table::fread("./severe/net.network.txt",data.table=F)  %>%
  filter(Node2 %in% zx_hub$id) %>%
  arrange(Node2,Node1) %>%
  mutate(severe=1) %>%
  mutate(id=paste0(paste0(Node1,MolName),Node2)) %>%
  left_join(SwissPrediction,by="id") %>%
  distinct(Node1,Node2,MolName, .keep_all = TRUE) %>%
  select(drug=drug_English,drug_py=drug,MolId=Node1,MolName,Target=Node2,TargetFullName=Target,Uniprot_ID) %>%
  mutate(id=paste0(paste0(drug,MolId),MolName)) %>%
  left_join(TCMSP,by="id") %>%
  select(drug,drug_py,MolId,MolName,MOlecule_url,InChIKey,Target,TargetFullName,Uniprot_ID) %>%
  filter(drug %in% zx_drug$drug)  %>%
  distinct() %>%
  arrange(Target,drug,MolId,MolName)

table(zx_vina$drug_py) %>% as.data.frame()
table(zx_vina$drug) %>% as.data.frame()
table(zx_vina$Target) %>% as.data.frame()

write.table(file="qx.vina.txt", qx_vina, sep="\t", quote=F, col.names=T, row.names=F)
write.table(file="zx.vina.txt", zx_vina, sep="\t", quote=F, col.names=T, row.names=F)

