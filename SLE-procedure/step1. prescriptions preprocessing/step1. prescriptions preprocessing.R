rm(list = ls()) 

setwd("C:/Users/Administrator/Desktop/SLE/procedure/step1. prescriptions preprocessing")

library(tidyverse)
library(dplyr)
library("plyr")
#去除重复的方子和药物
#Step1: 总共353个方子,去掉一个方子里重复的药物
drug0 <- read.table("drug0.csv",sep=",",header=T,check.names=F) 

drug<-as.data.frame(matrix(nrow = 1,ncol = nrow(drug0)))
for (i in 1:ncol(drug0)) {
  a = unique(drug0[,i])  %>%  as_tibble() %>% filter(value!="") %>% arrange(value) %>% t() %>% as.data.frame()
  drug=rbind.fill(drug,a)
}

drug <- drug[-1,]

rownames(drug)<-paste("pre",1:nrow(drug),sep="")
colnames(drug)<-paste("herb",1:ncol(drug),sep="")

#Step2: 去掉无法分类的79个方子，剩余274个方子
drug1 <- drug[-which(rownames(drug) %in% c("pre1","pre14","pre21","pre26","pre34","pre35","pre52","pre56",
                  "pre58","pre59","pre64","pre65","pre66","pre68","pre71","pre79",
                  "pre81","pre84","pre90","pre94","pre95","pre96","pre98","pre101",
                  "pre102","pre117","pre118","pre119","pre125","pre126","pre128",
                  "pre144","pre145","pre147","pre152","pre156","pre158","pre165",
                  "pre170","pre177","pre178","pre197","pre207","pre209","pre211",
                  "pre212","pre235","pre247","pre249","pre250","pre251","pre258",
                  "pre260","pre262","pre264","pre265","pre266","pre267","pre279",
                  "pre280","pre282","pre286","pre290","pre295","pre297","pre301",
                  "pre302","pre316","pre317","pre321","pre322","pre325","pre326",
                  "pre332","pre335","pre337","pre338","pre339","pre341")),]

#Step3: 去掉重复的方子33个，剩余241个方子
drug2 <- drug1 %>% distinct()
drug2[is.na(drug2)] <- ""
#Step4: 去掉中药数小于4种药的方子2个，最终剩余239个方子
drug3 <- drug2 %>% filter(herb4 != "") #删掉中药数小于4种药的方子
drug3_out <- cbind(pre=rownames(drug3),drug3)
write.table(drug3_out,file="drug.txt",sep="\t",quote=F,col.names=T,row.names = F)


#Step5:将数据库细分为两个database：轻型和重型,轻型103个方子，重型136个，后面分别寻找核心药物
drug_qx <- drug3[which(rownames(drug3) %in% c("pre3","pre7","pre8","pre10","pre11","pre15",
                                             "pre22","pre25","pre27","pre30","pre31","pre36",
                                             "pre41","pre42","pre44","pre47","pre48","pre55",
                                             "pre57","pre60","pre61","pre74","pre77","pre82",
                                             "pre85","pre86","pre88","pre92","pre99","pre109",
                                             "pre110","pre111","pre114","pre115","pre116","pre120",
                                             "pre121","pre124","pre127","pre131","pre133","pre136",
                                             "pre138","pre140","pre141","pre142","pre148","pre150",
                                             "pre153","pre154","pre155","pre157","pre161","pre162",
                                             "pre164","pre166","pre171","pre173","pre176","pre179",
                                             "pre185","pre186","pre189","pre193","pre195","pre198",
                                             "pre199","pre200","pre201","pre202","pre210","pre213",
                                            " pre215","pre217","pre220","pre221","pre223","pre225",
                                            "pre226","pre230","pre233","pre234","pre236","pre237",
                                            "pre239","pre241","pre242","pre252","pre253","pre255",
                                            "pre256","pre259","pre261","pre268","pre273","pre274",
                                            "pre275","pre276","pre292","pre294","pre296","pre299",
                                            "pre300","pre303","pre304","pre305","pre308","pre311",
                                            "pre312","pre313","pre315","pre323","pre329","pre336",
                                            "pre344","pre348","pre349","pre350","pre352","pre353")),]
drug_qx_out <- cbind(pre=rownames(drug_qx),drug_qx)
write.table(drug_qx_out,file="drug_qx.txt",sep="\t",quote=F,col.names=T,row.names = F)


drug_zx <- drug3[-which(rownames(drug3) %in% rownames(drug_qx)),]
drug_zx_out <- cbind(pre=rownames(drug_zx),drug_zx)
write.table(drug_zx_out,file="drug_zx.txt",sep="\t",quote=F,col.names=T,row.names = F)


