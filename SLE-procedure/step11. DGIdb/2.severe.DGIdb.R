rm(list = ls()) #清除缓存

setwd("C:/Users/Administrator/Desktop/SLE/procedure/step11. DGIdb")
library(statnet)
library(circlize)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

zx_gene <- read.table("severe.txt", header=F, sep="\t", check.names=F)
zx_interaction0 <- read.table("interactions.csv", header=T, sep=",", check.names=F) %>%
  filter(Name %in% zx_gene$V1) %>%
  arrange(Name,desc(`Interaction Score`)) %>%
  group_by(Name) %>%
  filter(row_number()<=4) %>%
  select(-`Interaction Type & Directionality`) %>%
  spread(key = "Name",
         value = `Interaction Score`)  %>% 
  arrange(Drug) %>%
  column_to_rownames("Drug") 

# 按照所有列排序
zx_interaction<- zx_interaction0[do.call(order, lapply(zx_interaction0,function(x)-x)),] %>%
  as.matrix()

rownames(zx_interaction) <- rownames(zx_interaction0)
colnames(zx_interaction) <- colnames(zx_interaction0)

nm <- unique(unlist(dimnames(zx_interaction)))
group <- c(rep("Drug",nrow(zx_interaction)),rep("Gene",ncol(zx_interaction)))
names(group) <- nm

grid.col = NULL
grid.col[rownames(zx_interaction)] = brewer.pal(nrow(zx_interaction),"Set3") 
grid.col[colnames(zx_interaction)] = brewer.pal(nrow(zx_interaction),"Set2") 


pdf(file="Figure 7J.severe_gene_drug.circ.pdf", width=10, height=10)
chordDiagram(zx_interaction, annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col,group=group)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .2, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex=0.8)
  circos.axis(h = "top", labels.cex = 0.6, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

