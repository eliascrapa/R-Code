library(CellChat)
library(patchwork)
library(ggpubr)
library(NMF)
library(ggalluvial)
library(ggplot2)

options(stringsAsFactors = FALSE)

setwd("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/OldAll")


data.dir <- './CellChatFiles'
dir.create(data.dir)
setwd(data.dir)

data.dir <- './comparison_CTRvsTNT'
dir.create(data.dir)
setwd(data.dir)

#netcentrality score
cellchat_CTR <- netAnalysis_computeCentrality(cellchat_CTR, slot.name = "netP") # the slot 'netP' means the inferred inte
cellchat_TNT <- netAnalysis_computeCentrality(cellchat_TNT, slot.name = "netP") # the slot 'netP' means the inferred inte
#Merging and comparing multiple sets


object.list <- list(CTR = cellchat_CTR,TNT = cellchat_TNT)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
#> An object of class CellChat created from a merged object with multiple datasets 

#data.dir <- './CellChat_Comparison/'
#dir.create(data.dir)
#setwd(data.dir)


#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

#gg1 + gg2
arrange <- ggarrange(gg1, gg2, ncol = 2, nrow = 1)
ggsave(filename=paste0("CTRvsTNTNumber_of_Interaction_Comparison.png"), plot= arrange, width = 10, height = 8, units = 'in', dpi = 300)


##Compare the number of interactions and interaction strength among different cell populations
#Differential number of interactions or interaction strength among different cell populations
#where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset
#compared to the first one.

par(mfrow = c(1,2), xpd=TRUE)
png("CTRvsTNT_Network.png", width = 960, height = 960, units = "px")
print(netVisual_diffInteraction(cellchat, weight.scale = T))
dev.off()

png("CTRvsTNT_Network_weight.png", width = 960, height = 960, units = "px")
print(netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight"))
dev.off()

#PDFs of same plots
pdf("CTRvsTNT_Network.pdf")
print(netVisual_diffInteraction(cellchat, weight.scale = T))
dev.off()

pdf("CTRvsTNT_Network_weight1122.pdf")
print(netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight"))
dev.off()



#differential number of interactions or interaction strength in a greater details using a heatmap
#In the colorbar, red(or blue) represents increased(or decreased) signaling in the second dataset 
#compared to the first one.

pdf("CTRvsTNT_VisualHeatmap.pdf")
par(mfrow = c(1,2), xpd=TRUE)

print(netVisual_heatmap(cellchat))
dev.off()

pdf("CTRvsTNT_VisualHeatmap_weight.pdf")
print(netVisual_heatmap(cellchat, measure = "weight"))
dev.off()


######DOES NOT WORK SOMEHOW
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
arrange <- ggarrange(gg1, gg2, ncol = 2, nrow = 1)
ggsave(filename=paste0("HeatmapInteraction_Comparison.pdf"), plot= arrange, width = 10, height = 10, units = 'in', dpi = 300)
ggsave(filename=paste0("HeatmapInteraction.png"), plot= gg1)#, width = 10, height = 10, device= "svg")#,units = 'in', device= "jpg")
ggsave(filename=paste0("HeatmapInteraction_weight.png"), plot= gg2, width = 10, height = 10, units = 'in', dpi = 300, device= "png")
#########



##################
#Number of Interactions

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  png(paste0("CTRvsTNT_Interactions,",i,".png"), width = 960, height = 960, units = "px") 
  print(netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i])))
  dev.off()
  }


###################################
#Differential number of interactions or interaction strength among different cell types

group.cellType <- c(rep("LEU", 5), rep("CM", 4), rep("undet",3), rep("ET", 5), rep("VSMC", 3), rep("FB", 3), rep("undetr", 2))
group.cellType <- factor(group.cellType, levels = c("LEU", "CM", "undet","ET","VSMC","FB","undetr"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  png(paste0("Number of interactions -_celltypes_labelededge",i,".png"), width = 960, height = 960, units = "px") 
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 14, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  dev.off()
}

##########################################################################################################
#differential number of interactions or interaction strength between any two cell types using circle plot
pdf("CTRvsTNT_DiffinterationNumberandStrength.pdf")
par(mfrow = c(1,2))#, xpd=TRUE)
plot(netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T))
plot(netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T))
dev.off()



#################
#Compare the major sources and targets in 2D space

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways

pdf("CTRvsTNT_MajorSourcesandTargesIn2DSpace.pdf")
plot(patchwork::wrap_plots(plots = gg))
dev.off()

############
#specific signaling changes of certain celltypes


for (n in 1:1){
pdf(paste0("CTRvsTNT_SignalingChanges_Scatterplot_C",n,".pdf"))
chosencluster <- paste0("C_",n)
print(netAnalysis_signalingChanges_scatter(cellchat, idents.use = chosencluster))#, signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
rm(chosencluster)
dev.off()
}
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "C9", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1

patchwork::wrap_plots(plots = list(gg1,gg2))

