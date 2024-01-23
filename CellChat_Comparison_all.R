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

data.dir <- './comparison_all'
dir.create(data.dir)
setwd(data.dir)

#netcentrality score
cellchat_CTR <- netAnalysis_computeCentrality(cellchat_CTR, slot.name = "netP") # the slot 'netP' means the inferred inte
cellchat_TNT <- netAnalysis_computeCentrality(cellchat_TNT, slot.name = "netP") # the slot 'netP' means the inferred inte
cellchat_MHC <- netAnalysis_computeCentrality(cellchat_MHC, slot.name = "netP") # the slot 'netP' means the inferred inte
#Merging and comparing multiple sets


object.list <- list(CTR = cellchat_CTR,TNT = cellchat_TNT,MHC = cellchat_MHC)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
#> An object of class CellChat created from a merged object with multiple datasets 

#data.dir <- './CellChat_Comparison/'
#dir.create(data.dir)
#setwd(data.dir)


#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")

#gg1 + gg2
arrange <- ggarrange(gg1, gg2, ncol = 2, nrow = 1)
ggsave(filename=paste0("Number_of_Interaction_Comparison.png"), plot= arrange, width = 10, height = 8, units = 'in', dpi = 300)


##Compare the number of interactions and interaction strength among different cell populations
#Differential number of interactions or interaction strength among different cell populations
#where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset
#compared to the first one.

par(mfrow = c(1,2), xpd=TRUE)
png("Network2.png", width = 960, height = 960, units = "px")
print(netVisual_diffInteraction(cellchat, weight.scale = T))
#ggsave("Network.png", plot = plot, device = "pdf")
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

pdf("VisualHeatmap.pdf")
par(mfrow = c(1,2,3), xpd=TRUE)

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


for (n in 1:25){
pdf(paste0("CTRvsTNT_SignalingChanges_Scatterplot_C",n,".pdf"))
chosencluster <- paste0("C",n)
print(netAnalysis_signalingChanges_scatter(cellchat, idents.use = chosencluster))#, signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
rm(chosencluster)
dev.off()
}
#gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "C9", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
#patchwork::wrap_plots(plots = list(gg1,gg2))


########################
#Identify signaling groups based on their functional similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)



########################
#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf("CTRvsTNT_Signalinggroups_basedon_structural similatries.pdf")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 4, top.label = T)
dev.off()
#> 2D visualization of signaling networks from datasets 1 2
pdf("CTRvsTNT_Signalinggroups_basedon_structural similatries_Pairwisezoomin.pdf")
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
#> 2D visualization of signaling networks from datasets 1 2
dev.off()
#> Compute the distance of signaling networks between datasets 1 2v

#Compute and visualize the pathway distance in the learned joint manifold
#larger means bigger difference between the datasets

pdf("CTRvsTNT_Difference_in_Pathwaydistance_higher_Score more different between both groups.pdf")
rankSimilarity(cellchat, type = "functional")
dev.off()


##############
#Identify and visualize the conserved and context-specific signaling pathways
# Compare the overall information flow of each signaling pathway
#This bar graph can be plotted in a stacked mode or not. Significant signaling pathways were 
#ranked based on differences in the overall information flow within the inferred networks between 
#NL and LS skin. The top signaling pathways colored red are enriched in NL skin, and these colored 
#green were enriched in the LS skin.


gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
pdf("CTRvsTNT_Informationflow_of_each_pathways_comparison.pdf")
gg1 + gg2
dev.off()



#############################
###########################################
#Compare outgoing (or incoming) signaling associated with each cell population

library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.7.1.1010
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 
#Heatmap Outgoing Signaling Comparison
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 7, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 7, height = 15)
pdf("CTRvsTNT_Outgoing_Signalpattern_Comparison_Heatmap.pdf")#, paper = "a3")
draw(ht1 + ht2, ht_gap = unit(0.2, "cm"))
dev.off()

#Heatmap Incoming Signaling Comparison
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 7, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 7, height = 15, color.heatmap = "GnBu")
pdf("CTRvsTNT_Incoming_Signalpattern_Comparison_Heatmap.pdf")
draw(ht1 + ht2, ht_gap = unit(0.2, "cm"))
dev.off()

#Heatmap Overall Signaling Comparison

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 7, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 7, height = 15, color.heatmap = "OrRd")
pdf("CTRvsTNT_Overall_Signalpattern_Comparison_Heatmap.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


#############################
############################################
#####################################
#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
#Identify dysfunctional signaling by comparing the communication probabities

for (sourcecluster in sourceclustervector)
  
  for (targetcluster in targetclustervector)
    
    sourcecluster <- 9
targetcluster <- c(21:23)

pdf(paste0("CTRvsTNT_ligand-receptor_pairs_Cluster_",sourcecluster,"to_Cluster",paste0(targetcluster, collapse = ","),".pdf"))
netVisual_bubble(cellchat, sources.use = sourcecluster, targets.use = targetcluster,  comparison = c(1, 2), angle.x = 45)
dev.off()
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat, sources.use = sourcecluster, targets.use = targetcluster,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in TNT", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = sourcecluster, targets.use = targetcluster,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in TNT", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf(paste0("CTRvsTNT_Increased_Decreased_Signaling_Comparison_by_communication probabities",sourcecluster,"to_Cluster",paste0(targetcluster, collapse = ","),".pdf"))
gg1 + gg2
dev.off()

##############
#Identify dysfunctional signaling by using differential expression analysis

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "TNT"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "TNT",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "CTR",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


#showing upgulated and down-regulated signaling ligand-receptor pairs using bubble plot 
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = sourcecluster, targets.use = targetcluster, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = sourcecluster, targets.use = targetcluster, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pdf(paste0("CTRvsTNT_Up-Downregulated_Signaling_Comparison_by_differentialGeneExpression",sourcecluster,"to_Cluster",paste0(targetcluster, collapse = ","),".pdf"))
gg1 + gg2
dev.off()

#Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf(paste0("CTRvsTNT_Up-Downregulated_Signaling_Comparison_by_differentialGeneExpression_Chordplot",sourcecluster,"to_Cluster",paste0(targetcluster, collapse = ","),".pdf"))
netVisual_chord_gene(object.list[[2]], sources.use = sourcecluster, targets.use = targetcluster, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], sources.use = sourcecluster, targets.use = targetcluster, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()



# compare all the interactions sending from Inflam.FIB to DC cells
pdf(paste0("CTRvsTNT_Up-Interactionpathways",sourcecluster,"to_Cluster",paste0(targetcluster, collapse = ","),".pdf"))
#
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = sourcecluster, targets.use = targetcluster, lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}
dev.off()


# compare all the interactions sending from fibroblast to inflamatory immune cells
#par(mfrow = c(1, 2), xpd=TRUE)
pdf(paste0("CTRvsTNT_Up-Interactionpathways",sourcecluster,"to_Cluster",paste0(targetcluster, collapse = ","),".pdf"))
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(21:22), targets.use = c(2:5),  title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
}
dev.off()

#################################################
#######################################
#########################################
#Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram

pathways.show <- c("FGF") 

#circleplot
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

pdf(paste0("CTR_vs_TNT_Circleplot_",pathways.show,"_Pathway.pdf"), paper = "a4r")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


#Heatmap
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
pdf(paste0("CTR_vs_TNT_Heatmapt_",pathways.show,"_Pathway.pdf"), paper = "a4r")
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

# Chord diagram
pdf(paste0("CTR_vs_TNT_Chordplot_",pathways.show,"_Pathway.pdf"), paper = "a4r")

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.


# Chord diagram with grouped celltypes
group.cellType <- c(rep("LEU", 5), rep("CM", 4), rep("undet",3), rep("ET", 5), rep("VSMC", 3), rep("FB", 3), rep("undetr", 2))
names(group.cellType) <- levels(object.list[[1]]@idents)
pdf(paste0("CTR_vs_TNT_Chordplot_grouped_Celltypes_",pathways.show,"_Pathway.pdf"), paper = "a4r")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}
dev.off()

#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

##############################
#Part V: Compare the signaling gene expression distribution between different datasets
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("CTR", "TNT")) # set factor level
plotGeneExpression(cellchat, signaling = "FGF", split.by = "datasets", colors.ggplot = T)
#> The default behaviour of split.by has changed.
#> Separate violin plots are now plotted side-by-side.
#> To restore the old behaviour of a single split violin,
#> set split.plot = TRUE.
#>       
#> This message will be shown once per session.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.







