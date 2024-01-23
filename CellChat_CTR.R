#devtools::install_github("sqjin/CellChat")

## create multple objects in a loop
#for (i in 0:14){
#  assign(paste0("immune", i), FindMarkers(AllCells.combined, ident.1 = "VEH", ident.2 = "IMQ", verbose = TRUE, group.by="stim", subset.ident = i) )
#}

#for (i in 0:15){
#  marker_i <- FindConservedMarkers(AllCells.combined, ident.1 = i, grouping.var = "stim", verbose =TRUE)
#  filename <- paste0("immunes.", i,".csv")
#  write.csv(marker_i, filename)
#}

#Do it as a list
#fit_ <- list()
#for (i in 1:20){
#  fit_[[i]] <- lm(y ~ x1 + x2 + x3, data=mydata)
#} 

library(CellChat)
library(patchwork)
library(ggpubr)
library(NMF)
library(ggalluvial)
library(umap)
library(reticulate)
library("xlsx")

options(stringsAsFactors = FALSE)

setwd("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All")

data.dir <- './CellChatFiles'
dir.create(data.dir)
setwd(data.dir)

CTR <- subset(all, cells = colnames(all)[all@meta.data[, "groups"] == "CTR"])
TNT <- subset(all, cells = colnames(all)[all@meta.data[, "groups"] == "TNT"])
MHC <- subset(all, cells = colnames(all)[all@meta.data[, "groups"] == "MHC"])

vector <- list()
vector <- list("CTR","TNT","MHC")

#Produce each CellChatfile
for (n in vector) {

  object <- assign(paste0("subset_",n ), subset(all, cells = colnames(all)[all@meta.data[, "groups"] == n]))
  #cellchat <- paste0("cellchat_",n)
  #object <- paste0("subset_",n)

Idents(object = object) <- "ArchR_Clusters"
my_levels <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25")
Idents(object) <- factor(Idents(object), levels= my_levels)

data.input <- GetAssayData(object, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(object)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels


assign(paste0("cellchat_",n ), createCellChat(object = data.input, meta = meta, group.by = "group"))

#cellchat <- addMeta(cellchat, meta = meta, meta.name = "cluster")
#cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
#levels(cellchat@idents) # show factor levels of the cell labels

rm(data.input)
rm(object)
rm(meta)
rm(list=paste0("subset_",n ))

}

#Choose MouseDatabase
CellChatDB <- CellChatDB.mouse #CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

genotype <- list()
genotype <- list("CTR","MHC","TNT")

cellchat_MHC@DB <- CellChatDB.use
cellchat_CTR@DB <- CellChatDB.use

#list.mods=lapply(ls(pattern="cellchat_"),get)
#names(list.mods) <- genotype

#for (i in genotype )  {
# set the used database in the object

data.dir <- './CellchatCTR'
dir.create(data.dir)
setwd(data.dir)

cellchat_CTR@DB <- CellChatDB.use

###########################################################################


# subset the expression data of signaling genes for saving computation cost
cellchat_CTR <- subsetData(cellchat_CTR) # This step is necessary even if using the whole database)
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat_CTR <- identifyOverExpressedGenes(cellchat_CTR)
cellchat_CTR <- identifyOverExpressedInteractions(cellchat_CTR)
# project gene expression data onto PPI network (optional)
cellchat_CTR <- projectData(cellchat_CTR, PPI.mouse)
#print("done for_",i)
#}
#######Inference of cell-cell communication network##
cellchat_CTR <- computeCommunProb(cellchat_CTR, population.size = TRUE,)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_CTR <- filterCommunication(cellchat_CTR, min.cells = 10)



##############################################################################
###########Extract the inferred cellular communication network as a data frame

#data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
rec_lig_CTR <- subsetCommunication(cellchat_CTR)
#Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways. 
sig_path_CTR <- subsetCommunication(cellchat_CTR,slot.name = "netP")
#gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
spec_pathway_CTR <- subsetCommunication(cellchat_CTR, signaling = c("FGF", "TGFb")) 

write.xlsx(rec_lig_CTR, file= "CTR_ligand_receptor_cellcell_communication.xlsx")
write.xlsx(sig_path_CTR, file= "CTR_singnaling_pathway_communication.xlsx")
write.xlsx(spec_pathway_CTR, file= "CTR_inferred_cellcell_communication_FGF_TGFb.xlsx")

rm(rec_lig_CTR)
rm(sig_path_CTR)
rm(spec_pathway_CTR)

#Infer the cell-cell communication at a signaling pathway level
cellchat_CTR <- computeCommunProbPathway(cellchat_CTR)


#Calculate the aggregated cell-cell communication network
#USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
cellchat_CTR <- aggregateNet(cellchat_CTR)

#Visualize  the aggregated cell-cell communication network. For example, showing the number of interactions 
#or the total interaction strength (weights) between any two cell groups using circle plot.


groupSize <- as.numeric(table(cellchat_CTR@idents))
par(mfrow = c(1,2), xpd=TRUE)
png("CTR_Network.png", width = 960, height = 960, units = "px")
netVisual_circle(cellchat_CTR@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
png("CTR_Network_Interaction_weigthed_strength.png", width = 960, height = 960, units = "px")
netVisual_circle(cellchat_CTR@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


#############################################################################
# 

data.dir <- './InfluenceofCluster'
dir.create(data.dir)
setwd(data.dir)

mat <- cellchat_CTR@net$weight

par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  
  Cluster_Influence <- paste0("CTR_Influence of Cluster_",i,".png")
    
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  png(Cluster_Influence, width = 480, height = 480, units = "px")
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    dev.off()
}

setwd("../")

#######################
#######Individual Plotting of Pathways and Auto and Paracrine Function

data.dir <- './AutoandParacrineFunction'
dir.create(data.dir)
setwd(data.dir)

pathways <- c("FGF","TGFb","PERIOSTIN")
vertex.receiver = seq(6,9)

for (k in pathways ) {
pathways.show <- k

print(pathways.show)


# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
 # a numeric vector. 
png(paste0("CTR_",k,"_signaling_from_Cardiomyocyte_clusters.png"), width = 1440, height = 960, units = "px")
netVisual_aggregate(cellchat_CTR, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

# Circle plot
png(paste0("CTR_",k,"_signaling_from_Cardiomyocyte_clusters_Circle_Plot.png"), width = 1440, height = 960, units = "px")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_CTR, signaling = pathways.show, layout = "circle")
dev.off()

# Chord diagram
png(paste0("CTR_",k,"_signaling_from_Cardiomyocyte_clusters_Chord_diagramm.png"), width = 1440, height = 960, units = "px")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_CTR, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
dev.off()


# Heatmap
png(paste0("CTR_",k,"_TGFb_signaling_from_Cardiomyocyte_clusters_Heatmap.png"), width = 1920, height = 1920, units = "px")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_CTR, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()



}

setwd("../")
# Chord diagram seperated
#group.cellType <- c(rep("Leukocytes", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
#names(group.cellType) <- levels(cellchat@idents)
#netVisual_chord_cell(cellchat_CTR, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.



##########################
######Compute the contribution of each ligand-receptor pair to the overall signaling pathway 
# and visualize cell-cell communication mediated by a single ligand-receptor pair

data.dir <- './Contributionof_LR_to_SignalP'
dir.create(data.dir)
setwd(data.dir)

pathways <- c("FGF","TGFb","PERIOSTIN")
vertex.receiver = seq(6,9)

for (k in pathways) {
  pathways.show <- k
  
  print(pathways.show)


#Bargraph
png(paste0("CTR_Contribution_of_L-R_pair_to_pathway_",k,"_from_C",vertex.receiver,".png"), width = 600, height = 600, units = "px")
print(netAnalysis_contribution(cellchat_CTR, signaling = pathways.show))
dev.off()


pairLR.CXCL <- extractEnrichedLR(cellchat_CTR, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
png(paste0("CTR_Contribution_of_top_L-R_pair_to_pathway_",k,"_from_C",vertex.receiver,"_HierarchyPlot.png"), width = 600, height = 600, units = "px")
print(netVisual_individual(cellchat_CTR, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver, layout = "hierarchy"))
dev.off()
# Circle plot
png(paste0("CTR_Contribution_of_top_L-R_pair_to_pathway_",k,"_from_C",vertex.receiver,"_CirclePlot.png"), width = 600, height = 600, units = "px")
print(netVisual_individual(cellchat_CTR, signaling = pathways.show, pairLR.use = LR.show, layout = "circle"))
dev.off()
# Chord diagram
png(paste0("CTR_Contribution_of_top_L-R_pair_to_pathway_",k,"_from_C",vertex.receiver,"_ChordPlot.png"), width = 600, height = 600, units = "px")
print(netVisual_individual(cellchat_CTR, signaling = pathways.show, pairLR.use = LR.show, layout = "chord"))
dev.off()
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

}
setwd("../")
######################################################################

data.dir <- './Pathway-Loop Fibroblasts'
dir.create(data.dir)
setwd(data.dir)
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat_CTR@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(21,23)
setwd("./Pathway-Loop Fibroblasts/")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat_CTR, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat_CTR, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
setwd("../")

##################################################
#manual interaction plot
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
gg <- netVisual_bubble(cellchat_CTR, sources.use = 8, targets.use = c(21:23), remove.isolate = FALSE)
#> Comparing communications on a single object
ggsave(filename=paste0("CTR_Interaction_8to212223_bubbleplot.pdf"), plot=gg, width = 3, height = 10, units = 'in', dpi = 300)



# show all the significant interactions (L-R pairs) associated with certain signaling pathways
gg <- netVisual_bubble(cellchat_CTR, sources.use = 9, targets.use = c(13:17), signaling = c("FGF","CXCL","TGFb", "PTPRM","LAMININ","VEGF"), remove.isolate = FALSE)
#> Comparing communications on a single object
ggsave(filename=paste0("CTR_Singnalingways_C9-13-17_FGF_CXCL_TGFb_Ptprm_LAMA_VEGF_bubbleplot.pdf"), plot=gg, width = 3, height = 8, units = 'in', dpi = 300)


# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB - CHORD-Plot

png("CTR_sig_Int_L-R pair_C8 - 22_ChordPlot.png", width = 1000, height = 1000, units = "px")
netVisual_chord_gene(cellchat_TNT, sources.use = 8, targets.use = c(22), lab.cex = 1,legend.pos.y = 30)
dev.off()


# show all the interactions received by Inflam.DC
png("CTR_Interactionsreceviced_C8 - C22_ChordPlot.png", width = 1000, height = 1000, units = "px")
netVisual_chord_gene(cellchat_CTR, sources.use = c(8), targets.use = 22, legend.pos.x = 15, lab.cex = 1.1)
dev.off()

png("CTR_Signalingpathways_TGFb_FGF_ChordPlot.png", width = 1000, height = 1000, units = "px")
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat_CTR, sources.use = c(8,9), targets.use = c(21:23), signaling = c("TGFb","FGF"),legend.pos.x = 8)
#> Note: The second link end is drawn out of sector 'CXCR4 '.
#> Note: The first link end is drawn out of sector 'CXCL12 '.
dev.off()


png("CTR_Signalingpathways_TGFb_FGF_VEGF_C89_c21_C23_ChordPlot.png", width = 1000, height = 1000, units = "px")
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat_CTR, sources.use = c(8,9), targets.use = c(21:23), signaling = c("VEGF","TGFb","FGF"),legend.pos.x = 8)
#> Note: The second link end is drawn out of sector 'CXCR4 '.
#> Note: The first link end is drawn out of sector 'CXCL12 '.
dev.off()

#######
#Plot the signaling gene expression distribution using violin/dot plot
png("CTR_ViolonPlot_TGFb.png", width = 1000, height = 1000, units = "px")
plotGeneExpression(cellchat_CTR, signaling = "TGFb")
dev.off()
#> Registered S3 method overwritten by 'spatstat':
#>   method     from
#>   print.boxx cli

pathways.show <- list()
pathways.show <- c("VEGF")
# Compute the network centrality scores
cellchat_CTR <- netAnalysis_computeCentrality(cellchat_CTR, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
png(paste0("CTR_Heatmap_",pathways.show,"_majorsignalingroles.png"), width = 1200, height = 1200, units = "px")
netAnalysis_signalingRole_network(cellchat_CTR, signaling = pathways.show, width = 32, height = 10, font.size = 18, font.size.title = 25)
dev.off()


#Visualize the dominant senders (sources) and receivers (targets) in a 2D space

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_CTR)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat_CTR, signaling = c("TGFb", "FGF","VEGF"))
#> Signaling role analysis on the cell-cell communication network from user's input
#gg1 + gg2
arrange <- ggarrange(gg1, gg2, ncol = 2, nrow = 1)
ggsave(filename=paste0("CTR_Senderreceiverin2D-Space.png"), plot= arrange, width = 10, height = 8, units = 'in', dpi = 300)


#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_CTR, pattern = "outgoing", width = 50, height = 40, font.size = 18, font.size.title = 25)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_CTR, pattern = "incoming", width = 50, height = 40, font.size = 18, font.size.title = 25)
png(paste0("CTR_Outgoing_Signaling.png"), width = 1700, height = 1300, units = "px")
ht1
dev.off()
png(paste0("CTR_INComing_Signaling.png"), width = 1700, height = 1300, units = "px")
ht2
dev.off()


# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat_CTR, signaling = c("TGFB", "VEGF", "FGF"))
png(paste0("CTR_SignalingroleAnalysis.png"), width = 1200, height = 1200, units = "px")
print(ht)
dev.off()

#################################################################################
# Identify and visualize incoming and outgoing communication pattern of secreting cells

direct <- c("outgoing","incoming")
data.dir <- './incomingoutgoing_pattern/'
dir.create(data.dir)
setwd(data.dir)

for (i in direct) { 

#selectK(cellchat_CTR, pattern = i)

nPatterns = 5

png(paste0("CTR_",i,"_communication_Pattern-n5.png"), width = 1200, height = 1000, units = "px")
cellchat_CTR <- identifyCommunicationPatterns(cellchat_CTR, pattern = i, k = nPatterns, width = 12, height = 21, font.size = 9)
dev.off()
print(paste0("1",i))

# river plot_patter_Riverplot
png(paste0("CTR_",i,"_communication__RiverPlot-n5.png"), width = 1700, height = 1300, units = "px")
print(netAnalysis_river(cellchat_CTR, pattern = i,font.size = 6.5, font.size.title = 18))
dev.off()

#> Please make sure you have load `library(ggalluvial)` when running this function
print(paste0("2",i))

# dot plot
png(paste0("CTR_",i,"_communication_DotPlot-n5.png"), width = 1000, height = 800, units = "px")
print(netAnalysis_dot(cellchat_CTR, pattern = i, font.size = 10))
dev.off()
print(paste0("3",i))
}

setwd("../")

################################################
##Manifold and classification learning analysis of signaling networks#
#Further, CellChat is able to quantify the similarity between all significant signaling pathways and then 
# group them based on their cellular communication network similarity. Grouping can be done either based on 
#the functional or structural similarity.

direct <- c("functional","structural")
data.dir <- './functional_and_structural_pattern/'
dir.create(data.dir)
setwd(data.dir)

for (i in direct) { 


#Identify signaling groups based on their functional and structual similarity
cellchat_CTR <- computeNetSimilarity(cellchat_CTR, type = i)
cellchat_CTR <- netEmbedding(cellchat_CTR, type = i)
#> Manifold learning of the signaling networks for a single dataset
cellchat_CTR <- netClustering(cellchat_CTR, type = i)
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
png(paste0("CTR_",i,"_similarity_group.png"), width = 1000, height = 800, units = "px")
print(netVisual_embedding(cellchat_CTR, type = i, label.size = 3.5))
dev.off()

png(paste0("CTR_",i,"__similarity_group_zoomin.png"), width = 1000, height = 800, units = "px")
print(netVisual_embeddingZoomIn(cellchat_CTR, type = i, nCol = 2))
dev.off()


}
setwd("../")


saveRDS(cellchat_CTR, file = "cellchat_CTR.rds")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
###################################################################################
###################################################################################
###################################################################################
