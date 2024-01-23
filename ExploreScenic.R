#http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/Tutorials_JupyterNotebooks/SCENIC_tutorial_2-ExploringOutput.html



# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
#library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(stringr)
library(xlsx)

setwd("/Users/eliascrapa/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/SCENIC/")
vsnDir <- "."


scenicLoomPath <- file.path(vsnDir, "output/10K 5w dataset.loom")
motifEnrichmentFile <- file.path(vsnDir, "output/Step2_MotifEnrichment.tsv")
file.exists(scenicLoomPath)
file.exists(motifEnrichmentFile)
list.files()




library(SCopeLoomR)
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name="MotifRegulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom)
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
cellClusters <- get_clusterings(loom)

close_loom(loom)

data.dir <- './results10k'
dir.create(data.dir)
setwd(data.dir)


cellClusters <- cellInfo %>% select("ArchR_Clusters","ClusterOneControl","groups")

cellClusters <- cellClusters[str_order((cellClusters$ClusterOneControl), numeric = TRUE),]

selectedClustering<- "ClusterOneControl" # select resolution

# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedClustering]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[ , str_sort(colnames(regulonActivity_byCellType_Scaled), numeric = TRUE)]

row.subsections <- c(15,12,9,15,9,9,6)
order <- str_sort(colnames(regulonActivity_byCellType_Scaled), numeric = TRUE)
options(repr.plot.width=8, repr.plot.height=22) # To set the figure size in Jupyter
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity", column_order = order, #row_order = sort(rownames(regulonActivity_byCellType_Scaled)), #order(colnames(as.matrix(regulonActivity_byCellType_Scaled))
                                   row_names_gp=grid::gpar(fontsize=6))) 
                                   #column_split = data.frame(rep(c("Leukocytes", "B", "C", "D", "E","F","G"), row.subsections)))) # row font size
pdf("Heatmap_Regulon_ClusterOneControl.pdf", width = 15, height = 20)
hm
dev.off()

regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later



###show exact values 14:

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)

viewTable(topRegulators, options = list(pageLength = 10))

write.xlsx2(topRegulators, file = paste0("Top_Regulators_",selectedClustering,".xlsx"))

#Cell-type specific regulators 16/17

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[(str_order(colnames(regulonAUC),numeric = TRUE)), selectedClustering])
## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss, order_rows = TRUE)
plotly::ggplotly(rssPlot$plot)


pdf(paste0("Rss_Plot_Regulon_",selectedClustering,".pdf"), width = 15, height = 20)
rssPlot
dev.off()

rssPlot_table <- as.data.frame(rssPlot[["df"]])
write.xlsx(rssPlot_table, file = paste0("RSSPlot-Table_",selectedClustering,".xlsx"))

#does not want to run in loop or 
n<- list()
n <- c(1:25)
sapply(n, FUN = makeplots)

makeplots <- function(n){
pdf(paste0("Rss_Marker_",selectedClustering,"_C",n,"_x_TNT.pdf"), width = 4, height = 6)  
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
plotRSS_oneSet(rss, setName = paste0("C",n,"_x_TNT")) # cluster ID
dev.off()

pdf(paste0("Rss_Marker_",selectedClustering,"_C",n,"_x_MHC.pdf"), width = 4, height = 6)  
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
plotRSS_oneSet(rss, setName = paste0("C",n,"_x_MHC")) # cluster ID
dev.off()

pdf(paste0("Rss_Marker_",selectedClustering,"_C",n,"_x_CTR.pdf"), width = 4, height = 6)  
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
plotRSS_oneSet(rss, setName = paste0("C",n,"_x_CTR")) # cluster ID
dev.off()
}

