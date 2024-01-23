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
library(tidyverse)
library(ggplot2)

setwd("/Users/eliascrapa/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/SCENIC/500bp/")
vsnDir <- "."


scenicLoomPath <- file.path(vsnDir, "output/500bpscenic.loom")
motifEnrichmentFile <- file.path(vsnDir, "output/Step2_MotifEnrichment.tsv")
file.exists(scenicLoomPath)
file.exists(motifEnrichmentFile)
list.files()




library(SCopeLoomR)
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
#exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name="MotifRegulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom)
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
cellClusters <- get_clusterings(loom)

#motifenrichmentfile
motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=1)[-3,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")

close_loom(loom)

data.dir <- './results500bp'
dir.create(data.dir)
setwd(data.dir)

load("/Users/eliascrapa/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/SCENIC/500bp/afterSCENICs.RData")

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

#add row_order for alphabetical Rows in Heatmap
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity", #for Alphabetic -> column_order = order, row_order = sort(rownames(regulonActivity_byCellType_Scaled)), #order(colnames(as.matrix(regulonActivity_byCellType_Scaled))
                                   row_names_gp=grid::gpar(fontsize=6))) 
                                   #column_split = data.frame(rep(c("Leukocytes", "B", "C", "D", "E","F","G"), row.subsections)))) # row font size
pdf(paste0("Heatmap_Regulon_",selectedClustering,"_alphabetic.pdf"), width = 15, height = 20)
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

write.xlsx2(topRegulators, file = paste0("11Top_Regulators_",selectedClustering,".xlsx"))

#Show all regulators based on Clusters
for (n in 1:25){

CX_topRegulators <- topRegulators %>% filter(str_detect(CellType, paste0("C",n)))
CX_topRegulators$Regulon <- as.character(CX_topRegulators$Regulon)
CX_topRegulators <-CX_topRegulators %>%  arrange((Regulon))
CX_topRegulators$Regulon <- as.factor(CX_topRegulators$Regulon)
CX_topRegulators <- CX_topRegulators[rev(order(CX_topRegulators$RelativeActivity)),]
CX_topRegulators$zscore <- scale(CX_topRegulators$RelativeActivity)
CX_topRegulators$Regulon <- gsub(r"{\s*\([^\)]+\)}", "", CX_topRegulators$Regulon)
CX_topRegulators$Regulon <- gsub("_.*","", CX_topRegulators$Regulon)

CX_topRegulators <- CX_topRegulators %>% filter(CellType == paste("") %>% mutate(CellType = case_when(CellType == paste0("C",n,"_x_TNT") ~ paste0("TnT C",n," – Macrophage"),
                                                                     CellType == paste0("C",n,"_x_MHC") ~ paste0("MyHC C",n," – Macrophage"),
                                                                     CellType == paste0("C",n,"_x_CTR") ~ paste0("CTR C",n," – Macrophage"))))

CXregplot <- ggplot(CX_topRegulators, aes(fill=CellType, y=RelativeActivity, x= reorder(forcats::fct_rev(Regulon), RelativeActivity))) + 
  geom_bar(position=position_dodge(width = 0.5), stat="identity") +
  coord_flip() +
  xlab("Active Regulons") + ylab("Relative Regulon Activity") +
  ggtitle(paste0("Regulon Activity Comparison")) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(fill = "Groups") + theme(legend.title.align=0.5)

ggsave(filename = paste0(selectedClustering,"C",n,"chosenRegulons_overlapping.pdf"), CXregplot, device = "pdf", width = 2, height = 3 )

}



#Chose certain TF an compare them
##################################

data.dir <- './ClusterOneControl_chosenTF'
dir.create(data.dir)
setwd(data.dir)

n =8L
name <- "CM"
#chosenTFs <- c("Gabpa", "Atf6", "Bach1", "Smarca4", "Jun", "Smarcc2") #Cluster 9
#chosenTFs <- c("Foxo1", "Pbx3", "Msx1", "Klf4", "Klf2", "Bcl6b", "Tcf4", "Irf1") #Cluster 14 
#chosenTFs <- c("Gata2", "Glix3", "Klf10", "Foxp1", "Sp4", "Klf2", "Klf4") #Cluster 16 
#chosenTFs <- c("Pbx3", "Lhx2", "Srebf1", "Foxl1", "Klf12","Otx2", "Msx1") #Cluster 20 
#chosenTFs <- c("Atf5", "Foxp2", "Srebf1", "Creb3l2", "Egr2", "Foxo1", "Srebf2") #Cluster 22 
#chosenTFs <- c("Pbx3", "Irf2", "Irf5", "Irf4", "Rel", "Stat2", "Irf9", "Mitf","Nyfa","Nr3c1") #Cluster 2
#chosenTFs <- c("Smarcc1", "Smarcc2","Smarca4","Bach1","Jun","Gabpa","Atf6", "Nfe2l2","Nr3c1")#, "Jund", "Bach2", "Fosb", "Junb", "Fosl1","Nfe2","Nfe2l2") #C9
#chosenTFs <- c("Cebpd","Cebpb","Cebpg","Foxp2","Nr3c1","Pbx1","Gli3","Tef") #C21
#chosenTFs <- c("Smarcc1", "Smarcc2","Smarca4","Bach1","Jun","Atf6", "Nfe2l2","Nr3c1","Pbx1","Foxp1","Fos","Fosb") # C8

{


  
  CX_topRegulators <- topRegulators %>% filter(str_detect(CellType, paste0("C",n,"_")))
  CX_topRegulators$Regulon <- as.character(CX_topRegulators$Regulon)
  CX_topRegulators <-CX_topRegulators %>%  arrange((Regulon))
  CX_topRegulators$Regulon <- as.factor(CX_topRegulators$Regulon)
  CX_topRegulators <- CX_topRegulators[rev(order(CX_topRegulators$RelativeActivity)),]
  CX_topRegulators$zscore <- scale(CX_topRegulators$RelativeActivity)
  CX_topRegulators$Regulon <- gsub(r"{\s*\([^\)]+\)}", "", CX_topRegulators$Regulon)
  CX_topRegulators$Regulon <- gsub("_.*","", CX_topRegulators$Regulon)
  CX_choRegulators <- CX_topRegulators %>% filter (Regulon %in% chosenTFs)
  
  CX_choRegulators <- CX_choRegulators %>% mutate(CellType = case_when(CellType == paste0("C",n,"_x_TNT") ~ paste0("TnT C",n," - ",name),
                                                                       CellType == paste0("C",n,"_x_MHC") ~ paste0("MyHC C",n," - ",name),
                                                                       CellType == paste0("C",n,"_x_CTR") ~ paste0("CTR C",n," - ",name)))
  
  CXregplot <- CX_choRegulators %>%
    mutate(CellType = fct_relevel(CellType, paste0("CTR C",n," - ",name), paste0("TnT C",n," - ",name), paste0("MyHC C",n," - ",name))) %>%
    ggplot( aes(fill=CellType, y=RelativeActivity, x= reorder(forcats::fct_rev(Regulon), RelativeActivity))) + 
    geom_bar(position=position_dodge(width = 0.5), stat="identity") +
    coord_flip() +
    #scale_fill_brewer(palette="Set1") +
    #scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99"))+ 
    scale_fill_manual(values=c("#CC9900", "#006666", "#006699"))+ 
    xlab("Active Regulons") + ylab("Relative Regulon Activity") +
    ggtitle(paste0("Regulon Activity Comparison")) +
    theme(plot.title = element_text(hjust = 0.5)) + labs(fill = "Groups") + theme(legend.title.align=0.5, legend.text = element_text(size=8)) + #, legend.key.width = unit(1, 'in'))  
    theme(legend.position = c(0.86, 0.15))
  ggsave(filename = paste0(selectedClustering,"_C",n,"_TF-",paste0(chosenTFs, collapse ="-"),"_overlapping.pdf"), CXregplot, device = "pdf", width = 6, height = 5)
  file <- paste0(selectedClustering,"_C",n,"_TF-",paste0(chosenTFs, collapse ="-"),"_overlapping.pdf")
  system2('open', args = c('-a Preview.app', file), wait = FALSE)

}

setwd('..')


viewTable(topRegulators, options = list(pageLength = 10))

write.xlsx2(topRegulators, file = paste0("Top_Regulators_",selectedClustering,".xlsx"))

#Cell-type specific regulators 16/17

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[(str_order(colnames(regulonAUC),numeric = TRUE)), selectedClustering])
## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)


pdf(paste0("Rss_Plot_Regulon_",selectedClustering,".pdf"), width = 15, height = 20)
rssPlot
dev.off()

rssPlot_table <- as.data.frame(rssPlot[["df"]])
write.xlsx(rssPlot_table, file = paste0("RSSPlot-Table_",selectedClustering,".xlsx"))


for (n in 1:25){
  
  CX_RssRegulators <- rssPlot_table %>% filter(str_detect(cellType, paste0("C",n)))
  CX_RssRegulators$Topic <- as.character(CX_RssRegulators$Topic)
  CX_RssRegulators <-CX_RssRegulators %>%  arrange((Topic))
  CX_RssRegulators$Topic <- as.factor(CX_RssRegulators$Topic)
  
  
  CXregplot <- ggplot(CX_RssRegulators, aes(fill=cellType, y=Z, x= forcats::fct_rev(Topic))) + 
    geom_bar(position=position_dodge(width = 0.5), stat="identity") +
    coord_flip() 
  
  ggsave(filename = paste0("./RSSRegulatorComparisonplot/C",n,"RSSregPlot_overlapping.pdf"), CXregplot, device = "pdf", width = 6, height = 22 )
  
}

#does not want to run in loop or 
n<- list()
n <- c(1:25)
sapply(n, FUN = makeplotsArchR)

makeplots <- function(n) {
pdf(paste0("Rss_Marker_",selectedClustering,"_C",n,"_x_TNT.pdf"), width = 4, height = 6)  
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
print(plotRSS_oneSet(rss, setName = paste0("C",n,"_x_TNT"))) # cluster ID
dev.off()

pdf(paste0("Rss_Marker_",selectedClustering,"_C",n,"_x_MHC.pdf"), width = 4, height = 6)  
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
print(plotRSS_oneSet(rss, setName = paste0("C",n,"_x_MHC"))) # cluster ID
dev.off()

pdf(paste0("Rss_Marker_",selectedClustering,"_C",n,"_x_CTR.pdf"), width = 4, height = 6)  
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
print(plotRSS_oneSet(rss, setName = paste0("C",n,"_x_CTR"))) # cluster ID
dev.off()
}

makeplotsArchR <- function(n) {
pdf(paste0("./RSS_Marker/Rss_Marker_",selectedClustering,"C",n,".pdf"), width = 4, height = 6)  
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
print(plotRSS_oneSet(rss, setName = paste0("C",n))) # cluster ID
dev.off()
}

## Explore Regulons 32

length(regulons)
sum(lengths(regulons)>=10)
viewTable(cbind(nGenes=lengths(regulons)), options=list(pageLength=10))

numberofgenesinregulon <- cbind(nGenes=lengths(regulons))
write.xlsx(numberofgenesinregulon, file = "NumberofGenesperRegulon500bp.xlsx")
grep("Bach1", names(regulons), value=T) # paste0("^","EBF1","_")
Smarcc2REgulonGenes <- as.data.frame(str_sort(regulons[["Smarcc2_extended"]]))
Bach1REgulonGenes <- as.data.frame(str_sort(regulons[["Bach1_extended"]]))

JunREgulonGenes<- as.data.frame(str_sort(regulons[["Jun"]]))
#TNT highest
gene <- "Slc39a11"
gene <- "Hbb-bs"
gene <- "Hbbt2"
gene <- "Tnnt2"

gene <- c("Slc39a11","Tnnt2")
gene <- "Prkg1"
gene <- "Dmd"
gene <- "Prkg1"
gene <- "Sorbs2"
gene <- "Nppa"
#MHC highest
gene <- "Sorbs2"
gene <- "Pdlim5"
gene <- "Celf2"
gene <- "Dmd"
gene <- "Ctnna3"
gene <- "Pcdh7"
gene <- "Slc8a1"
gene <- "Myh7"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]

#UpGenes <- read.xlsx(paste0("/Users/eliascrapa/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_MHC_vs_CONTROLS_C",n,".xlsx"), 1)
UpGenes <- read.xlsx(paste0("/Users/eliascrapa/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsMHC/RNA_ClusterMarker_TNT_vs_MHC_C",n,".xlsx"), 1)


UpGenes <- UpGenes[rev(order(UpGenes$avg_log2FC)),]
UpGenes <- UpGenes %>% filter(p_val_adj < 0.01)

genes <- UpGenes  %>% slice_max(n=3, order_by=avg_log2FC) %>% select(names)
genes <- unlist(genes)

dim(regulons_incidMat)
genes <- c("Slc39a11", "Tnnt2")#, "Hbb-bs","Prkg1","Dmd","Sorbs2","Nppa") 
incidMat_subset <- regulons_incidMat[,genes]
incidMat_subset <- incidMat_subset[rowSums(incidMat_subset)>25000,]

incidMat_subset



rowSums(incidMat_subset)
order(incidMat_subset)

library(matrixStats)
rowCounts(incidMat_subset, value=TRUE, na.rm=FALSE)

mat <- incidMat_subset[order(rowCounts(incidMat_subset, value=TRUE, na.rm=FALSE)),]

mat[rowSums(mat = TRUE) >= 2, ]

incidMat_subsetmat[rowSums(mat != "") = TRUE, ]

rowCounts(incidMat_subset, value=TRUE, na.rm=FALSE)

#41
tableSubset <- motifEnrichment[Ar=="Smarcc2"]
viewMotifs(tableSubset, colsToShow = c("logo", "NES", "TF" ,"Annotation"), options=list(pageLength=5)) 
