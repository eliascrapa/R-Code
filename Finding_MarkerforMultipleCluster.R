library(tidyverse)
library(xlsx)
library(sf)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggpubr)
library(pheatmap)
library(openxlsx)

data.dir <- './ClusterMarker'
dir.create(data.dir)
setwd(data.dir)

load("/Users/eliascrapa/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/recent.RData")

#new.cluster.ids <- c("C1","C2 – Macrophages","C3 – Monocytes","C4 – CD-8+ T-Cells","C5 – B-Cells","C6 – CM_Atrial",
#                     "C7 – CM","C8 – CM","C9 – CM_ventricular","C10 – EpiC","C11 – EpiC","C12 – EpiC",
#                     "C13 – EC","C14 – EC","C15 – EC","C16 – EC","C17 – ET_Lymphatic","C18 – VSMC",
#                     "C19 – VSMC","C20 – VSMC","C21 - Fib","C22 - Fib","C23 - Fib","C24","C25")

#Add grouped name to Seurat Object
new.cluster.ids <- c("unsp","Macrophages-Monocytes","Macrophages-Monocytes","T-Cells","B-Cells","Cardiomyocyte",
                     "Cardiomyocyte","Cardiomyocyte","Cardiomyocyte","Epicardial Cells","Epicardial Cells","Epicardial Cells",
                     "Endothelial Cells","Endothelial Cells","unsp","Endothelial Cells","Endothelial Cells","VSMC",
                     "VSMC","VSMC","Fibroblasts","Fibroblasts","Fibroblasts","unsp","unsp") # 
                     
Idents(object = all) <- "ArchR_Clusters"
my_levels <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25")
Idents(all) <- factor(Idents(all), levels= my_levels)
names(new.cluster.ids) <- levels(all)
all <- RenameIdents(all, new.cluster.ids)
all@meta.data[["CellsGrouped"]]<- Idents(object = all)
DimPlot(all, reduction = "harmony_archr", label = TRUE, pt.size = 0.4, label.size = 6) + NoLegend()
Idents(object = all) <- "ArchR_Clusters"


cellgrouping <- c("ArchR_Clusters","CellsGrouped")
markers_combined_list <- list()
topgenes_binom_up_list <- list()
markers_up_all_maxp_list <- list()
markers_up_lfc1 <- list()

for (c in cellgrouping){
  
  ggsave(filename = paste0(c," - UMAP Cellgroups.pdf"),
         plotReducedDim(
           all_sce,
           'HARMONY_ARCHR',
           ncomponents = 2,
           percentVar = NULL,
           colour_by = c
         ))
  
  Idents(object = all) <- c
  
  all_sce <- as.SingleCellExperiment(all, assay = "RNA")
  
  for (i in 1:length(levels(as.factor(all_sce@colData@listData[[c]])))) {

    #save overview UMAP plot with respective cellgrouping



    # testing for the alternative hypothesis that LFC > 1
    markers_up_lfc1 <- findMarkers(
      all_sce, 
      groups = factor(all_sce@colData@listData[[c]]), # clusters to compare
      block = all_sce@colData@listData[["orig.ident"]],    # covariates in statistical model
      test.type = "t",   # t-test (default)
      direction = "up", # test for up-regulated genes only
      lfc = 1, # null hypothesis log-fold-change = 1
      pval.type = "any" # ranking of p-values based on any comparison (default)
    )
    
    cluster_markers_up_lfc1 <- markers_up_lfc1@listData[[i]]
    cluster_markers_up_lfc1[cluster_markers_up_lfc1$Top <= 3, ]

    
    markers_up_lfc1@listData <-       lapply(markers_up_lfc1@listData, function(x) {
      x <- cbind(row.names(x),x)
      
    }) 
    
    write.xlsx(markers_up_lfc1@listData, file = paste0(c," - ",i," - markers_up_lfc1.xlsx"))
    
    # Heatmaps ----
    
    # select some top genes for cluster 8
    cluster_top10 <- cluster_markers_up_lfc1[cluster_markers_up_lfc1$Top <= 10, ]
    
    # heatmap of expression values
    ggsave(filename = paste0(c," - ",names(markers_up_lfc1@listData[i])," - Heatmap of expression values.pdf"),
   # pdf(file=paste0(c," - ",names(markers_up_lfc1@listData[i])," - Heatmap of expression values.pdf"))
           plotHeatmap(all_sce, 
                features = rownames(cluster_top10),
                order_columns_by = c("ArchR_Clusters","CellsGrouped", "orig.ident"))
  )#  dev.off()
    
   
# ranking based on the maximum p-value across all pairwise comparisons
markers_up_all_maxp_list <- findMarkers(
  all_sce, 
  groups = factor(all_sce@colData@listData[[c]]), # clusters to compare
  block = all_sce@colData@listData[["orig.ident"]],    # covariates in statistical model
  test.type = "t",   # t-test (default)
  direction = "up", # test for up-regulated genes only
  lfc = 0, # null hypothesis log-fold-change = 1
  pval.type = "all" # ranking of p-values based on all comparisons
)



# fetching top markers for cluster 8
cluster_markers_up_all <- markers_up_all_maxp_list@listData[[i]]
cluster_markers_up_all[1:5,]

clustermarker <- rownames(cluster_markers_up_all[1:5,])


#Plot Genes 
p1 <- plotExpression(all_sce, features = clustermarker[1], x = "ident") + theme(axis.text.x = element_text(angle = 45, 
                                                                                                     hjust = 1))                                                                                             
p2 <- plotReducedDim(
  all_sce,
  'HARMONY_ARCHR',
  ncomponents = 2,
  percentVar = NULL,
  colour_by = clustermarker[1],
  by_exprs_values = "logcounts"
)

p3 <- plotExpression(all_sce, features = clustermarker[2], x = "ident") + theme(axis.text.x = element_text(angle = 45, 
                                                                                                           hjust = 1))                                                                                             
p4 <- plotReducedDim(
  all_sce,
  'HARMONY_ARCHR',
  ncomponents = 2,
  percentVar = NULL,
  colour_by = clustermarker[2],
  by_exprs_values = "logcounts"
)

p5 <- plotExpression(all_sce, features = clustermarker[3], x = "ident") + theme(axis.text.x = element_text(angle = 45, 
                                                                                                           hjust = 1))                                                                                             
p6 <- plotReducedDim(
  all_sce,
  'HARMONY_ARCHR',
  ncomponents = 2,
  percentVar = NULL,
  colour_by = clustermarker[3],
  by_exprs_values = "logcounts"
)

p7 <- plotExpression(all_sce, features = clustermarker[4], x = "ident") + theme(axis.text.x = element_text(angle = 45, 
                                                                                                           hjust = 1))                                                                                             
p8 <- plotReducedDim(
  all_sce,
  'HARMONY_ARCHR',
  ncomponents = 2,
  percentVar = NULL,
  colour_by = clustermarker[4],
  by_exprs_values = "logcounts"
)

p9 <- plotExpression(all_sce, features = clustermarker[5], x = "ident") + theme(axis.text.x = element_text(angle = 45, 
                                                                                                           hjust = 1))                                                                                             
p10 <- plotReducedDim(
  all_sce,
  'HARMONY_ARCHR',
  ncomponents = 2,
  percentVar = NULL,
  colour_by = clustermarker[5],
  by_exprs_values = "logcounts"
)

ggsave(filename = paste0(c," - ",names(markers_up_lfc1@listData[i])," - Expression Plots Top5.pdf"),
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol=4, nrow=3),width = 16,
height = 10
)




# Alternative testing strategies ----
# Wilcoxon rank-sum test
markers_wilcox_up <- findMarkers(
  all_sce, 
  groups = factor(all_sce@colData@listData[[c]]), # clusters to compare
  block = all_sce@colData@listData[["orig.ident"]],    # covariates in statistical model
  test.type = "wilcox",   # t-test (default)
  direction = "up"
)

vsmc_markers_wilcox_up <- markers_wilcox_up[[i]]
vsmc_markers_wilcox_up[1:10,]
head(vsmc_markers_wilcox_up)


# make a heatmap of AUC values
# we use a custom colour palette that diverges around 0.5
# we optionally do not cluster rows to keep genes in their ranking order
ggsave(filename = paste0(c," - ",names(markers_wilcox_up@listData[i])," - Heatmap Wilcox Rank AUC.pdf"),
       pheatmap(vsmc_markers_wilcox_up[vsmc_markers_wilcox_up$Top <= 6, 4:12],
                breaks = seq(0.7, 1, length.out = 21),
                color = viridis::cividis(21), 
                cluster_rows = FALSE)
)


# Binomial test of proportions
markers_binom_up <- findMarkers(
  all_sce, 
  groups = factor(all_sce@colData@listData[[c]]), # clusters to compare
  block = all_sce@colData@listData[["orig.ident"]],    # covariates in statistical model
  test.type = "binom",   # t-test (default)
  direction = "up"
)

# make a heatmap of expression values for top genes in cluster 8
cluster_markers_binom_up <- markers_binom_up[[i]]
ggsave(filename = paste0(c," - ",names(markers_binom_up@listData[i])," - Heatmap Binomial Up Genes.pdf"),
       plotExpression(all_sce, 
               x = "CellsGrouped",
               features = rownames(cluster_markers_binom_up)[1:4]), width = 18,
       height = 10,
)

# Combining multiple tests ----

markers_combined <- multiMarkerStats(
  t = findMarkers(
    all_sce, 
    groups = factor(all_sce@colData@listData[[c]]), # clusters to compare
    block = all_sce@colData@listData[["orig.ident"]],    # covariates in statistical model
    direction = "up"
    
  ),
  wilcox = findMarkers(
    all_sce, 
    groups = factor(all_sce@colData@listData[[c]]), # clusters to compare
    block = all_sce@colData@listData[["orig.ident"]],    # covariates in statistical model
    test = "wilcox",
    direction = "up"
  ),
  binom = findMarkers(
    all_sce, 
    groups = factor(all_sce@colData@listData[[c]]), # clusters to compare
    block = all_sce@colData@listData[["orig.ident"]],    # covariates in statistical model
    test = "binom",
    direction = "up"
  )
)



markers_up_all_maxp_list@listData <-       lapply(markers_up_all_maxp_list@listData, function(x) {
  x <- cbind(row.names(x),x)
  }) 
openxlsx::write.xlsx(markers_up_all_maxp_list@listData, file = paste0(c," - ",i," - markers_up_all_maxp_list.xlsx"), overwrite = TRUE, rowNames = TRUE)


     
       markers_combined@listData <-       lapply(markers_combined@listData, function(x) {
           x <- cbind(row.names(x),x)
         
       }) 
    
      
       
# the first few rows and columns of the combined results table
write.xlsx(markers_combined@listData, file = paste0(c," - ",i," - Multiple-Tests_combined.xlsx"), overwrite = TRUE, rowNames = TRUE)
print(paste0("wrote combinedtable ",c," - ","C",i))


}

}
