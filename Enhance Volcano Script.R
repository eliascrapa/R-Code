# VolcanoPlotscript

#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
gc()

library(tidyverse)
options(java.parameters = "-Xmx8000m")

library(EnhancedVolcano)
library(magrittr)
library("xlsx")



setwd("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/")

data.dir <- './VP'
dir.create(data.dir)
setwd(data.dir)

geno <- c('TNT','MHC')

for (i in geno) {
#Volcanoplot MHC and TNT
for (n in 7:25){
  
  skip_to_next <- FALSE
  
  
  
  input <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_",i,"_vs_CONTROLS_C",n,".xlsx")
 
  
  title <- paste0("Venn Diagram Cluster ",n)
  Volcano <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/VP/",i,"_Cluster_",n,".png")
  
  Cluster <- read.xlsx(input, 1)

  
  sapply(Cluster, class) 

  Cluster[,3:6] <- sapply(Cluster[,3:6],as.numeric)

  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  
 #####Cluster <- dplyr::filter(Cluster, p_val_adj < 0.05)
 
  #MHC <- dplyr::filter(MHC, p_val_adj < 0.05)
 #TNT_sel <-dplyr::select(TNT, names, avg_log2FC)
 #MHC_sel <-dplyr::select(MHC, names, avg_log2FC)

 if (!between(n, 6, 9)) {
   Cluster <- subset(Cluster, names != "Tnnt2")
   Cluster <- subset(Cluster, names != "Hbb-bs")
   Cluster <- subset(Cluster, names != "Myh6")
   #Cluster <- subset(Cluster, names != (grep(Cluster$names,"^m")))
   Cluster <- Cluster[!grepl("^m",Cluster$names),]
   Cluster <- Cluster[!grepl("Ttn",Cluster$names),]
 } else {
   
 }


 
  rownames(Cluster) <- Cluster$names

print(Cluster)

  
# png device
png(Volcano)
tryCatch(
print(EnhancedVolcano(Cluster,
                      lab = rownames(Cluster),
                      x = 'avg_log2FC',
                      y = 'p_val',
                      pCutoff = 10e-6,
                      FCcutoff = 0.5,
                      caption = 'FC cutoff, 0.5; p-value cutoff, 10e-8',
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      #typeConnectors = "closed",
                      #endsConnectors = "first",
                      #lengthConnectors = unit(5, "npc")
                      )
      
      #xlim = -1.5,
      #ylim = 0.1
      
      
),

error = function(e) { skip_to_next <<- TRUE})
  
if(skip_to_next) { next }  
  
# Close device
dev.off()


}


}


################################################################################


suppressPackageStartupMessages({
  library(Seurat)
  library(venn)
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(enrichR)
  library(rafalib)
})




Idents(object = all) <- all@meta.data$ClusterOneControl
table(all@active.ident)

# plot this clustering
plot_grid(ncol = 3, DimPlot(all, label = T) + NoAxes(), DimPlot(all, group.by = "orig.ident") + 
            NoAxes(), DimPlot(all, group.by = "type") + NoAxes())


# Compute differentiall expression
markers_genes <- FindAllMarkers(alldata, logfc.threshold = 0.2, test.use = "DESeq2", 
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, 
                                assay = "RNA")
