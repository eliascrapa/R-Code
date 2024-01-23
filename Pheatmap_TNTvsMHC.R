#install.packages("pheatmap")
library("xlsx")
library("tidyverse")
library("magrittr")

# load package
library(pheatmap)
library(RColorBrewer)

setwd("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/")

data.dir <- './TNTvsMHC/'
dir.create(data.dir)
setwd(data.dir)

data.dir <- './Heatmaps_TNTvsMHC'
dir.create(data.dir)
setwd(data.dir)


#Heatmapcomparison MHC and TNT
for (n in 1:25){
  
  

  tnt_input <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsMHC/RNA_ClusterMarker_TNT_vs_MHC_C",n,".xlsx")
  filem <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsCTR_C",n,".rnk")  
  title <- paste0("Venn Diagram Cluster ",n)

  TNT <- read.xlsx(tnt_input, 1)

  
  sapply(TNT, class) 
 
  #change columns 3-6 to numeric values!!!!! change if other imputformat
  TNT[,3:6] <- sapply(TNT[,3:6],as.numeric)
 
  
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  
####Filtering out Cardiomyocytegenes in Non-CM, + filtering mito genes  
  TNT <- TNT[!grepl("^mt-",TNT$names),]
  if (!between(n, 6, 9)) {
    #TNT_logFC_sel <- subset(TNT_logFC_sel, names != (grep(TNT_logFC_sel$names,"^m")))
    TNT <- TNT[!grepl("Tnnt2",TNT$names),]
    TNT <- TNT[!grepl("Hbb-bs",TNT$names),]
    TNT <- TNT[!grepl("Myh6",TNT$names),]
    TNT <- TNT[!grepl("^mt-",TNT$names),]
    TNT <- TNT[!grepl("Ttn",TNT$names),]
  } else {
    
  }
  
  
  
  #filter for adj p-value <0.05
  TNT_sel <- dplyr::filter(TNT, p_val_adj < 0.05)

  #remove other colums in table
  TNT_sel <-dplyr::select(TNT_sel, names, avg_log2FC)#,p_val_adj)


  #add top 30 and lowest 30 genes, order genes - if loop checks for more then 60 genes because else slice_ produced doubbles
  if (nrow(TNT_sel) >60){
  TNTtop10 <- rbind((TNT_sel %>% slice_min(TNT_sel$avg_log2FC,n = 30)),(TNT_sel %>% slice_max(TNT_sel$avg_log2FC,n = 30))) 
  } else{ TNTtop10 <- TNT_sel}
  TNTtop10 <- TNTtop10[rev(order(TNTtop10$avg_log2FC)),]
  
  
  rownames(TNTtop10) <- TNTtop10$names
  TNTtop10 <- as.matrix(TNTtop10[2])

#output
if (nrow(TNTtop10) >0 ){

gg <- pheatmap(TNTtop10, display_numbers = F, cluster_rows = FALSE, cluster_cols = FALSE,scale = "none", cellwidth = 25, cellheight = 7.8, fontsize_row = 7.8, main = paste0("Cluster",n," TnT(+) vs MHC(-)"))
ggsave(paste0("TNTvsMHC_Heatmap_Cluster_",n,".pdf"),gg, width=2, height=4, units="in", scale=2)

} else {}


  
}
