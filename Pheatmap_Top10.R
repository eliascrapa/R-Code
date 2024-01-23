#install.packages("pheatmap")
library("xlsx")
library("tidyverse")
library("magrittr")

# load package
library(pheatmap)
library(RColorBrewer)

setwd("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/")
data.dir <- './Heatmaps'
dir.create(data.dir)
setwd(data.dir)
data.dir <- './10genes'
dir.create(data.dir)
setwd(data.dir)
#Heatmapcomparison MHC and TNT
for (n in 2:2){
  
  

  tnt_input <- paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_TNT_vs_CONTROLS_C",n,".xlsx")
  mhc_input <- paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_MHC_vs_CONTROLS_C",n,".xlsx")
  filem <-  paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsCTR_C",n,".rnk")  
  title <- paste0("Venn Diagram Cluster ",n)

  TNT <- read.xlsx(tnt_input, 1)
  MHC <- read.xlsx(mhc_input, 1)
  
  sapply(TNT, class) 
  sapply(MHC, class) 
  #change columns 3-6 to numeric values!!!!! change if other imputformat
  TNT[,3:6] <- sapply(TNT[,3:6],as.numeric)
  MHC[,3:6] <- sapply(MHC[,3:6],as.numeric)
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
  
  MHC <- MHC[!grepl("^mt-",MHC$names),]
  if (!between(n, 6, 9)) {
    #TNT_logFC_sel <- subset(TNT_logFC_sel, names != (grep(TNT_logFC_sel$names,"^m")))
    MHC <- MHC[!grepl("Tnnt2",MHC$names),]
    MHC <- MHC[!grepl("Hbb-bs",MHC$names),]
    MHC <- MHC[!grepl("Myh6",MHC$names),]
    MHC <- MHC[!grepl("^mt-",MHC$names),]
    MHC <- MHC[!grepl("Ttn",MHC$names),]
  } else {
    
  }
  
  #filter for adj p-value <0.05
  TNT_sel <- dplyr::filter(TNT, p_val_adj < 0.05)
  MHC_sel <- dplyr::filter(MHC, p_val_adj < 0.05)
  #remove other colums in table
  TNT_sel <-dplyr::select(TNT_sel, names, avg_log2FC)#,p_val_adj)
  MHC_sel <-dplyr::select(MHC_sel, names, avg_log2FC)#,p_val_adj)

  #add top 30 and lowest 30 genes, order genes - if loop checks for more then 60 genes because else slice_ produced doubbles
  if (nrow(TNT_sel) >60){
  TNTtop10 <- rbind((TNT_sel %>% slice_min(TNT_sel$avg_log2FC,n = 10)),(TNT_sel %>% slice_max(TNT_sel$avg_log2FC,n = 10))) 
  } else{ TNTtop10 <- TNT_sel}
  TNTtop10 <- TNTtop10[rev(order(TNTtop10$avg_log2FC)),]
  
  if (nrow(MHC_sel) > 60){
  MHCtop10 <- rbind((MHC_sel %>% slice_min(MHC_sel$avg_log2FC,n = 10)),(MHC_sel %>% slice_max(MHC_sel$avg_log2FC,n = 10))) 
  } else{ MHCtop10 <- MHC_sel}
  MHCtop10 <- MHCtop10[rev(order(MHCtop10$avg_log2FC)),]
  
  #remove other rows in original table for following joining, Columnnames are set
  TNT <-dplyr::select(TNT, names, avg_log2FC)#,p_val_adj)
  MHC <-dplyr::select(MHC, names, avg_log2FC)
  
  TNTtop10JMHC <- left_join(TNTtop10, MHC, by = "names")  %>% set_colnames(c("names","TNT","MHC"))
  MHCtop10JTNT <- left_join(MHCtop10, TNT, by = "names") %>% set_colnames(c("names","MHC","TNT"))

  #namesrow is moved to rownames and removed and converted to matrix
  rownames(TNTtop10JMHC) <- TNTtop10JMHC$names
  TNTtop10JMHC <- as.matrix(TNTtop10JMHC[2:3])
  colnames(TNTtop10JMHC) <- c("TnT","MyHC")
  
  rownames(MHCtop10JTNT) <- MHCtop10JTNT$names
  MHCtop10JTNT <- as.matrix(MHCtop10JTNT[2:3])
  colnames(MHCtop10JTNT) <- c("MyHC","TnT")
  

#output
if (nrow(TNTtop10JMHC) >0 ){
pdf(paste0("TNT_Heatmap_Cluster_",n,".pdf"))   #, width = 480, height = 1240)
pheatmap(TNTtop10JMHC, display_numbers = F, cluster_rows = FALSE, cluster_cols = FALSE,scale = "none", cellwidth = 25, cellheight = 18, fontsize_row = 16)
dev.off()
} else {}

if (nrow(MHCtop10JTNT) >0 ){
pdf(paste0("MHC_Heatmap_Cluster_",n,".pdf"))
pheatmap(MHCtop10JTNT, display_numbers = F, cluster_rows = FALSE, cluster_cols = FALSE,scale = "none", cellwidth = 25, cellheight = 18, fontsize_row = 16)
dev.off()
} else {}
 
  
}
