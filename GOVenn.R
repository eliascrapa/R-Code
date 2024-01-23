#Make Venn Diagramm
library("xlsx")
library(GOplot)


setwd("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/")
data.dir <- './Venn'
dir.create(data.dir)
setwd(data.dir)
#VEnnDiagramm MHC and TNT
for (n in 1:1){
  
  skip_to_next <- FALSE
  
 
  
  tnt_input <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_TNT_vs_CONTROLS_C",n,".xlsx")
  mhc_input <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_MHC_vs_CONTROLS_C",n,".xlsx")
  filem <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsCTR_C",n,".rnk")  
  title <- paste0("Venn Diagram Cluster ",n)
  vennoutput <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/Venn/Venndiagram_TNTvsMHC_Cluster_",n,".pdf")
  
  TNT <- read.xlsx(tnt_input, 1)
  MHC <- read.xlsx(mhc_input, 1)
  
  sapply(TNT, class) 
  sapply(MHC, class) 
  TNT[,3:6] <- sapply(TNT[,3:6],as.numeric)
  MHC[,3:6] <- sapply(MHC[,3:6],as.numeric)
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  
  TNT <- dplyr::filter(TNT, p_val_adj < 0.05)
  MHC <- dplyr::filter(MHC, p_val_adj < 0.05)
  TNT_sel <-dplyr::select(TNT, names, avg_log2FC)
  MHC_sel <-dplyr::select(MHC, names, avg_log2FC)
  
pdf(vennoutput, width = 960, height = 720)
tryCatch(
print(  
mylist <- GOVenn(TNT_sel,MHC_sel, title = title, label = c("TNT","MHC"), lfc.col = c('red','white','lightblue'), circle.col= c('#ffc000','#1f74ad'),plot = F) # ,, , 
),


error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }  
dev.off()

mylist <- GOVenn(TNT_sel,MHC_sel, title = title, label = c("TNT","MHC"), lfc.col = c('red','white','lightblue'), circle.col= c('#ffc000','#1f74ad'),plot = F) # ,, , 

lapply(mylist[["table"]], write, paste0(title,".txt"), append=TRUE, ncolumns=1000)

}
