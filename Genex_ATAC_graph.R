#Make Venn Diagramm
library("xlsx")
library(GOplot)
library(ggplot2)
library(tidyverse)


setwd("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/")
data.dir <- './GenexCluster'
dir.create(data.dir)
setwd(data.dir)
#VEnnDiagramm MHC and TNT
my_list <-list()
for (n in 23:25){
  
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
 
  
  df <- cbind(length(TNT_sel$avg_log2FC[TNT_sel$avg_log2FC>=0]),length(TNT_sel$avg_log2FC[TNT_sel$avg_log2FC<0]))
  df <- rbind(df,(cbind(length(MHC_sel$avg_log2FC[MHC_sel$avg_log2FC>=0]),length(MHC_sel$avg_log2FC[MHC_sel$avg_log2FC<0]))))
  df <- cbind(c("TNT","MHC"),n,df)
  colnames(df) <- c("Genotype", "Cluster", "Upregulated", "Downregulated")
  
  my_list[[length(my_list) + 1]] <- df
  
}

saveRDS(my_list, file="GeneexpressionPerGEnotypeacrossCluster.RDS")

df2 <- as.data.frame(do.call(rbind, (my_list)))
df2$Upregulated <- as.numeric(df2$Upregulated)
df2$Downregulated <- as.numeric(df2$Downregulated)
df2$Cluster <- as.factor(as.numeric(df2$Cluster))





ggplot(data=df2, aes(x=Cluster, y=Upregulated, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()



## Using your df.m data frame
ggplot(df2, aes(Cluster), ylim(-200:2500)) + 
  geom_bar(data = df2, 
           aes(y = Upregulated, fill = Genotype), stat = "identity", position = "dodge") +
  geom_bar(data = df2, 
           aes(y = -Downregulated, fill = Genotype), stat = "identity", position = "dodge") + 
  geom_hline(yintercept = 0,colour = "grey90")+
  theme_minimal()+
  scale_fill_brewer(palette="Paired")


gg1 <- last_plot() + 
  geom_text(data = df2, 
            aes(Cluster, Upregulated, group=Genotype, label=Upregulated),
            position = position_dodge(width=0.9), vjust = -1.3, size=2.3) +
  geom_text(data = df2, 
            aes(Cluster, -Downregulated, group=Genotype, label=Downregulated),
            position = position_dodge(width=0.9), vjust = +2.2, size=2.3) +
  coord_cartesian(ylim = c(-200, 2500))


gg1
ggsave(file="Genex-PLot.pdf",gg1, width=12, height=9, units="in")


####################################################

ATAC <- as.data.frame(read.xlsx2("//Users/eliascrapa/ArchR/All/Features_ClusterOneControl.xlsx",1))
colnames(ATAC) <- c("Cluster","Genotype","TotalFeatures","Upregulated","Percentage Upregulated","Downregulated","Percentage Downregulated")
ATAC <- ATAC %>% mutate_at('Cluster','TotalFeatures','Upregulated','Percentage Upregulated','Downregulated','Percentage Downregulated', as.numeric) %>% str()

cols.num <- c("Cluster","TotalFeatures","Upregulated","Percentage Upregulated","Downregulated","Percentage Downregulated")
ATAC[cols.num] <- sapply(ATAC[cols.num],as.numeric)
ATAC$Cluster <- as.factor(ATAC$Cluster)

ggplot(ATAC, aes(Cluster), ylim(-95000, 95000)) + 
  geom_bar(data = ATAC, 
           aes(y = Upregulated, fill = Genotype), stat = "identity", position = "dodge") +
  geom_bar(data = ATAC, 
           aes(y = -Downregulated, fill = Genotype), stat = "identity", position = "dodge") + 
  geom_hline(yintercept = 0,colour = "grey90")+
  theme_minimal()+
  scale_fill_brewer(palette="Paired")


gg2 <- last_plot() + 
  geom_text(data = ATAC, 
            aes(Cluster, Upregulated, group=Genotype, label=Upregulated),
            position = position_dodge(width=0.9), vjust = -1.3, size=2.3) +
  geom_text(data = ATAC, 
            aes(Cluster, -Downregulated, group=Genotype, label=Downregulated),
            position = position_dodge(width=0.9), vjust = +2.2, size=2.3) +
  coord_cartesian(ylim = c(-95000, 50000))


gg2
ggsave(file="ATAC-PLot.pdf",gg2, width=12, height=9, units="in")
