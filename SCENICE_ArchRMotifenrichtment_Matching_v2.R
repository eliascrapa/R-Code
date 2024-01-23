#Comparing Cardiac Hypertrophy in IPA Analysis. 

library(xlsx)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)

library("magrittr")

# load package
library(pheatmap)
library(RColorBrewer)
library(fuzzyjoin)




setwd("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/SCENIC/500bp/")
data.dir <- './MotifenrichmentSCENIC/newMATCH-Zbigger0.3_maxdist0.18_3'
dir.create(data.dir)
setwd(data.dir)
tfmatches <- list()

for (n in 1:25){
  
  tabledfTNT <- read.csv2(paste0("/Users/eliascrapa/ArchR/All/additional/CHROMVARdev/Chromvar_MotifsC",n,"_x_TNT.csv"), sep = ",", header = T)
  tabledfMHC <- read.csv2(paste0("/Users/eliascrapa/ArchR/All/additional/CHROMVARdev/Chromvar_MotifsC",n,"_x_MHC.csv"), sep = ",", header = T)
  
  SCENICtable <- read.xlsx(paste0("/Users/eliascrapa/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/SCENIC/500bp/results500bp/Top_Regulators_ClusterOneControl.xlsx"), 1)
  
  #SCENIC Data Manipulation
  #Calculate Z-Score from Relative Activity
  SCENICtable$relactzscore <- SCENICtable$RelativeActivity
  #SCENICtable$doublezscore <- scale(SCENICtable$RelativeActivity)
  SCENICtable <-SCENICtable %>% mutate(TF = Regulon ) # Create duplicate of column with new name
  SCENICtable$TF <- gsub(r"{\s*\([^\)]+\)}", "", SCENICtable$TF) # remove Parcenthesis with number of Genes
  #SCENICtable$Regulon_short <-abbreviate(SCENICtable$Regulon_short, minlength = 5, method = c("left.kept"), use.classes = F)
  SCENICtable <- SCENICtable %>% dplyr::select("CellType","TF", everything()) 
  SCENICtable$TF <- gsub("_.*","", SCENICtable$TF)
  
  # subset 
  SCENICtableTNT <- SCENICtable %>% filter(CellType == paste0("C",n,"_x_TNT"))
  SCENICtableTNT <- SCENICtableTNT %>% mutate(substring = substr(SCENICtableTNT$TF,1,3))
  SCENICtableMHC <- SCENICtable %>% filter(CellType == paste0("C",n,"_x_MHC"))
  SCENICtableMHC <- SCENICtableMHC %>% mutate(substring = substr(SCENICtableMHC$TF,1,3))
  
  #TNT
  tabledfTNT$name<- gsub("_.*","",tabledfTNT$name)
  colnames(tabledfTNT) <- c("TF_number", "group", "group_name", "seqnames", "idx","TF","Log2FC","FDR","zscore")
  tabledfTNT <- tabledfTNT %>% select("TF","zscore","FDR","group_name","idx") 
  tabledfTNT[,2:3] <- lapply(tabledfTNT[,2:3], as.numeric)
  tabledfTNT <- tabledfTNT[order(tabledfTNT$zscore),] 
  tabledfTNT <- tabledfTNT %>% distinct(TF, .keep_all = TRUE)
  tabledfTNT <- tabledfTNT %>% filter(FDR < 0.01)
  
  tabledfTNT <- tabledfTNT %>% dplyr::select("group_name", everything()) 
  tabledfTNT <- tabledfTNT %>% rename("CellType" = "group_name")
  tabledfTNT <- tabledfTNT %>% dplyr::select("CellType","TF", everything()) 
  tabledfTNT <- tabledfTNT %>% mutate(substring = substr(tabledfTNT$TF,1,3))
  
  
  #fuzzy joining to also match slightly differently call TFs
  TNTjoinedstring <-  stringdist_join(tabledfTNT, SCENICtableTNT, 
                  by = "TF",
                  #by.y = "Regulonshort",
                  mode = "left",
                  ignore_case = FALSE, 
                  method = "jw", 
                  #p=0.25,
                  max_dist = 0.2, #this distance determined by manual review of the matches
                  distance_col = "dist") %>%
    group_by(TF.y) %>%
    slice_min(order_by = dist, n = 1)
  
  TNTjoinedstring <- TNTjoinedstring %>% mutate(combzscore = (zscore + relactzscore)/2) # calculate aggregated Z-Score
  
  TNTjoinedstring <- TNTjoinedstring %>% filter(substring.x==substring.y) %>% filter(combzscore >0.3)
  
  TNTjoinedMAT <- TNTjoinedstring %>% select(TF.x,TF.y,combzscore)
  TNTjoinedMAT <- TNTjoinedMAT[rev(order(TNTjoinedMAT$combzscore)),]
  TNTjoinedMAT$TF.y <- as.factor(TNTjoinedMAT$TF.y)
  TNTjoinedMAT$TF.x <- as.factor(TNTjoinedMAT$TF.x)
  
  ggp <- ggplot(TNTjoinedMAT, aes(TF.y,TF.x)) +                           # Create heatmap with ggplot2
    geom_tile(aes(fill = combzscore, y = reorder(TF.y,combzscore),x = reorder(TF.x,combzscore) )) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Enriched Motifs in Chromatin accessibility data") + ylab("Active Regulons")+
    ggtitle(paste0("Aggregated Regulon-/Motifenrichment- Z-Scores for respective TFs TnT - Cluster ",n)) + 
    theme(legend.position = c(.9, .15)) + guides(fill=guide_colourbar(title="Combined Z-Score", direction ="horizontal", barwidth = 10, ticks = FALSE, title.position="top", title.hjust=0.5))  +
    theme(text=element_text(size=17))
  
  
  pdf(paste0("Aggregated Z-Score Motifaccesibility and Regulonactivity_C",n,"_x_TnT.pdf"), width = 14, height = 6)  
  #options(repr.plot.width=14, repr.plot.height=4) # To set the figure size in Jupyter
  print(ggp)                                                               # Print heatmap
  dev.off()
  
  
  #MHC
  tabledfMHC$name<- gsub("_.*","",tabledfMHC$name)
  colnames(tabledfMHC) <- c("TF_number", "group", "group_name", "seqnames", "idx","TF","Log2FC","FDR","zscore")
  tabledfMHC <- tabledfMHC %>% select("TF","zscore","FDR","group_name","idx") 
  tabledfMHC[,2:3] <- lapply(tabledfMHC[,2:3], as.numeric)
  tabledfMHC <- tabledfMHC[order(tabledfMHC$zscore),] 
  tabledfMHC <- tabledfMHC %>% distinct(TF, .keep_all = TRUE)
  tabledfMHC <- tabledfMHC %>% filter(FDR < 0.01)
  
  tabledfMHC <- tabledfMHC %>% dplyr::select("group_name", everything()) 
  tabledfMHC <- tabledfMHC %>% rename("CellType" = "group_name")
  tabledfMHC <- tabledfMHC %>% dplyr::select("CellType","TF", everything()) 
  tabledfMHC  <- tabledfMHC  %>% mutate(substring = substr(tabledfMHC $TF,1,3))
  
  #fuzzy joining to also match slightly differently call TFs
  MHCjoinedstring <-  stringdist_join(tabledfMHC, SCENICtableMHC, 
                                      by = "TF",
                                      #by.y = "Regulonshort",
                                      mode = "left",
                                      ignore_case = FALSE, 
                                      method = "jw", 
                                      max_dist = 0.2, #this distance determined by manual review of the matches
                                      distance_col = "dist") %>%
    group_by(TF.y) %>%
    slice_min(order_by = dist, n = 1)
  
  MHCjoinedstring <- MHCjoinedstring %>% mutate(combzscore = (zscore + relactzscore)/2) # calculate aggregated Z-Score
  MHCjoinedstring <- MHCjoinedstring %>% filter(substring.x==substring.y) %>% filter(combzscore >0.3)
  
  MHCjoinedMAT <- MHCjoinedstring %>% select(TF.x,TF.y,combzscore)
  MHCjoinedMAT <- MHCjoinedMAT[rev(order(MHCjoinedMAT$combzscore)),]
  MHCjoinedMAT$TF.y <- as.factor(MHCjoinedMAT$TF.y)
  MHCjoinedMAT$TF.x <- as.factor(MHCjoinedMAT$TF.x)
  
  ggp <- ggplot(MHCjoinedMAT, aes(TF.y,TF.x)) +                           # Create heatmap with ggplot2
    geom_tile(aes(fill = combzscore, y = reorder(TF.y,combzscore),x = reorder(TF.x,combzscore) )) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Enriched Motifs in Chromatin accessibility data") + ylab("Active Regulons")+
    ggtitle(paste0("Aggregated Regulon-/Motifenrichment- Z-Scores for respective TFs MyHC - Cluster ",n)) +
    theme(legend.position = c(.9, .15)) + guides(fill=guide_colourbar(title="Combined Z-Score", direction ="horizontal", barwidth = 10, ticks = FALSE, title.position="top", title.hjust=0.5))  +
    theme(text=element_text(size=17))
  
pdf(paste0("Aggregated Z-Score Motifaccesibility and Regulonactivity_C",n,"_x_MyHC.pdf"), width = 14, height = 6)  
  #options(repr.plot.width=14, repr.plot.height=4) # To set the figure size in Jupyter
  print(ggp)                                                           # Print heatmap
  dev.off()
  
  
#tfmatches[[n]] <- if_else(nrow(TNTjoinedstring) >0 & nrow(MHCjoinedstring) >0, rbind(TNTjoinedstring,MHCjoinedstring),if_else(nrow(TNTjoinedstring) >0 & nrow(MHCjoinedstring) == 0,TNTjoinedstring,if_else(nrow(TNTjoinedstring) ==0 & nrow(MHCjoinedstring) >0,MHCjoinedstring,NA_real_)))

  
  
}

matchTFs <- sapply(tfmatches,rbind)  
  
  
