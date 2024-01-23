library(ggplot2)
library(xlsx)
library(ggpubr)
library(tidyverse)

setwd("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/")
data.dir <- './Markermotifs_Chromvar_norm'
dir.create(data.dir)
setwd(data.dir)
  

for (n in 1:25){
  
  tabledfTNT <- read.csv2(paste0("/Users/eliascrapa/ArchR/All/additional/CHROMVARdev/Chromvar_MotifsC",n,"_x_TNT.csv"), sep = ",", header = T)
  tabledfMHC <- read.csv2(paste0("/Users/eliascrapa/ArchR/All/additional/CHROMVARdev/Chromvar_MotifsC",n,"_x_MHC.csv"), sep = ",", header = T)
  
  #TNT
  tabledfTNT$name<- gsub("_.*","",tabledfTNT$name)
  colnames(tabledfTNT) <- c("TF_number", "group", "group_name", "seqnames", "idx","TF","Log2FC","FDR","zscore")
  tabledfTNT <- tabledfTNT %>% select("TF","zscore","FDR","group_name","idx") 
  tabledfTNT[,2:3] <- lapply(tabledfTNT[,2:3], as.numeric)
  tabledfTNT <- tabledfTNT[order(tabledfTNT$zscore),] 
  tabledfTNT <- tabledfTNT %>% distinct(TF, .keep_all = TRUE)
  tabledfTNT <- tabledfTNT %>% slice_max(tabledfTNT$'zscore',n = 30, with_ties =F)
  tabledfTNT <- tabledfTNT %>% drop_na() 

  
  
  #MHC
  tabledfMHC$name<- gsub("_.*","",tabledfMHC$name)
  colnames(tabledfMHC) <- c("TF_number", "group", "group_name", "seqnames", "idx","TF","Log2FC","FDR","zscore")
  tabledfMHC <- tabledfMHC %>% select("TF","zscore","FDR","group_name","idx") 
  tabledfMHC[,2:3] <- lapply(tabledfMHC[,2:3], as.numeric)
  tabledfMHC <- tabledfMHC[order(tabledfMHC$zscore),] 
  tabledfMHC <- tabledfMHC %>% distinct(TF, .keep_all = TRUE)
  tabledfMHC <- tabledfMHC %>% slice_max(tabledfMHC$'zscore',n = 30, with_ties =F)
  tabledfMHC <- tabledfMHC %>% drop_na()
  
  tabledfTNT$same<- ifelse(sapply(tabledfTNT$TF, `%in%`, (dplyr::intersect(tabledfTNT$TF,tabledfMHC$TF))),"matching", "different")
  tabledfMHC$same<- ifelse(sapply(tabledfMHC$TF, `%in%`, (dplyr::intersect(tabledfTNT$TF,tabledfMHC$TF))),"matching", "different")
  
  
  pdf(paste0("MotifMarkers_C",n,"_MHCvsallCTR.pdf"))#, width = 1200, height = 1800)
  print(ggbarplot(tabledfMHC, x = "TF", y = "zscore",
                  fill = "same", # change fill c0olor by mpg_level
                  color = "white",            # Set bar border colors to white
                  palette = "jama",            # jco journal color palett. see ?ggpar
                  sort.val = "asc",          # Sort the value in descending order
                  sort.by.groups = F,     # Don't sort inside each group
                  x.text.angle = 90,          # Rotate vertically x axis texts
                  ylab = "z-score",
                  legend.title = paste0("C",n,"_MyHc"),
                  rotate = T,
                  ggtheme = theme_minimal(),
                  font.tickslab  = 11,
                  position = position_stack(reverse = TRUE),
                  xlim = c(0,3)
                  
                  
  ))
  dev.off()
  
  pdf(paste0("MotifMarkers_C",n,"_TNTvsallCTR.pdf"))#, width = 1200, height = 1800)
  print(ggbarplot(tabledfTNT, x = "TF", y = "zscore",
                  fill = "same",             # change fill c0olor by mpg_level
                  (reverse = TRUE), 
                  color = "white",            # Set bar border colors to white
                  palette = "jama",            # jco journal color palett. see ?ggpar
                  sort.val = "asc",          # Sort the value in descending order
                  sort.by.groups = F,     # Don't sort inside each group
                  x.text.angle = 90,          # Rotate vertically x axis texts
                  ylab = "z-score",
                  legend.title = paste0("C",n,"_TnT"),
                  rotate = T,
                  ggtheme = theme_minimal(),
                  font.tickslab  = 11,
                  position = position_stack(reverse = TRUE),
                  xlim = c(0,3)
                  
  ))
  dev.off()
  
  # }
  
}

#-------------------------------------------------------------------------





`%notin%` <- function(x,y) !(x %in% y) 
`%isin%` <- function(x,y) (x %in% y) 

cluster <- list()
cluster <- c("TNT_C2","TNT_C6","TNT_C9","TNT_C14","TNT_C16","TNT_C19","TNT_C20","TNT_C21","TNT_C22","MHC_C2","MHC_C6","MHC_C9","MHC_C14","MHC_C16","MHC_C20","MHC_C22")
cluster <- c("C6_x_MHC","C9_x_MHC","C14_x_MHC","C16_x_MHC","C20_x_MHC","C21_x_MHC","C22_x_MHC","C6_x_TNT","C9_x_TNT","C14_x_TNT","C16_x_TNT","C20_x_TNT","C21_x_TNT","C22_x_TNT")
cluster <- c("C9_x_MHC")


genotype <- list()
genotype <- c("TNT","MHC")

for (q in genotype){
  
  for (n in 1:25){
    
    tabledf <- read.csv2(paste0("/Users/eliascrapa/ArchR/All/additional/motifs/Table_Markermotifs_up_Enriched_C",n,"_x_",q,"vs_allCTR.csv"), sep = ",", header = T)
    
    #names(tabledf) <- as.matrix(tabledf[1, ])
    #tabledf <- tabledf[-1, ]
    #tabledf[] <- lapply(tabledf, function(x) type.convert(as.character(x)))
    tabledf$TF<- gsub("_.*","",tabledf$TF)
    colnames(tabledf) <- c("Number", "TF", "mlog10Padj", "rank")
    #IF (nrow(tabledf %>% filter(tabledf$mlog10Padj > 1.3)) = 0){ next
    #}ELSE {
    tabledf <- tabledf %>% filter(tabledf$mlog10Padj > 1.3) 
    tabledf <- tabledf[order(tabledf$rank),]
    #tabledf2 <- head(tabledf, 40)
    
    
    tabledf <- tabledf %>% slice_min(tabledf$rank,n = 30)
    tabledf$mlog10Padj <- as.numeric(tabledf$mlog10Padj)
    tabledf <- tabledf[order(-tabledf$rank),]
    #tabledf2$posneg <- ifelse(tabledf2$zscore > 0,"pos", "neg")
    
    pdf(paste0("MotifMarkers_C",n,"_",q,"vsallCTR.pdf"))#, width = 1200, height = 1800)
    print(ggbarplot(tabledf, x = "TF", y = "mlog10Padj",
                    fill = "#EFC000FF",           # change fill color by mpg_level
                    color = "white",            # Set bar border colors to white
                    palette = "jco",            # jco journal color palett. see ?ggpar
                    #sort.val = "desc",          # Sort the value in descending order
                    sort.by.groups = F,     # Don't sort inside each group
                    x.text.angle = 90,          # Rotate vertically x axis texts
                    ylab = "-log10 adj P-Value",
                    legend.title = paste0("C",n,"_",q),
                    rotate = T,
                    ggtheme = theme_minimal(),
                    font.tickslab  = 11
                    
    ))
    dev.off()
    
    
    # }
    
  }
}





ggpar(ggbarplot(tabledfMHC, x = "TF", y = "zscore",
          fill = "same", # change fill c0olor by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jama",            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = F,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "z-score",
          legend.title = paste0("C",n,"_MyHc"),
          rotate = T,
          ggtheme = theme_minimal(),
          font.tickslab  = 11,
          position = position_stack(reverse = TRUE)
          ), ylim = c(0,3), orientation = c("horizontal")
          

)

p1 <- (ggplot(tabledfMHC, aes(x = reorder(TF, zscore), y = zscore))+
  geom_col(aes(fill = same), width = 0.7)) + coord_flip() +
  scale_y_continuous(name = "Z- Score", limits = c(0,2)) + # breaks, labels,  trans)+ 
  xlab("TF") 

p2 <- (ggplot(tabledfTNT, aes(x = reorder(TF, zscore), y = zscore))+
         geom_col(aes(fill = same), width = 0.7)) + coord_flip() +
  scale_y_continuous(name = "Z- Score", limits = c(0,2)) + # breaks, labels,  trans)+ 
  xlab("TF") 

p1+p2
