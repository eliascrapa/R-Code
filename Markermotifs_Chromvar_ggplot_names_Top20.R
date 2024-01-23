library(ggplot2)
library(xlsx)
library(ggpubr)
library(tidyverse)
require(gridExtra)
library("ggsci")
library("patchwork")
#library(stringr)
library(stringi)


library(grid)
library(tidyverse)
library(shadowtext)
library(extrafont)

setwd("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/")
data.dir <- './Markermotifs_Chromvar_norm_Top20'
dir.create(data.dir)
setwd(data.dir)
  

for (n in 1:25){
  
  tabledfTNT <- read.csv2(paste0("/Users/eliascrapa/Box/ArchR/All/Baserun/additional/CHROMVARdev/Chromvar_MotifsC",n,"_x_TNT.csv"), sep = ",", header = T)
  tabledfMHC <- read.csv2(paste0("/Users/eliascrapa/Box/ArchR/All/Baserun/additional/CHROMVARdev/Chromvar_MotifsC",n,"_x_MHC.csv"), sep = ",", header = T)
  
  #TNT
  tabledfTNT$name<- gsub("_.*","",tabledfTNT$name)
  colnames(tabledfTNT) <- c("TF_number", "group", "group_name", "seqnames", "idx","TF","Log2FC","FDR","zscore")
  tabledfTNT <- tabledfTNT %>% select("TF","zscore","FDR","group_name","idx") 
  tabledfTNT[,2:3] <- lapply(tabledfTNT[,2:3], as.numeric)
  tabledfTNT <- tabledfTNT[order(tabledfTNT$zscore),] 
  tabledfTNT <- tabledfTNT %>% distinct(TF, .keep_all = TRUE)
  tabledfTNT <- tabledfTNT %>% slice_max(tabledfTNT$'zscore',n = 30, with_ties =F)
  tabledfTNT <- tabledfTNT %>% drop_na() 
  tabledfTNT$TF <-abbreviate(tabledfTNT$TF, minlength = 8, dot = FALSE, strict = FALSE, method = c("left.kept"), named = TRUE)

  
  
  #MHC
  tabledfMHC$name<- gsub("_.*","",tabledfMHC$name)
  colnames(tabledfMHC) <- c("TF_number", "group", "group_name", "seqnames", "idx","TF","Log2FC","FDR","zscore")
  tabledfMHC <- tabledfMHC %>% select("TF","zscore","FDR","group_name","idx") 
  tabledfMHC[,2:3] <- lapply(tabledfMHC[,2:3], as.numeric)
  tabledfMHC <- tabledfMHC[order(tabledfMHC$zscore),] 
  tabledfMHC <- tabledfMHC %>% distinct(TF, .keep_all = TRUE)
  tabledfMHC <- tabledfMHC %>% slice_max(tabledfMHC$'zscore',n = 30, with_ties =F)
  tabledfMHC <- tabledfMHC %>% drop_na()
  tabledfMHC$TF <-abbreviate(tabledfMHC$TF, minlength = 8, dot = FALSE, strict = FALSE, method = c("left.kept"), named = TRUE)
  
  
  tabledfTNT$same<- ifelse(sapply(tabledfTNT$TF, `%in%`, (dplyr::intersect(tabledfTNT$TF,tabledfMHC$TF))),"matching", "different")
  tabledfMHC$same<- ifelse(sapply(tabledfMHC$TF, `%in%`, (dplyr::intersect(tabledfTNT$TF,tabledfMHC$TF))),"matching", "different")
  
   maxvalue <- plyr::round_any((max(union(tabledfTNT$zscore,tabledfMHC$zscore))),0.5,f = ceiling)
  
  print (maxvalue)
  
  tabledfTNT <- tabledfTNT %>% mutate_at(c('TF'), as.factor)
  #tabledfTNT$rank = seq(length(tabledfTNT$TF)) * 0.9
  #tabledfTNT$rank  <- as.numeric(rownames(tabledfTNT))
  #tabledfTNT$TF<- gsub(" ", "", tabledfTNT$TF, fixed = TRUE)
  #tabledfTNT$TF<- stri_replace_all_charclass(tabledfTNT$TF, "\\p{WHITE_SPACE}", "")
  
  tabledfTNT <- tabledfTNT[1:20,]
  tabledfMHC <- tabledfMHC[1:20,]
  
 p1 <- (ggplot(tabledfTNT, aes(x = reorder(TF, zscore), y = zscore))+
           geom_col(aes(fill = same))) + coord_flip() +
    ggtitle("TnT") +
   scale_y_continuous(name = "Z-Score", limits = c(0,maxvalue)) + # breaks, labels,  trans)+ 
    xlab("TF") + 
    #scale_fill_jco()  + 
   scale_fill_manual(values = c("different" ="#3c8cdc", "matching"="#e8c142"))+
    theme(legend.text=element_text(size=18),
          axis.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(), #remove x axis labels
          axis.ticks.y=element_blank(), #remove x axis ticks
          legend.title = element_blank(), 
          plot.title = element_text(hjust = 0.95, vjust = - 107, size = 20))#  #, plot.title = element_text(vjust = - 10)#+ #vlegend.position = "none") 
  #scale_fill_discrete(breaks = c("matching", "different"))
 
 
 if(nrow(tabledfTNT)>0) {
   p1 = p1 + geom_text(aes(label=TF, y= 0.05), colour = "black", nudge_x = +0.0, hjust=-0.0,  size=8)
 }
  
  p2 <- (ggplot(tabledfMHC, aes(x = reorder(TF, zscore), y = zscore))+
          geom_col(aes(fill = same))) + 
    coord_flip() +
    ggtitle("MyHC") +
    scale_y_continuous(name = "Z-Score", limits = c(0,maxvalue)) + # breaks, labels,  trans)+ 
    xlab("TF") + 
    #scale_fill_jco()  + 
    scale_fill_manual(values = c("different" ="#3c8cdc", "matching"="#e8c142"))+
    theme(legend.text=element_text(size=18),
          axis.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(), #remove x axis labels
          axis.ticks.y=element_blank(), #remove x axis ticks
          legend.title = element_blank(), 
          plot.title = element_text(hjust = 0.95, vjust = - 107, size = 20))#  #, plot.title = element_text(vjust = - 10)#+ #vlegend.position = "none") 
  #scale_fill_discrete(breaks = c("matching", "different"))
  
    
    if(nrow(tabledfMHC)>0) {
      p2 = p2 + geom_text(aes(label=TF, y= 0.05), colour = "black", nudge_x = +0.0, hjust=-0.0,  size=8)
    }
  
  design="11112223"
    
   # pp<- p1 + p2 + guide_area() + plot_layout(design=design, guides = "collect", widths = 50)# & theme(legend.position= c(.7, .8))#c("top"))
    pp<- p1 + p2 + plot_layout(guides = "collect", widths = c(50,50), heights = c(50,1)) & theme(legend.position= ("top"))
    pp
     ggsave(pp, device="pdf", file=paste0("MotifMarkers_C",n,"_TNT_MHCvsallCTR.pdf"), height = 8 , width = 12)
  
  # }
  
}

  #-------------------------------------------------------------------------


