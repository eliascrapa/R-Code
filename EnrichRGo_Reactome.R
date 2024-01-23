#https://ycl6.github.io/GO-Enrichment-Analysis-Demo/4_enrichR.html


library("clusterProfiler")
library("enrichplot")
library("org.Mm.eg.db")
library("ggplot2")
library(ggpubr)
library(xlsx)
library(tidyverse)
library("enrichR")


setwd("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/")

data.dir <- './EnrichRGo_Reactome/'
dir.create(data.dir)
setwd(data.dir)

data.dir <- './MyHC/'
dir.create(data.dir)
setwd(data.dir)
#TNT Export RNK File for Genesetenrichment
for (n in 20:20){

#source <- paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsMHC/RNA_ClusterMarker_TNT_vs_MHC_C",n,".xlsx")
#source <- paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_TNT_vs_CONTROLS_C",n,".xlsx")
source <- paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_MHC_vs_CONTROLS_C",n,".xlsx")
upgroup <-  paste0("Enriched in MyHC ")
downgroup <-  paste0("Enriched in CTR ")


TNTvsMHC <- read.xlsx(source, 1)

sapply(TNTvsMHC, class) 
#TNTvsMHC[,3:6] <- sapply(TNTvsMHC[,3:6],as.numeric)

TNT_logFC <- dplyr::filter(TNTvsMHC, p_val_adj < 0.05)
TNT_logFC_sel <- dplyr::select(TNT_logFC, , names, avg_log2FC,p_val_adj)  #select(TNT_logFC, -c(1,3,5,6))
#remove_rownames(TNT_logFC_sel)
TNT_logFC_sel <- dplyr::arrange(TNT_logFC_sel, desc(avg_log2FC))


if (!between(n, 6, 9)) {
  #TNT_logFC_sel <- subset(TNT_logFC_sel, names != (grep(TNT_logFC_sel$names,"^m")))
  TNT_logFC_sel <- TNT_logFC_sel[!grepl("Tnnt2",TNT_logFC_sel$names),]
  TNT_logFC_sel <- TNT_logFC_sel[!grepl("Hbb-bs",TNT_logFC_sel$names),]
  TNT_logFC_sel <- TNT_logFC_sel[!grepl("Myh6",TNT_logFC_sel$names),]
  TNT_logFC_sel <- TNT_logFC_sel[!grepl("^mt-",TNT_logFC_sel$names),]
  TNT_logFC_sel <- TNT_logFC_sel[!grepl("Ttn",TNT_logFC_sel$names),]
} else {
  
}

#names(TNT_logFC_sel) <- NULL
#write.table(TNT_logFC_sel, file = filem, sep = '\t',row.names = F, quote = FALSE) 


up.idx <- which(TNT_logFC_sel$p_val_adj < 0.05 & TNT_logFC_sel$avg_log2FC > 0) # FDR < 0.05 and logFC > 0
dn.idx <- which(TNT_logFC_sel$p_val_adj < 0.05 & TNT_logFC_sel$avg_log2FC < 0) # FDR < 0.05 and logFC < 0


length(up.idx)
length(dn.idx)

if (length(up.idx) >0) {
  
  if (length(dn.idx) >0) {
  
all.genes <- TNT_logFC_sel$names
up.genes <- TNT_logFC_sel[up.idx,]$names
dn.genes <- TNT_logFC_sel[dn.idx,]$names

up.genes
dn.genes

# Use fromType = "ENSEMBL" if your input identifier is Ensembl gene ID
up.genes.df = as.data.frame(up.genes)
dn.genes.df = as.data.frame(dn.genes)

colnames(dn.genes.df) <- "names"
colnames(up.genes.df) <- "names"

head(up.genes.df, 10)
head(dn.genes.df, 10)

#load databases
dbs <- listEnrichrDbs()
dbs <- dbs[order(dbs$libraryName),]

class(dbs)
dim(dbs)
head(dbs)

#Check which databases are available
dbs$libraryName
dbs[grep("2021",dbs$libraryName),]$libraryName 

dbs_go <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")
dbs_pw <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse")#, "BioPlanet_2019")
dbs_dd <- c("Reactome_2016")

#Enrichment of GO Terms

upEnriched_go <- enrichr(genes = up.genes.df$names, databases = dbs_go)
dnEnriched_go <- enrichr(genes = dn.genes.df$names, databases = dbs_go)

#class(upEnriched_go)
#names(upEnriched_go)

# View top 5 terms in the first element of the list
#head(upEnriched_go[[1]], 15)


#Pathways
upEnriched_pw <- enrichr(genes = up.genes.df$names, databases = dbs_pw)
dnEnriched_pw <- enrichr(genes = dn.genes.df$names, databases = dbs_pw)

#head(upEnriched_pw[[3]], 10)

#Diseases/Drugs analysis

upEnriched_dd <- enrichr(genes = up.genes.df$names, databases = dbs_dd)
dnEnriched_dd <- enrichr(genes = dn.genes.df$names, databases = dbs_dd)

#head(dnEnriched_dd[[2]], 10)




#upEnriched_go[[1]][order(upEnriched_go[["GO_Molecular_Function_2021"]][["Adjusted.P.value"]]),] %>% filter(upEnriched_go[["GO_Molecular_Function_2021"]][["Adjusted.P.value"]] < "0.25")


plotlist <- list()


#GOdatabase <- paste0(gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_go)[l]))))


#GO-Terms
GOdatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_go)[1]))))
goup  <- plotEnrich(upEnriched_go[[1]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup) 
godown <- plotEnrich(dnEnriched_go[[1]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup) 
write.xlsx((upEnriched_go[[1]][order(upEnriched_go[[1]][["Adjusted.P.value"]]),] %>% filter(upEnriched_go[[1]][["Adjusted.P.value"]] < "0.25")), file = paste0(GOdatabase," - ",upgroup,"Cluster",n," FDR under 0.25.xlsx"))
write.xlsx((dnEnriched_go[[1]][order(dnEnriched_go[[1]][["Adjusted.P.value"]]),] %>% filter(dnEnriched_go[[1]][["Adjusted.P.value"]] < "0.25")), file = paste0(GOdatabase," - ",downgroup,"Cluster",n," FDR under 0.25.xlsx"))
plot <- ggarrange(goup,godown, labels =c("A","B"))
p1 <-annotate_figure(plot, top = text_grob(GOdatabase, color = "red", face = "bold", size = 17)) 

GOdatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_go)[2]))))
goup  <- plotEnrich(upEnriched_go[[2]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup) 
godown <- plotEnrich(dnEnriched_go[[2]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup) 
write.xlsx((upEnriched_go[[2]][order(upEnriched_go[[2]][["Adjusted.P.value"]]),] %>% filter(upEnriched_go[[2]][["Adjusted.P.value"]] < "0.25")), file = paste0(GOdatabase," - ",upgroup,"Cluster",n," FDR under 0.25.xlsx"))
write.xlsx((dnEnriched_go[[2]][order(dnEnriched_go[[2]][["Adjusted.P.value"]]),] %>% filter(dnEnriched_go[[2]][["Adjusted.P.value"]] < "0.25")), file = paste0(GOdatabase," - ",downgroup,"Cluster",n," FDR under 0.25.xlsx"))


plot <- ggarrange(goup,godown, labels =c("A","B"))

plot <- ggarrange(goup,godown, labels =c("A","B"))
p2<- annotate_figure(plot, top = text_grob(GOdatabase, color = "red", face = "bold", size = 17)) 

GOdatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_go)[3]))))
goup  <- plotEnrich(upEnriched_go[[3]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup) 
godown <- plotEnrich(dnEnriched_go[[3]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup) 
write.xlsx((upEnriched_go[[3]][order(upEnriched_go[[3]][["Adjusted.P.value"]]),] %>% filter(upEnriched_go[[3]][["Adjusted.P.value"]] < "0.25")), file = paste0(GOdatabase," - ",upgroup,"Cluster",n," FDR under 0.25.xlsx"))
write.xlsx((dnEnriched_go[[3]][order(dnEnriched_go[[3]][["Adjusted.P.value"]]),] %>% filter(dnEnriched_go[[3]][["Adjusted.P.value"]] < "0.25")), file = paste0(GOdatabase," - ",downgroup,"Cluster",n," FDR under 0.25.xlsx"))

plot <- ggarrange(goup,godown, labels =c("A","B"))
p3<- annotate_figure(plot, top = text_grob(GOdatabase, color = "red", face = "bold", size = 17)) 




#Pathwayanalysis

Pathwaydatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_pw)[1]))))
pwup  <- plotEnrich(upEnriched_pw[[1]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup) 
pwdown <- plotEnrich(dnEnriched_pw[[1]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup) 
write.xlsx((upEnriched_pw[[1]][order(upEnriched_pw[[1]][["Adjusted.P.value"]]),] %>% filter(upEnriched_pw[[1]][["Adjusted.P.value"]] < "0.25")), file = paste0(Pathwaydatabase," - ",upgroup,"Cluster",n," FDR under 0.25.xlsx"))
write.xlsx((dnEnriched_pw[[1]][order(dnEnriched_pw[[1]][["Adjusted.P.value"]]),] %>% filter(dnEnriched_pw[[1]][["Adjusted.P.value"]] < "0.25")), file = paste0(Pathwaydatabase," - ",downgroup,"Cluster",n," FDR under 0.25.xlsx"))

plot <- ggarrange(pwup,pwdown, labels =c("A","B"))
p4 <-annotate_figure(plot, top = text_grob(Pathwaydatabase, color = "red", face = "bold", size = 17)) 

Pathwaydatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_pw)[2]))))
pwup  <- plotEnrich(upEnriched_pw[[2]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup) 
pwdown <- plotEnrich(dnEnriched_pw[[2]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup) 
write.xlsx((upEnriched_pw[[2]][order(upEnriched_pw[[2]][["Adjusted.P.value"]]),] %>% filter(upEnriched_pw[[2]][["Adjusted.P.value"]] < "0.25")), file = paste0(Pathwaydatabase," - ",upgroup,"Cluster",n," FDR under 0.25.xlsx"))
write.xlsx((dnEnriched_pw[[2]][order(dnEnriched_pw[[2]][["Adjusted.P.value"]]),] %>% filter(dnEnriched_pw[[2]][["Adjusted.P.value"]] < "0.25")), file = paste0(Pathwaydatabase," - ",downgroup,"Cluster",n," FDR under 0.25.xlsx"))

plot <- ggarrange(pwup,pwdown, labels =c("A","B"))
p5 <-annotate_figure(plot, top = text_grob(Pathwaydatabase, color = "red", face = "bold", size = 17)) 

Pathwaydatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_dd)[1]))))
ddup  <- plotEnrich(upEnriched_dd[[1]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup) 
dddown <- plotEnrich(dnEnriched_dd[[1]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup) 
write.xlsx((upEnriched_dd[[1]][order(upEnriched_dd[[1]][["Adjusted.P.value"]]),] %>% filter(upEnriched_dd[[1]][["Adjusted.P.value"]] < "0.25")), file = paste0(Pathwaydatabase," - ",upgroup,"Cluster",n," FDR under 0.25.xlsx"))
write.xlsx((dnEnriched_dd[[1]][order(dnEnriched_dd[[1]][["Adjusted.P.value"]]),] %>% filter(dnEnriched_dd[[1]][["Adjusted.P.value"]] < "0.25")), file = paste0(Pathwaydatabase," - ",downgroup,"Cluster",n," FDR under 0.25.xlsx"))
 
plot <- ggarrange(ddup,dddown, labels =c("A","B"))
p6 <-annotate_figure(plot, top = text_grob(Pathwaydatabase, color = "red", face = "bold", size = 17)) 


plotlist <- list(p1,p2,p3,p4,p5,p6)

ggexport(plotlist = plotlist, filename = paste0("Cluster ",n," - GO and Pathway Analysis - Enriched in ",upgroup,"Cluster",n," FDR under 0.25.pdf"), width=20, height=11, res=300)



#d1 <- plotEnrich(upEnriched_dd[[1]], showTerms = 10, numChar = 30, y = "Count", orderBy = "Combined.Score", title = upgroup)
#d2 <- plotEnrich(dnEnriched_dd[[1]], showTerms = 10, numChar = 30, y = "Count", orderBy = "Combined.Score", title = downgroup)
#d1+d2


} else {
  

  
}
  
  
} else {
  
}

}










