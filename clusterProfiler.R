library("clusterProfiler")
library("enrichplot")
library("org.Mm.eg.db")
library("ggplot2")
library(xlsx)
library(tidyverse)


setwd("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/")

data.dir <- './clusterProfiler/'
dir.create(data.dir)
setwd(data.dir)
#TNT Export RNK File for Genesetenrichment
for (n in 1:25){

source <- paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsMHC/RNA_ClusterMarker_TNT_vs_MHC_C",n,".xlsx")
source <- paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_TNT_vs_CONTROLS_C",n,".xlsx")
filem <-  paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/RNK/TNTvsCTR_C",n,".rnk")  


TNTvsMHC <- read.xlsx(source, 1)

sapply(TNTvsMHC, class) 
#TNTvsMHC[,3:6] <- sapply(TNTvsMHC[,3:6],as.numeric)

TNT_logFC <- dplyr::filter(TNTvsMHC, p_val_adj < 0.05)
TNT_logFC_sel <-dplyr::select(TNT_logFC, , names, avg_log2FC,p_val_adj)  #select(TNT_logFC, -c(1,3,5,6))
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

all.genes <- TNT_logFC_sel$names
up.genes <- TNT_logFC_sel[up.idx,]$names
dn.genes <- TNT_logFC_sel[dn.idx,]$names

up.genes
dn.genes


#BP: Biological Process
#CC: Cellular Component
#MF: Molecular Function

ontology <- "BP"

outTitle <- paste0("clusterProfiler_GO-", ontology, "_ORA_simplify")
outTitle

# Use fromType = "ENSEMBL" if your input identifier is Ensembl gene ID
all.genes.df = bitr(all.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

head(all.genes.df, 10)

up.genes.df = bitr(up.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

head(up.genes.df, 10)

dn.genes.df = bitr(dn.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

head(dn.genes.df, 10)


upEGO = enrichGO(gene = up.genes.df$ENTREZID, universe = all.genes.df$ENTREZID, 
                 OrgDb = "org.Mm.eg.db", ont = ontology, pvalueCutoff = 0.8, 
                 pAdjustMethod = "BH", qvalueCutoff = 0.8, readable = TRUE)
upEGO

dnEGO = enrichGO(gene = dn.genes.df$ENTREZID, universe = all.genes.df$ENTREZID,
                 OrgDb = "org.Mm.eg.db", ont = ontology, pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
dnEGO


upSimGO = simplify(upEGO, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", 
                   semData = NULL)
nrow(upSimGO)

dnSimGO = simplify(dnEGO, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", 
                   semData = NULL)
nrow(dnSimGO)

png(paste0(outTitle, "_up.png"), width = 9, height = 6, units = "in", res = 300)
barplot(upSimGO, showCategory = 20) + 
  ggtitle(paste0("GO-", ontology," ORA of up-regulated genes")) + 
  xlab("Enriched terms") + ylab("Count")
invisible(dev.off())

}






















#MHC Export RNK File for Genesetenrichment
for (n in 1:25){
  
  source <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_MHC_vs_CONTROLS_C",n,".xlsx")
  filem <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/RNK/MHCvsCTR_C",n,".rnk")  
  
  
  TNTvsCTR_C <- read.xlsx(source, 1)
  
  sapply(TNTvsCTR_C, class) 
  TNTvsCTR_C[,3:6] <- sapply(TNTvsCTR_C[,3:6],as.numeric)
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  
  TNT_logFC <- dplyr::filter(TNTvsCTR_C, p_val_adj < 0.05)
  TNT_logFC_sel <-dplyr::select(TNT_logFC, , names, avg_log2FC)  #select(TNT_logFC, -c(1,3,5,6))
  #remove_rownames(TNT_logFC_sel)
  TNT_logFC_sel <- dplyr::arrange(TNT_logFC_sel, desc(avg_log2FC))
  
  
  if (!between(n, 6, 9)) {
    TNT_logFC_sel <- subset(TNT_logFC_sel, names != "Tnnt2")
    TNT_logFC_sel <- subset(TNT_logFC_sel, names != "Hbb-bs")
    TNT_logFC_sel <- subset(TNT_logFC_sel, names != "Myh6")
    #TNT_logFC_sel <- subset(TNT_logFC_sel, names != (grep(TNT_logFC_sel$names,"^m")))
    TNT_logFC_sel <- TNT_logFC_sel[!grepl("^mt-",TNT_logFC_sel$names),]
    TNT_logFC_sel <- TNT_logFC_sel[!grepl("Ttn",TNT_logFC_sel$names),]
  } else {
    
  }
  
  names(TNT_logFC_sel) <- NULL
  write.table(TNT_logFC_sel, file = filem, sep = '\t',row.names = F, quote = FALSE) 
}


  
}
# Gensetenrichment all
my_sequence = c(1,2,3,4,6,8,9,12,13,14,15,16,17,19,20,21,22,23,25)
my_sequence = c(1,2,3)

database <- c('geneontology_Biological_Process_noRedundant','geneontology_Cellular_Component_noRedundant','geneontology_Molecular_Function_noRedundant','pathway_KEGG','pathway_Reactome','pathway_Wikipathway','network_Transcription_Factor_target','network_miRNA_target','phenotype_Mammalian_Phenotype_Ontology','community-contributed_MuscleGeneSets_Duddy_2017')


database <- readGmt("/Users/eliascrapa/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/enrichmentdatabase/FROMKEGG/geneset.gmt", cache = NULL)

genotype  <- c('TNT','MHC')
for (k in genotype){

for (i in database){
  
for (n in 1:1 ){
  #1:25
  
  source <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/RNK/",k,"vsCTR_C",n,".rnk")  
  project <- paste0(i,"_",k,"vsCTR_Cluster_",n)
  
  
  
  outputDirectory <- getwd()
  
  skip_to_next <- FALSE
  
 

  
  tryCatch(gsea_res <- WebGestaltR::WebGestaltR(enrichMethod = "GSEA",
                                       organism="mmusculus",
                                       enrichDatabaseType = "genesymbol",
                                       interestGeneType = "genesymbol",
                                       interestGeneFile= source,
                                       sigMethod="top",
                                       topThr=10,
                                       minNum = 3,
                                       maxNum = 5000,
                                       enrichDatabase=i,
                                       outputDirectory=outputDirectory,
                                       projectName = project),
                                       
                                       error = function(e) { skip_to_next <<- TRUE})
  
           if(skip_to_next) { next }  
  
}

}
  

}
prepareGseaInput#name


#1                  geneontology_Biological_Process
#2      geneontology_Biological_Process_noRedundant
#3                  geneontology_Cellular_Component
#4      geneontology_Cellular_Component_noRedundant
#5                  geneontology_Molecular_Function
#6      geneontology_Molecular_Function_noRedundant
#7                                     pathway_KEGG
#8                                  pathway_Panther
#9                                 pathway_Reactome
#10                             pathway_Wikipathway
#11                                   network_CORUM
#12                                  network_CORUMA
#13                      network_Kinase_phosphosite
#14                           network_Kinase_target
#15                             network_PPI_BIOGRID
#16                                network_PTMsigDB
#17             network_Transcription_Factor_target
#18                            network_miRNA_target
#19          phenotype_Mammalian_Phenotype_Ontology
#20             chromosomalLocation_CytogeneticBand
#21         community-contributed_5htGeneSets_Conte
#22 community-contributed_MuscleGeneSets_Duddy_2017

GMTall <- c()
for (q in 1:63) {
  GMT <- as.data.frame(read.table(paste0("/Users/eliascrapa/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/enrichmentdatabase/geneset-",q,".gmt"), sep="\t", header =F))
print(GMT)
  GMTall <- rbind.fill(GMTall,GMT)
}
write.table(GMTall, file="GMTAll.gmt", sep="\t",
         row.names = F, col.names = F, na = "NATT", quote = FALSE)

genotype  <- c('TNT','MHC')
for (k in genotype){
    
    for (n in 1:1 ){
      #1:25
      
      source <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/RNK/",k,"vsCTR_C",n,".rnk")  
      project <- paste0(i,"_",k,"vsCTR_Cluster_",n)
      
      
      
      outputDirectory <- getwd()
      
      skip_to_next <- FALSE
  # rank<-  as.data.frame(read_table(source,col_names = c('gene','score')))
        
  #input <-    prepareInputMatrixGsea(rank, GMTall)
      
      
      tryCatch(gsea_res <- WebGestaltR::WebGestaltR(enrichMethod = "GSEA",
                                                    organism="hsapiens",
                                                    enrichDatabaseType = "genename",
                                                    interestGeneType = "genesymbol",
                                                    interestGeneFile= source,
                                                    sigMethod="top",
                                                    topThr=10,
                                                    minNum = 3,
                                                    maxNum = 5000,
                                                    enrichDatabaseFile="GMTall.gmt",
                                                    outputDirectory=outputDirectory,
                                                    projectName = project,
                                                    gseaPlotFormat = c("png", "svg")
                                                    ),
                        
               
               error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) { next }  
      
    }
    
  }
  
  






#TNT Gensetenrichment
my_sequence = c(1,2,3,4,6,8,9,12,13,14,15,16,17,19,20,21,22,23,25)
my_sequence = c(16,17,19,20,21,22,23,25)
for (n in my_sequence){
  
  
  source <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsCTR_C",n,".rnk")  
  project <- paste0("phenotype_Mammalian_Phenotype_Ontology_TNTvsCTR_Cluster_",n)
  
  
  rankFile <- system.file("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsCTR_C9.rnk", package="WebGestaltR")
  outputDirectory <- getwd()
  
  gsea_res <- WebGestaltR::WebGestaltR(enrichMethod = "GSEA",
                                       organism="mmusculus",
                                       enrichDatabaseType = "genesymbol",
                                       interestGeneType = "genesymbol",
                                       interestGeneFile= source,
                                       sigMethod="top",
                                       topThr=10,
                                       minNum = 3,
                                       maxNum = 5000,
                                       enrichDatabase="phenotype_Mammalian_Phenotype_Ontology",
                                       outputDirectory=outputDirectory,
                                       projectName = project
                                       
  )
  
}


