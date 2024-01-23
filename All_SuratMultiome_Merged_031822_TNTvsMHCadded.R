#BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))
#BiocManager::install("scDblFinder")
#install.packages('stringr')
#BiocManager::install("MAST")
#devtools::install_github('xzhoulab/iDEA')
#install.packages('GOplot')
#if (!requireNamespace('BiocManager', quietly = TRUE))
#install.packages('BiocManager')

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)
library(hdf5r)
library(patchwork)
library(stringr)
library(rtracklayer)
SeqinfoForUCSCGenome("mm10")
library(scDblFinder)
library(SingleCellExperiment)
library(harmony)
library("xlsx")
library(MAST)
library(GOplot)
library(EnhancedVolcano)
library(magrittr)
library(dittoSeq)
library(SoupX)
library(tidyverse)
library(DropletUtils)

#Set Analysisfolder::
analysisfolder <- "~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All"
#setwd("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All")
setwd(analysisfolder)


# We’ll create a Seurat object based on the gene expression data, and then add in 
# theATAC-seq data as a second assay. You can explore the Signac getting started 
# vignette for more information on the creation and processing of a ChromatinAssay object.



#install.packages("SoupX")


#####################DATAAdding loop

# the 10x hdf5 file contains both data types. 


#DATAADDING with SoupX Cleanup and Cleanup
{

setwd("/Users/eliascrapa/DATA")

sampleinputlist <- list()
#sampleinputlist <- as.list(c("CTR-TNT","MHC","CTR-MHC"))
sampleinputlist <- as.list(c("CTR-TNT","CTR-MHC"))
#sampleinputlist <- as.list("TNT")
for (k in sampleinputlist) {
  
  #Adding Counts and ATAC Data TNT
  data.dir <- paste0(analysisfolder,'/',k)
  dir.create(data.dir)
  
  #Soup-X Contamination Filtering  
  # Load data and estimate soup profile
  # extract RNA 
  sc = load10X(paste0("./",k,"/"))   
  
  
  
  Samplecluster <- read.csv(paste0("./",k,"/analysis/clustering/gex/graphclust/clusters.csv"))
  sc = setClusters(sc, setNames(Samplecluster$Cluster,Samplecluster$Barcode)) #setNames(clusters, barcodes)
  
  UMAPGEX <- read.csv(paste0("./",k,"/analysis/dimensionality_reduction/gex/umap_projection.csv"))
  UMAPGEX$Cluster <- Samplecluster$Cluster
  UMAPGEX$Annotation <- Samplecluster$Cluster
  sc = setDR(sc, UMAPGEX[,c("UMAP.1","UMAP.2")])
  
  
  #visual sanity check
  library(ggplot2)
  dd = UMAPGEX#[colnames(sc$toc), ]
  mids = aggregate(cbind(UMAP.1,UMAP.2) ~ Annotation, data = dd, FUN = mean)
  gg = ggplot(dd, aes(UMAP.1, UMAP.2)) + geom_point(aes(colour = Annotation), size = 0.2) + ggtitle("TNT") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
  plot(gg)
  
  gg = ggplot(dd, aes(UMAP.1, UMAP.2))+ geom_point(size = 0.2)
  plot(gg)
  
  dd$Tnnt2 = sc$toc["Tnnt2", ]
  gg = ggplot(dd, aes(UMAP.1, UMAP.2)) + geom_point(aes(colour = Tnnt2 > 0))
  plot(gg)
  dd$'Hbb-bs' = sc$toc["Hbb-bs", ]
  gg = ggplot(dd, aes(UMAP.1, UMAP.2)) + geom_point(aes(colour = 'Hbb-bs' > 0))
  plot(gg)
  
  gg = plotMarkerMap(sc, "Tnnt2", pointSize = 0.5)
  plot(gg)
  gg = plotMarkerMap(sc, "Hbb-bs", pointSize = 0.5)
  plot(gg)
  gg = plotMarkerDistribution(sc)
  plot(gg)
  
  print(paste0(k))
  png(paste0(analysisfolder,"/",k,"/",k,"Soupestimation_",k,".png"))
  sc = autoEstCont(sc)
  dev.off()
  
  sc = setContaminationFraction(sc, 0.2)
  #Genes with highest expression in background.
  BackgroundgenesSample <- head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 100)
  write.csv(BackgroundgenesSample, file = paste0(analysisfolder,"/",k,"/",k,"Backgroundgenes_",k,".csv"))
  
  
  backgroundlist <- list()
  backgroundlist <- as.list(head(rownames(BackgroundgenesSample), n = 25))
  
  
  #adjust Counts based on SoupX, save Counts to SparseMatrix and H5-file
  adj.rna_counts <- adjustCounts(sc, roundToInt = T)
  
  
  DropletUtils:::write10xCounts(paste0("./",k,"/soupX_TNT_filt"), adj.rna_counts,
                                overwrite = TRUE)
  
  DropletUtils:::write10xCounts(paste0("./",k,"/soupX_TNT_filt.h5"), adj.rna_counts, type = ("HDF5"),
                                #barcodes = colnames(x),
                                #gene.id = rownames(x),
                                #gene.symbol = gene.id,
                                #gene.type = "Gene Expression",
                                overwrite = TRUE,
                                genome = "MM10",
                                version = ("3"),
                                chemistry = "10x Multiome v1",
                                original.gem.groups = '4d4',
                                library.ids = "TNT"
  )
  
  
  for (l in backgroundlist){
    png(paste0(analysisfolder,"/",k,"/",l,"_ChangeExpressionSoupX.png"), width = 1000, height = 800, units = "px")
    print(l)
    gg = plotChangeMap(sc, adj.rna_counts, paste(l))
    plot(gg)
    dev.off()
    
    
  }
  
  
  inputdata.10x <- Read10X_h5(paste0("./",k,"/filtered_feature_bc_matrix.h5"))
  # extract RNA and ATAC data
  #rna_counts <- inputdata.10x$`Gene Expression`
  atac_counts <- inputdata.10x$Peaks
  rm(inputdata.10x)
  
  # Create Seurat object
  assign(paste0(k),CreateSeuratObject(counts = adj.rna_counts))
  
  
  
  
  
  # Now add in the ATAC-seq data
  # we'll only use peaks in standard chromosomes
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
  seqlevels(annotations) <- ucsc.levels
  genome(annotations) <- "mm10"
  
  frag.file <- paste0("./",k,"/atac_fragments.tsv.gz")
  assign(paste0(k,"chrom_assay"), CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = frag.file,
    min.cells = 10,
    annotation = annotations
  ))
  
  
  rm(backgroundlist)
  rm(dd)
  rm(gg)
  rm(l)
  rm(Samplecluster)
  rm(sc)
  rm(mids)
  rm(UMAPGEX)
  rm(atac_counts)
  rm(adj.rna_counts)
  rm(grange.counts)
  rm(grange.use)
  rm(BackgroundgenesSample)
  rm(frag.file)
  rm(ucsc.levels)
  rm(data.dir)
  rm(annotations)
  
}                   
rm(k)                

CTR.TNT <- `CTR-TNT`
CTR.MHC <- `CTR-MHC`
CTR.TNTchrom_assay <- `CTR-TNTchrom_assay`
CTR.MHCchrom_assay <- `CTR-MHCchrom_assay`

TNT[["percent.mt"]] <- PercentageFeatureSet(TNT, pattern = "^mt-")
TNT[["ATAC"]] <- TNTchrom_assay

CTR.TNT[["percent.mt"]] <- PercentageFeatureSet(CTR.TNT, pattern = "^mt-")
CTR.TNT[["ATAC"]] <- CTR.TNTchrom_assay

MHC[["percent.mt"]] <- PercentageFeatureSet(MHC, pattern = "^mt-")
MHC[["ATAC"]] <- MHCchrom_assay

CTR.MHC[["percent.mt"]] <- PercentageFeatureSet(CTR.MHC, pattern = "^mt-")
CTR.MHC[["ATAC"]] <- CTR.MHCchrom_assay


rm(`TNTchrom_assay`)

rm(`CTR-TNTchrom_assay`)

rm(`CTR.TNTchrom_assay`)

rm(`MHCchrom_assay`)

rm(`CTR-MHCchrom_assay`)
rm(`CTR.MHCchrom_assay`)

rm(CTR-MHC)
rm(CTR-TNT)


ctr.mhc <- CTR.MHC
rm(CTR.MHC)
ctr.tnt <- CTR.TNT
rm(CTR.TNT)
mhc <- MHC
rm(MHC)
tnt <- TNT
rm(TNT)

}



#SAVE Projectevniroment
save.image(file='/Volumes/415-5899036/SingleCellSequencing/Seurat/AllDataAdded.RData')


#RNA and ATAC Pre-Analysis
{
  
  #CTR.MHC
  {
# RNA analysis
DefaultAssay(ctr.mhc) <- "RNA"
ctr.mhc <- SCTransform(ctr.mhc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(ctr.mhc) <- "ATAC"
ctr.mhc <- RunTFIDF(ctr.mhc)
ctr.mhc <- FindTopFeatures(ctr.mhc, min.cutoff = 'q0')
ctr.mhc <- RunSVD(ctr.mhc)
ctr.mhc <- RunUMAP(ctr.mhc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#Integrated Analyis
ctr.mhc <- FindMultiModalNeighbors(ctr.mhc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
ctr.mhc <- RunUMAP(ctr.mhc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
ctr.mhc <- FindClusters(ctr.mhc, graph.name = "wsnn", algorithm = 3, verbose = TRUE, resolution=1.2) #was resolution=1.2
DefaultAssay(ctr.mhc) <- "SCT"
ctr.mhc <- FindVariableFeatures(ctr.mhc, selection.method = "vst", nfeatures = 2000)
ctr.mhc <- RunPCA(ctr.mhc, features = VariableFeatures(object = ctr.mhc), verbose = TRUE)

ctr.mhc



  }
  
  #MHC
  {
    # RNA analysis
    DefaultAssay(mhc) <- "RNA"
    mhc <- SCTransform(mhc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    
    
    # ATAC analysis
    # We exclude the first dimension as this is typically correlated with sequencing depth
    DefaultAssay(mhc) <- "ATAC"
    mhc <- RunTFIDF(mhc)
    mhc <- FindTopFeatures(mhc, min.cutoff = 'q0')
    mhc <- RunSVD(mhc)
    mhc <- RunUMAP(mhc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    
    #Integrated Analyis
    mhc <- FindMultiModalNeighbors(mhc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
    mhc <- RunUMAP(mhc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    mhc <- FindClusters(mhc, graph.name = "wsnn", algorithm = 3, verbose = TRUE, resolution=1.2) #was resolution=1.2
    DefaultAssay(mhc) <- "SCT"
    mhc <- FindVariableFeatures(mhc, selection.method = "vst", nfeatures = 2000)
    mhc <- RunPCA(mhc, features = VariableFeatures(object = mhc), verbose = TRUE)
    
    mhc

  }
  
  #CTR.TNT
  {
    # RNA analysis
    DefaultAssay(ctr.tnt) <- "RNA"
    ctr.tnt <- SCTransform(ctr.tnt, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    
    
    # ATAC analysis
    # We exclude the first dimension as this is typically correlated with sequencing depth
    DefaultAssay(ctr.tnt) <- "ATAC"
    ctr.tnt <- RunTFIDF(ctr.tnt)
    ctr.tnt <- FindTopFeatures(ctr.tnt, min.cutoff = 'q0')
    ctr.tnt <- RunSVD(ctr.tnt)
    ctr.tnt <- RunUMAP(ctr.tnt, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    
    #Integrated Analyis
    ctr.tnt <- FindMultiModalNeighbors(ctr.tnt, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
    ctr.tnt <- RunUMAP(ctr.tnt, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    ctr.tnt <- FindClusters(ctr.tnt, graph.name = "wsnn", algorithm = 3, verbose = TRUE, resolution=1.2) #was resolution=1.2
    DefaultAssay(ctr.tnt) <- "SCT"
    ctr.tnt <- FindVariableFeatures(ctr.tnt, selection.method = "vst", nfeatures = 2000)
    ctr.tnt <- RunPCA(ctr.tnt, features = VariableFeatures(object = ctr.tnt), verbose = TRUE)
    
    ctr.tnt
  
  }
  
  #TNT
  {
    # RNA analysis
    DefaultAssay(tnt) <- "RNA"
    tnt <- SCTransform(tnt, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    
    
    # ATAC analysis
    # We exclude the first dimension as this is typically correlated with sequencing depth
    DefaultAssay(tnt) <- "ATAC"
    tnt <- RunTFIDF(tnt)
    tnt <- FindTopFeatures(tnt, min.cutoff = 'q0')
    tnt <- RunSVD(tnt)
    tnt <- RunUMAP(tnt, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    
    #Integrated Analyis
    tnt <- FindMultiModalNeighbors(tnt, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
    tnt <- RunUMAP(tnt, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    tnt <- FindClusters(tnt, graph.name = "wsnn", algorithm = 3, verbose = TRUE, resolution=1.2) #was resolution=1.2
    DefaultAssay(tnt) <- "SCT"
    tnt <- FindVariableFeatures(tnt, selection.method = "vst", nfeatures = 2000)
    tnt <- RunPCA(tnt, features = VariableFeatures(object = tnt), verbose = TRUE)
    
    tnt
    
    
  }
  
  rm(mca_result)
  rm(rawcounts)
  rm(Seurat_Object_Diet)
  rm(plot1)
  rm(plot2)
  rm(ctr.mhc.sce)
}

DimPlot(tnt, reduction = "wnn.umap")  
DimPlot(ctr.tnt, reduction = "wnn.umap") 
DimPlot(mhc, reduction = "wnn.umap") 
DimPlot(ctr.mhc, reduction = "wnn.umap") 



#SAVE Projectevniroment
save.image(file='/Volumes/415-5899036/SingleCellSequencing/Seurat/Premerge.RData')
#Merging-Workflow
#Mergeprocedure #Idents
{
  #Adding Group to orig.ident
  levels(tnt@meta.data[["orig.ident"]]) <- rep("TNT",12764)
  levels(ctr.tnt@meta.data[["orig.ident"]]) <- rep("CTR-TNT",12140)
  levels(mhc@meta.data[["orig.ident"]]) <- rep("MHC",15827)
  levels(ctr.mhc@meta.data[["orig.ident"]]) <- rep("CTR-MHC",9448)
  
all <- merge(ctr.mhc, y = c(mhc, ctr.tnt, tnt), add.cell.ids = c("CTR-MHC", "MHC", "CTR-TNT", "TNT"), project = "5wHCM", merge.data = TRUE)
all
#Transfer Identy for active.ident to orig.ident
#all@meta.data[["orig.ident"]] <- all@active.ident
#all@meta.data[["orig.ident"]] <- gsub("Object","",all@meta.data[["orig.ident"]])
#cellsseurat <- as.data.frame(all@active.ident)

#SAVE Projectevniroment
save.image(file='/Volumes/415-5899036/SingleCellSequencing/Seurat/DataMerged.RData')

}

setwd(analysisfolder)

#Exporting Embeddings from ArchrR Object - ArchR Object has to be loaded
{
  #df_harmony <- proj@embeddings@listData[["UMAP_Harmony"]]@listData[["df"]]
  #rownames(df_harmony) <- gsub("#","_",rownames(df_harmony))
  #matrix_harmony <- as.matrix(df_harmony)
  
  #####################DATA has to be exstracted from single ArchR proj file -----
  #all the data to a dataframe + Exporting the Embedding and changing the Cellname to Match the seuratnaming
  df<- proj@embeddings@listData[["UMAP_Harmony"]]@listData[["df"]]
  df$ArchR_Clusters <- proj@cellColData@listData[["Clusters"]]
  #df$ClusterBySample <- proj@cellColData@listData[["ClusterBySample"]]
  #df <- df[ -c(1,2) ]
  rownames(df) <- gsub("#","_",rownames(df))
  #Add exported Metadata to Seurat
  all <- AddMetaData(all, metadata = df)
 
  
  #Adding Harmony-UMAP as a reduction to Seurat
  harmonyumpa <- as.data.frame(all@meta.data[["Harmony.UMAP_Dimension_1"]])
  harmonyumpa$"Harmony.UMAP_Dimension_2" <- (all@meta.data[["Harmony.UMAP_Dimension_2"]])
  names(harmonyumpa)[names(harmonyumpa) == 'all@meta.data[["Harmony.UMAP_Dimension_1"]]'] <- "Harmony.UMAP_Dimension_1"
  harmonyumpa <- as.matrix(harmonyumpa)
  all@reductions$harmony_archr <- CreateDimReducObject(embeddings = harmonyumpa, key = "archr_", assay="RNA")
  rownames(all[["harmony_archr"]]@cell.embeddings) <- Cells(all)
  
  rm(harmonyumpa)
  rm(df)
  
  #Change Idents to ArchR_Clusters and Sort them 
  Idents(object = all) <- "ArchR_Clusters"
  my_levels <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25")
  Idents(all) <- factor(Idents(all), levels= my_levels)
  
  
  #Subset Cells to remove Cells that are in Archr and not in Seurat
  all <- subset(
    x = all,
    subset = Harmony.UMAP_Dimension_1 != "NA"
  ) 
 
  all#all NA have to be set off!!
  
  all@meta.data$ClusterbySample <- paste0(all@meta.data$ArchR_Clusters, "_x_", all@meta.data$orig.ident)
  Idents(object = all) <- all@meta.data$ClusterbySample
  all
  
  #SAVE Projectevniroment
  save.image(file='/Volumes/415-5899036/SingleCellSequencing/Seurat/DataMerged_ClusterMetadata_ad.RData')

  
  all <- subset(
    x = all,
    subset = nCount_ATAC < 1.2e5 &
      nCount_ATAC > 1e3 &
      nCount_RNA < 30000 &
      nCount_RNA > 300 &
      percent.mt < 30
  )
  
  all
  #row.names(df) <- lapply(row.names(df), function(x) paste(x, sep="_", "1"))
  #df <- cbind(rownames(df), data.frame(df, row.names=NULL))
  #df$`rownames(df)` <- gsub(".*#","",df$`rownames(df)`)
  
  all@meta.data[["orig.ident"]] <- gsub("CTR-MHC","MHC-CTR",all@meta.data[["orig.ident"]])
  all@meta.data[["orig.ident"]] <- gsub("CTR-TNT","TNT-CTR",all@meta.data[["orig.ident"]])
}


#Data-Cleanup
{
#save.image(file='/Volumes/415-5899036/SingleCellSequencing/Seurat/AllDataMerged.RData')

rm (ctr.mhc)
rm (ctr.tnt) 
rm (mhc)
rm (tnt)
rm (cellsseurat)
save.image(file='/Volumes/415-5899036/SingleCellSequencing/Seurat/reductionsadded.RData')
}

########################################
#Adding New Clusternames
#Change Idents to ArchR_Clusters and Sort them 

Idents(object = all) <- "ArchR_Clusters"
my_levels <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25")
Idents(all) <- factor(Idents(all), levels= my_levels)

new.cluster.ids <- c("C1", "Macrophage_C2", "Monocytes_C3", "CD8TCells_C4", "BCells_C5", "atrialCM_C6", 
                     "CM2_C7", "FGF13+CM_C8", "ventricularCM_C9", "C10", "C11", 
                     "C12", "ET1_C13", "mainET2_C14", "ET3_C15", "venousET4_C16", 
                     "lymphaticET5_C17", "VSMC1_C18", "VSMC2_C19", "mainVSMC3_C20", "FB1_C21",
                     "FB2_C22", "FB3_C23", "CM5?_C24", "C25")
names(new.cluster.ids) <- levels(all)
all <- RenameIdents(all, new.cluster.ids)
DimPlot(all, reduction = "harmony_archr", label = TRUE, pt.size = 0.9, repel = TRUE) + NoLegend()



cellnames <- as.character(all@active.ident)
all <- AddMetaData(all, metadata = cellnames, col.name = 'clusternames')


Idents(object = all) <- "clusternames"
my_levels <- c("C1", "Macrophage_C2", "Monocytes_C3", "CD8TCells_C4", "BCells_C5", "atrialCM_C6", 
               "CM2_C7", "FGF13+CM_C8", "ventricularCM_C9", "C10", "C11", 
               "C12", "ET1_C13", "mainET2_C14", "ET3_C15", "venousET4_C16", 
               "lymphaticET5_C17", "VSMC1_C18", "VSMC2_C19", "mainVSMC3_C20", "FB1_C21",
               "FB2_C22", "FB3_C23", "CM5?_C24", "C25")
Idents(all) <- factor(Idents(all), levels= my_levels)



#We next perform pre-processing and dimensional reduction on both assays independently, using standard approaches for RNA and ATAC-seq data.

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(all) <- "ATAC"
all <- RunTFIDF(all)
all <- FindTopFeatures(all, min.cutoff = 'q0')
all <- RunSVD(all)
all <- RunUMAP(all, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# RNA analysis
DefaultAssay(all) <- "RNA"
all <- SCTransform(all, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

#save.image(file='/Volumes/415-5899036/SingleCellSequencing/Seurat/AllAssaysanalized.RData')

#PCA Analysis and Determination of Sensible Clustering
all <- RunPCA(all, features = VariableFeatures(object = all), verbose = TRUE)


#We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering
#all <- RunHarmony(all, reduction = "wnn.umap",  plot_convergence = TRUE, group.by.vars = "orig.ident", assay.use="SCT")
#all <- RunHarmony(all, reduction = "wnn.umap",  plot_convergence = TRUE, group.by.vars = "orig.ident", reduction.save = "harmony")  #project.dim = F

all <- FindMultiModalNeighbors(all, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
all <- RunUMAP(all, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
all <- FindClusters(all, graph.name = "wsnn", algorithm = 3, verbose = TRUE, resolution=1.2, max.clusters = 25) #was resolution=1.2


#############was here


{


DimPlot(all, reduction = "pca")
DimPlot(all, reduction = "harmony_archr")  
DimPlot(all, reduction = "wnn.umap")  

# Examine and visualize PCA results a few different ways
print(all[["pca"]], dims = 1:25, nfeatures = 5)
print(all[["harmony_archr"]], dims = 1:25, nfeatures = 5)

VizDimLoadings(all, dims = 1:2, reduction = "harmony_archr")

DimHeatmap(all, dims = 1, cells = 500, balanced = TRUE, reduction = "harmony_archr")
DimHeatmap(all, dims = 6:9, cells = 1000, balanced = TRUE)


}


#Finding differentially expressed features
####################################
########################################

Idents(object = all) <- "ArchR_Clusters"
my_levels <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25")
Idents(all) <- factor(Idents(all), levels= my_levels)

Idents(object = all) <- "orig.ident"


#Set groups - join both Controls fror DEG
all$groups <- plyr::mapvalues(
  x = all@meta.data[["orig.ident"]], 
  from = c("MHC-CTR", "MHC", "TNT-CTR", "TNT"), 
  to = c("CTR", "MHC", "CTR", "TNT")
)

#all$ClusterOneControl <- paste(all@meta.data[["ArchR_Clusters"]],"_x_",all@meta.data[["groups"]])

all@meta.data$ClusterOneControl <- paste0(all@meta.data$ArchR_Clusters, "_x_", all@meta.data$groups)
Idents(object = all) <- all@meta.data$ClusterOneControl
all

#all<- AddMetaData(all, metadata=groups)

FeaturePlot(all, features = c("Tnnt2", "Hbb-bs", "Col3a1", "Col1a1", "-1", "ACTA", "FGF13", 
                                          "FGF11", "Nppa"), min.cutoff = "q9", reduction ="harmony_archr")

#Idents(object = all) <- all@meta.data$ClusterbySample

DefaultAssay(all) <- "RNA"

all <- NormalizeData(object = all)
all <- ScaleData(all)


save.image(file='/Volumes/415-5899036/SingleCellSequencing/Seurat_clusteringdone.RData')

all[["RNA"]]@counts<-as.matrix(all[["RNA"]]@counts)+1



# Differential Geneexpression of TNT vs CTR_joined
data.dir <- './DGE/'
dir.create(data.dir)
setwd(data.dir)
for (n in 16:25){
  #Geneexpression TNTvsCT
  #  value <- paste0(groups[[n]])
  ident1<- paste0("C",n,"_x_TNT")  
  ident2<- paste0("C",n,"_x_CTR")
  filem <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/","RNA_ClusterMarker_TNT_vs_CONTROLS_C",n,".xlsx")  
  #Volcano <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/","Volcanoplot_TNT_vs_CONTROLS_C_",n,".png")
  
  #motifsenriched <- paste0("/Pairwise/New/","Markers-Motifs-Enriched_C",n,"_x_TNT","vs_CTR")
  #motifsUpfile <- paste0("//Users/eliascrapa/ArchR/All/cs17-out/Plots//Pairwise/New/","Table_Markermotifs_up_Enriched_C",n,"_x_TNT","vs_CTR.csv")
  #motifsDofile <- paste0("//Users/eliascrapa/ArchR/All/cs17-out/Plots/Pairwise/New/","Table_Markermotifs_down_Enriched_C",n,"_x_TNT","vs_CTR.csv")
  #print(value)
  
  # find all markers of cluster 9
  
  cluster.markers <- FindMarkers(all, ident.1= ident1, ident.2 = ident2 , min.cells.feature = 5, test.use = "DESeq2", slot = "counts", assay = "RNA")
  names <- rownames(cluster.markers)
  
  
  
  # png device
#  png(Volcano)
  
#  print(EnhancedVolcano(cluster.markers,
#                        lab = rownames(cluster.markers),
#                        x = 'avg_log2FC',
#                        y = 'p_val',
#                        #pCutoff = 0.05,
#                       FCcutoff = 0.5,
#                        xlim = -1.5,
#                        ylim = 0.1)
        
#  )
  # Close device
#dev.off()
  
  
  
  rownames(cluster.markers) <- NULL
  cluster.markers <- cbind(names,cluster.markers)
  
  #cluster.markers$rownames <- rownames(cluster.markers)
  write.xlsx(cluster.markers, file = filem)
  #head(cluster22.markers, n = 25, wt = avg_log2FC)   
  
  
}

# Differential Geneexpression of MHC vs CTR_joined
for (n in 11:25){
  #Geneexpression MH_CvsCT
  #  value <- paste0(groups[[n]])
  ident1<- paste0("C",n,"_x_MHC")  
  ident2<- paste0("C",n,"_x_CTR")
  filem <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/","RNA_ClusterMarker_MHC_vs_CONTROLS_C",n,".xlsx")  
  #Volcano <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/Volcanoplots/DESeq2/","Volcanoplot_MHC_vs_CONTROLS_C_",n,".png")
  
  # find all markers of cluster 9
  
  cluster.markers <- FindMarkers(all, ident.1= ident1, ident.2 = ident2 , min.cells.feature = 5, test.use = "DESeq2", slot = "counts", assay = "RNA")
  names <- rownames(cluster.markers)
  
  
  

  
  
  
  rownames(cluster.markers) <- NULL
  cluster.markers <- cbind(names,cluster.markers)
  
  #cluster.markers$rownames <- rownames(cluster.markers)
  write.xlsx(cluster.markers, file = filem)
  #head(cluster22.markers, n = 25, wt = avg_log2FC)   
  
  
  
}

# Differential Geneexpression of TNT vs MHC
data.dir <- './TNTvsMHC/'
dir.create(data.dir)
setwd(data.dir)
for (n in 9:9){

  #Geneexpression MH_CvsCT
  #  value <- paste0(groups[[n]])
  ident1<- paste0("C",n,"_x_TNT")  
  ident2<- paste0("C",n,"_x_MHC")
  filem <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsMHC/","RNA_ClusterMarker_TNT_vs_MHC_C",n,".xlsx")  
  #Volcano <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/Volcanoplots/DESeq2/","Volcanoplot_MHC_vs_CONTROLS_C_",n,".png")
  
  # find all markers of cluster 9
  
  cluster.markers <- FindMarkers(all, ident.1= ident1, ident.2 = ident2 , min.cells.feature = 5, test.use = "DESeq2", slot = "counts", assay = "RNA")
  names <- rownames(cluster.markers)
  
  
  
  
  
  
  
  rownames(cluster.markers) <- NULL
  cluster.markers <- cbind(names,cluster.markers)
  
  #cluster.markers$rownames <- rownames(cluster.markers)
  write.xlsx(cluster.markers, file = filem)
  #head(cluster22.markers, n = 25, wt = avg_log2FC)   
  
  
  
}
setwd("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All")


##TNT DE with TNT-CTR
Idents(object = all) <- all@meta.data[["ClusterbySample"]]
for (n in 1:1){
  #Geneexpression TNTvsCT
  #  value <- paste0(groups[[n]])
  ident1<- paste0("C",n,"_x_TNT")  
  ident2<- paste0("C",n,"_x_CTR-TNT")
  filem <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/","RNA_ClusterMarker_TNT_vs_CTR-TNT_C",n,".xlsx")  
  #Volcano <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/","Volcanoplot_TNT_vs_CONTROLS_C_",n,".png")
  
  #motifsenriched <- paste0("/Pairwise/New/","Markers-Motifs-Enriched_C",n,"_x_TNT","vs_CTR")
  #motifsUpfile <- paste0("//Users/eliascrapa/ArchR/All/cs17-out/Plots//Pairwise/New/","Table_Markermotifs_up_Enriched_C",n,"_x_TNT","vs_CTR.csv")
  #motifsDofile <- paste0("//Users/eliascrapa/ArchR/All/cs17-out/Plots/Pairwise/New/","Table_Markermotifs_down_Enriched_C",n,"_x_TNT","vs_CTR.csv")
  #print(value)
  
  # find all markers of cluster 9
  
  cluster.markers <- FindMarkers(all, ident.1= ident1, ident.2 = ident2 , min.cells.feature = 5, test.use = "DESeq2", slot = "counts", assay = "RNA")
  names <- rownames(cluster.markers)
  
  
  
  # png device
  #  png(Volcano)
  
  #  print(EnhancedVolcano(cluster.markers,
  #                        lab = rownames(cluster.markers),
  #                        x = 'avg_log2FC',
  #                        y = 'p_val',
  #                        #pCutoff = 0.05,
  #                       FCcutoff = 0.5,
  #                        xlim = -1.5,
  #                        ylim = 0.1)
  
  #  )
  # Close device
  #dev.off()
  
  
  
  rownames(cluster.markers) <- NULL
  cluster.markers <- cbind(names,cluster.markers)
  
  #cluster.markers$rownames <- rownames(cluster.markers)
  write.xlsx(cluster.markers, file = filem)
  #head(cluster22.markers, n = 25, wt = avg_log2FC)   
  
  
}

all$celltype <- Idents(all)

# png device
png("UMAP-Plot-RNA_harmony5.png", width = 1600, height = 1280, units = "px", pointsize = 12)
#We can visualize clustering based on gene expression, ATAC-seq, or WNN analysis. The differences are more subtle than in the previous analysis (you can explore the weights, which are more evenly split than in our CITE-seq example), but we find that WNN provides the clearest separation of cell states
p1 <- DimPlot(all, reduction = "umap.rna", group.by = "ArchR_Clusters", label = TRUE, label.size = 0.05, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(all, reduction = "umap.atac", group.by = "orig.ident", label = TRUE, label.size = 0.05, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(all, reduction = "wnn.umap", group.by = "ArchR_Clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")
p4 <- DimPlot(all, reduction = "harmony_archr", group.by = "ArchR_Clusters", label = TRUE, label.size = 5, repel = F) + ggtitle("Harmony")
p4
p1 + p2 + p3 + p4  & NoLegend() & theme(plot.title = element_text(hjust = 1))
# Close device
dev.off()

# png device
png("Plotdifferent_groups.png", width = 960, height = 480, units = "px", pointsize = 12)
p5 <- DimPlot(all, reduction = "harmony_archr", split.by = "orig.ident", label = F, label.size = 1, repel = T) + ggtitle("Harmony")
p5
# Close device
dev.off()

# png device
png("UMAP-Plot-ATAC.png")
p2  & NoLegend() & theme(plot.title = element_text(hjust = 1))
# Close device
dev.off()

# png device
png("UMAP-Plot-WNN.png")
p3  & NoLegend() & theme(plot.title = element_text(hjust = 1))
# Close device
dev.off()

#For example, the ATAC-seq data assists in the separation of CD4 and CD8 T cell states. This is due to the presence of multiple loci that exhibit differential accessibility between different T cell subtypes. For example, we can visualize ‘pseudobulk’ tracks of the CD8A locus alongside violin plots of gene expression levels, using tools in the Signac visualization vignette.

## to make the visualization easier, subset T cell clusters
celltype.names <- levels(all)
tcell.names <- grep("1|2|3", celltype.names,value = TRUE)
tcells <- subset(all, idents = tcell.names)

# png device
png("CoveragePlot.png")

CoveragePlot(tcells,region = '1', features = '1', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE)

# Close device
dev.off()

saveRDS(all, file = "../output/TNT.rds")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#Next, we will examine the accessible regions of each cell to determine enriched motifs. As described in the Signac motifs vignette, there are a few ways to do this, but we will use the chromVAR package from the Greenleaf lab. This calculates a per-cell accessibility score for known motifs, and adds these scores as a third assay (chromvar) in the Seurat object.

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(RPresto)

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(all) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(all), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
all <- SetAssayData(all, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
all <- RunChromVAR(
  object = all,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

#Finally, we explore the multimodal dataset to identify key regulators of each cell state. Paired data provides a unique opportunity to identify transcription factors (TFs) that satisfy multiple criteria, helping to narrow down the list of putative regulators to the most likely candidates. We aim to identify TFs whose expression is enriched in multiple cell types in the RNA measurements, but also have enriched accessibility for their motifs in the ATAC measurements.
# As an example and positive control, the CCAAT Enhancer Binding Protein (CEBP) family of proteins, including the TF CEBPB, have been repeatedly shown to play important roles in the differentiation and function of myeloid cells including monocytes and dendritic cells. We can see that both the expression of the CEBPB, and the accessibility of the MA0466.2.4 motif (which encodes the binding site for CEBPB), are both enriched in monocytes.


#returns MA0466.2
motif.name <- ConvertMotifID(all, name = 'GATA4')
gene_plot <- FeaturePlot(all, features = "sct_GATA4", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(all, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
# png device
png("GeneMotifPlot.png")

gene_plot | motif_plot

# Close device
dev.off()

# We’d like to quantify this relationship, and search across all cell types to find similar examples. To do so, we will use the presto package to perform fast differential expression. We run two tests: one using gene expression data, and the other using chromVAR motif accessibilities. presto calculates a p-value based on the Wilcox rank sum test, which is also the default test in Seurat, and we restrict our search to TFs that return significant results in both tests.
# presto also calculates an “AUC” statistic, which reflects the power of each gene (or motif) to serve as a marker of cell type. A maximum AUC value of 1 indicates a perfect marker. Since the AUC statistic is on the same scale for both genes and motifs, we take the average of the AUC values from the two tests, and use this to rank TFs for each cell type:

markers_rna <- presto:::wilcoxauc.Seurat(X = all, group_by = 'celltype', assay = 'data', seurat_assay = 'SCT')
markers_motifs <- presto:::wilcoxauc.Seurat(X = all, group_by = 'celltype', assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(all, id = motif.names)

# a simple function to implement the procedure above
topTFs <- function(celltype, padj.cutoff = 1e-2) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
    arrange(-motif.auc)
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 11, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)
}


# We can now compute, and visualize, putative regulators for any cell type. We recover well-established regulators, including TBX21 for NK cells, IRF4 for plasma cells, SOX4 for hematopoietic progenitors, EBF1 and PAX5 for B cells, IRF8 and TCF4 for pDC. We believe that similar strategies can be used to help focus on a set of putative regulators in diverse systems.
# identify top markers in NK and visualize
head(topTFs("NK"), 3)

motif.name <- ConvertMotifID(all, name = 'ALX3')
gene_plot <- FeaturePlot(all, features = "sct_ALX3", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(all, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')

# png device
png("GeneMotifPlot2.png")

gene_plot | motif_plot

# Close device
dev.off()

# identify top markers in pDC and visualize
head(topTFs("pDC"), 3)

motif.name <- ConvertMotifID(all, name = 'TCF4')
gene_plot <- FeaturePlot(all, features = "sct_TCF4", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(all, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot

# identify top markers in HSPC and visualize
head(topTFs("CD16 Mono"),3)

motif.name <- ConvertMotifID(all, name = 'SPI1')
gene_plot <- FeaturePlot(all, features = "sct_SPI1", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(all, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')


# png device
png("GeneMotifPlot3.png")

gene_plot | motif_plot

# Close device
dev.off()

# identify top markers in other cell types
head(topTFs("Naive B"), 3)

head(topTFs("HSPC"), 3)

head(topTFs("Plasma"), 3)


#saving sessioninfo
library(sessioninfo)
sess <- session_info()
writeLines(capture.output(sess), paste0("sessioninfo_ArchR_",sess[["platform"]][["date"]],".txt"))



all

DefaultAssay(all) <- "RNA"
all <- SCTransform(all, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
