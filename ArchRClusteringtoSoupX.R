#Workflow to adjust ArchR Clustering to SoupX


#manual loading of Soupdata
#tod <- Read10X_h5("/Users/eliascrapa/ArchR/DATA/TNT/raw_feature_bc_matrix.h5")
#toc <- Read10X_h5("/Users/eliascrapa/ArchR/DATA/TNT/filtered_feature_bc_matrix.h5")
#tod <-tod[["Gene Expression"]]  #just use the genexpression matrix - loose the Peakmatrix
#toc <-toc[["Gene Expression"]]  #just use the genexpression matrix - loose the Peakmatrix
#sc = SoupChannel(tod,toc)


#export clustering
Idents(object = TNT) <- "ArchR_Clusters"
TNTcluster_ <- data.frame(TNT@active.ident)

rownames(TNTcluster_) <- gsub("TNT_","",rownames(TNTcluster_))
rowname <- rownames(TNTcluster_) 
TNTcluster_ <-cbind(rowname,TNTcluster_)
rm(rowname)
rownames(TNTcluster_) <- NULL
colnames(TNTcluster_) <- c("cellIDs", "clusterIDs")
TNTcluster_ <- as.matrix(TNTcluster_)
TNTcluster_ <-structure(as.character(TNTcluster_$clusterIDs),names=as.character(TNTcluster_$cellIDs))
sc = setClusters(sc, setNames(TNTcluster_$clusterIDs,TNTcluster_$cellIDs)) #setNames(clusters, barcodes)

data(PBMC_sc)

data(PBMC_metaData)

tmpDir = tempdir(check = TRUE)
tmpDir = ("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/tmpDir")
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmpDir, "tod.tar.gz"))
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmpDir, "toc.tar.gz"))
untar(file.path(tmpDir, "tod.tar.gz"), exdir = tmpDir)
untar(file.path(tmpDir, "toc.tar.gz"), exdir = tmpDir)
