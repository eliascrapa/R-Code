#Soup-X Contamination Filtering

install.packages("SoupX")
library(SoupX)
# Load data and estimate soup profile
sc = load10X("Path/to/cellranger/outs/folder/")
# Estimate rho
sc = autoEstCont(sc)
# Clean the data
out = adjustCounts(sc)