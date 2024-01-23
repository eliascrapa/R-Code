library(ggplot2)
library(xlsx)
library(ggpubr)
library(tidyverse)

setwd("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/")
data.dir <- './IPA'
dir.create(data.dir)
setwd(data.dir)


cluster <- list()
cluster <- c("TNT_C2","TNT_C6","TNT_C9","TNT_C14","TNT_C16","TNT_C19","TNT_C20","TNT_C21","TNT_C22","MHC_C2","MHC_C6","MHC_C9","MHC_C14","MHC_C16","MHC_C20","MHC_C22")
cluster <- c("TNT_C9","MHC_C9")

for (n in cluster) {
  
tabledf <- read.xlsx(paste0("/Users/eliascrapa/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/IPA/",n,"_CanonicalPathway.xls"),1)

names(tabledf) <- as.matrix(tabledf[1, ])
tabledf <- tabledf[-1, ]
tabledf[] <- lapply(tabledf, function(x) type.convert(as.character(x)))
colnames(tabledf) <- c("Pathways", "-logpval", "ratio","zscore",  "Molecules")
tabledf <- tabledf[order(-tabledf$zscore),]
tabledf3 <- tabledf %>% filter(tabledf$`-logpval` > 1.3) %>%  filter(!is.na(zscore))
#tabledf2 <- head(tabledf, 40)

tabledf2 <- rbind((tabledf3 %>% slice_min(tabledf3$zscore,n = 8) %>% filter(zscore < 0)),(tabledf3 %>% slice_max(tabledf3$zscore,n = 10) %>% filter(zscore > 0)))

tabledf2$posneg <- ifelse(tabledf2$zscore > 0,"pos", "neg")

pdf(paste0("Pathwayanalysis_IPA_",n,".pdf"))#, width = 1200, height = 1800)
print(ggbarplot(tabledf2, x = "Pathways", y = "zscore",
          fill = "posneg",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = F,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "z-score",
          legend.title = paste0(n),
          rotate = T,
          ggtheme = theme_minimal(),
          font.tickslab  = 11
         # top = 10
))
dev.off()

}
