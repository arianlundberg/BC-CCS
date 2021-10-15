required.packages <- c("msigdbr", "limma", "viridis","scatterplot3d", "cowplot", "ggbiplot", "biomaRt", "reshape2", "scales", "M3C","gganatogram", "ggfortify",
                       "randomcoloR","forcats","maftools","rstatix","hrbrthemes","org.Hs.eg.db","RColorBrewer", "Rtsne", "plyr", "tidyverse", "stringr", "lubridate",
                       "data.table", "dplyr", "readxl", "foreach", "openxlsx", "GenomicRanges", "DESeq2","GenomeInfoDb", "doParallel", "ComplexHeatmap", "ggthemes",
                       "ggsci", "ggplot2", "ggformula", "ggstatsplot", "ggExtra", "ggpubr","DataCombine","doBy","gghalves")

lapply(required.packages,require, character.only = TRUE)

'%out%' <- Negate('%in%')
## function of coefficient variation
compute_cv <- function(x) sd(x) / mean(x)

### correlation function
corrFunc <- function(var1,var2,data){
    result=cor.test(data[,var1],data[,var2],method='spearman')
    data.frame(var1, var2, result[c("estimate","p.value","statistic","method")],
               stringsAsFactors=FALSE)
}

