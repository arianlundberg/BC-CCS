## loading required data
######
load(file="~/manuscript_data.RData")
source(file = "~/functions.R")

### select genes with higher correlation +/-0.3
high.corr.genes <-corrs.order[abs(corrs.order$estimate)>0.3,]
#"ADRB2"     "EPCAM"     "MELTF"     "MEX3A"     "PLEKHG4"   "LINC01224" "ZNF93"     "RAD9B"     "E2F7"      "TCF19"     "C19orf57"

eset.high.corr <- eSet.final[rownames(eSet.final)%in%high.corr.genes$var2,eSet.final$twoType%in% 'Tumour']
pData(eset.high.corr) <- pData(eset.high.corr)[c('twoClass.code','tissue.colors','CCS')]
pData(eset.high.corr) <- merge(pData(eset.high.corr),CCS_relative.df[c('CCS_change','rownames')],by.x='row.names',by.y='rownames',sort=F,all.x=T)
rownames(pData(eset.high.corr)) <- eset.high.corr$Row.names;pData(eset.high.corr) <- pData(eset.high.corr)[,-1]
eset.high.corr$CCS_change.log2 <- log2(eset.high.corr$CCS_change)
eset.high.corr$groups <- as.factor(ifelse(eset.high.corr$twoClass.code %in% c("HNSC", "KICH", "KIRP", "UCEC"),'Group 1',
                                          ifelse(eset.high.corr$twoClass.code %in% c("CESC", "OV", "UCS"),'Group 2',NA)))

### generate the corr score - negative correlation * -1

corr.exprs <- eSet.final[match(high.corr.genes$var2,rownames(eSet.final)),eSet.final$twoType%in% 'Tumour']
#### eset of only 2groups

eset.high.corr.2group <- eset.high.corr[,!is.na(eset.high.corr$groups)]

#### reorder matrix based on groups

pData(eset.high.corr.2group) <- with(pData(eset.high.corr.2group), pData(eset.high.corr.2group)[order(groups),])
exprs(eset.high.corr.2group) <- exprs(eset.high.corr.2group)[,rownames(pData(eset.high.corr.2group))]

colData.heatmap.2group <- data.frame(groups=eset.high.corr.2group$groups,
                                     Cancer.type=as.character(eset.high.corr.2group$twoClass.code),
                                     CCS=as.character(eset.high.corr.2group$CCS),
                                     CCS_relative=rescale(eset.high.corr.2group$CCS_change.log2,to = c(0,100)))
colData.heatmap.2group$CCS <- factor(colData.heatmap.2group$CCS,labels = c('Low','Inter','High'))

tmp <- rbind.data.frame(ARDB2=exprs(eset.high.corr.2group)[which(rownames(eset.high.corr.2group)=='ADRB2'),]*-1,exprs(eset.high.corr.2group)[which(rownames(eset.high.corr.2group)!='ADRB2'),])
eset.high.corr.2group$corr.score <- colSums(tmp)

### determining the best cutpoint based on maximizes the product of Sensitivity and Specificity or Accuracy Area

cutpoint_C3 <- OptimalCutpoints::optimal.cutpoints(data=pData(eset.high.corr.2group),X = 'corr.score' ,methods='MaxProdSpSe',status ='groups',tag.healthy = 'Group 2')
cutpoint_C3 <- cutpoint_C3$MaxProdSpSe$Global$optimal.cutoff$cutoff

colData.heatmap.2group$corr.score.bin <- as.factor(ifelse(eset.high.corr.2group$corr.score > cutpoint_C3,1,0))

### generating CC core score based on the genes with high correlation value in the whole samples

tmp.all <- rbind.data.frame(ARDB2=exprs(eset.high.corr)[which(rownames(eset.high.corr)=='ADRB2'),]*-1,exprs(eset.high.corr)[which(rownames(eset.high.corr)!='ADRB2'),])

eset.high.corr$corr.score <- colSums(tmp.all)

### using two groups cut-point to binarize the data in all samples
eset.high.corr$corr.score.bin <- as.factor(ifelse(eset.high.corr$corr.score > cutpoint_C3,1,0))

CCS_relative.color = circlize::colorRamp2(seq(0, 100, by=1), (plasma(101)))

Data_col.heatmap.2group = list('Cancer types'=c("ACC"="#C1A72F","BLCA"="#FAD2D9","LGG"="#D49DC7","GBM"="#D49DC7",
                                                "BRCA"="#ED2891","CESC"="#F6B667","COAD"="#9EDDF9","UCEC"="#FBE3C7",
                                                "ESCA"="#007EB5","HNSC"="#97D1A9","KICH"="#ED1C24","KIRC"="#ED1C24","KIRP"="#ED1C24","LIHC"="#CACCDB",
                                                "LUAD"="#D3C3E0","LUSC"="#D3C3E0","OV"="#D97D25","PAAD"="#6E7BA2","PRAD"="#7E1918","SKCM"="#BBD642",
                                                "STAD"="#00AEEF","THCA"="#F9ED32","UCS"="#F89420"),
                               Groups=c("Group 1"="tomato4","Group 2"='#F89420'),
                               'C3 score'=c('Low'='gray29','High'='gray87'),
                               'CCS relative'=CCS_relative.color,
                               CCS=c("Low"="Black","Inter"="gray","High"="gold"))

exprs.data.2group <- as.matrix(eset.high.corr.2group)

#######
colData.annot.2group = HeatmapAnnotation('Cancer types' = colData.heatmap.2group$Cancer.type,
                                         'CCS relative'=colData.heatmap.2group$CCS_relative,
                                         'C3 score' =ifelse(colData.heatmap.2group$corr.score.bin=='0','Low','High'),
                                         col = Data_col.heatmap.2group,show_legend =T,show_annotation_name = F,na_col='white',
                                         annotation_legend_param = list('CCS relative' = list(
                                             title = "CCS relative",direction='horizontal',
                                             at = seq(0,100,by=10),
                                             labels = seq(0,100,by=10)),
                                             'Cancer types'=list(ncol=3,title='Cancer types')))


##### all other cancers

eset.high.corr.all_others <- eset.high.corr[,is.na(eset.high.corr$groups)]
tmp.all_others <- rbind.data.frame(ARDB2=exprs(eset.high.corr.all_others)[which(rownames(eset.high.corr.all_others)=='ADRB2'),]*-1,exprs(eset.high.corr.all_others)[which(rownames(eset.high.corr.all_others)!='ADRB2'),])
eset.high.corr.all_others$corr.score <- colSums(tmp.all_others)

#### reorder matrix based on groups

pData(eset.high.corr.all_others) <- with(pData(eset.high.corr.all_others), pData(eset.high.corr.all_others)[order(twoClass.code),])
exprs(eset.high.corr.all_others) <- exprs(eset.high.corr.all_others)[,rownames(pData(eset.high.corr.all_others))]
colData.heatmap.all_others <- data.frame(Cancer.type=as.character(eset.high.corr.all_others$twoClass.code),
                                         CCS=as.character(eset.high.corr.all_others$CCS),
                                         CCS_relative=rescale(eset.high.corr.all_others$CCS_change.log2,to = c(0,100)))
colData.heatmap.all_others$CCS <- factor(colData.heatmap.all_others$CCS,labels = c('Low','Inter','High'))
colData.heatmap.all_others$corr.score.bin <- as.factor(ifelse(eset.high.corr.all_others$corr.score > cutpoint_C3,1,0))

Data_col.heatmap.all_others = list('Cancer types'=c("ACC"="#C1A72F","BLCA"="#FAD2D9","LGG"="#D49DC7","GBM"="#D49DC7",
                                                    "BRCA"="#ED2891","CESC"="#F6B667","COAD"="#9EDDF9","UCEC"="#FBE3C7",
                                                    "ESCA"="#007EB5","HNSC"="#97D1A9","KICH"="#ED1C24","KIRC"="#ED1C24","KIRP"="#ED1C24","LIHC"="#CACCDB",
                                                    "LUAD"="#D3C3E0","LUSC"="#D3C3E0","OV"="#D97D25","PAAD"="#6E7BA2","PRAD"="#7E1918","SKCM"="#BBD642",
                                                    "STAD"="#00AEEF","THCA"="#F9ED32","UCS"="#F89420"),
                                   Groups=c("Group 1"="tomato4","Group 2"='#F89420'),
                                   'C3 score'=c('Low'='gray29','High'='gray87'),
                                   'CCS relative'=CCS_relative.color,
                                   CCS=c("Low"="Black","Inter"="gray","High"="gold"))

exprs.data.all_others <- as.matrix(eset.high.corr.all_others)

#######
colData.annot.all_others = HeatmapAnnotation('Cancer types' = colData.heatmap.all_others$Cancer.type,
                                             'CCS relative'=colData.heatmap.all_others$CCS_relative,
                                             'C3 score' =ifelse(colData.heatmap.all_others$corr.score.bin=='0','Low','High'),
                                             col = Data_col.heatmap.all_others,show_legend =T,show_annotation_name = T,na_col='white',
                                             annotation_legend_param = list('CCS relative' = list(
                                                 title = "CCS relative",direction='horizontal',
                                                 at = seq(0,100,by=10),
                                                 labels = seq(0,100,by=10)),
                                                 'Cancer types'=list(ncol=3,title='Cancer types')))

### generating correlation bars (right annotations)
corr.data.all_others <- merge(t(exprs(eset.high.corr.all_others)),pData(eset.high.corr.all_others)['CCS_change.log2'],by='row.names',sort=F)
rownames(corr.data.all_others) <- corr.data.all_others$Row.names;corr.data.all_others <- corr.data.all_others[,-1]

corrbar.all_others <- cor(corr.data.all_others[,-12],corr.data.all_others$CCS_change.log2,use = 'complete.obs')

#### merge heatmaps into one figure
### Figure 5A

set.seed(1363)

heatmap_list <- Heatmap(matrix = exprs.data.2group,name = 'log2(TPM+1)',border = T,cluster_rows = T,heatmap_legend_param = list(direction='horizontal'),
                        cluster_columns = T,top_annotation = colData.annot.2group,row_names_gp = gpar(fontsize=10),row_dend_reorder = T,show_row_names = F,show_column_names = F,show_column_dend = T,
                        column_names_gp = gpar(fontsize=6),column_title="Two groups") +
    Heatmap(matrix = exprs.data.all_others,name = 'log2(TPM+1)',border = T,cluster_rows = T,heatmap_legend_param = list(direction='horizontal'),
            cluster_columns = T,top_annotation = colData.annot.all_others,row_names_gp = gpar(fontsize=10),row_dend_reorder = T,show_row_names = T,show_column_names = F,show_column_dend = T,
            right_annotation = rowAnnotation("Pearson's r"=anno_barplot(corrbar.all_others,border = T,width = unit(3,'cm'),
                                                                        ylim = c(-0.3,0.5),which = 'row',baseline = 0,height = unit(1.5,'cm'),axis = T,bar_width = 0.8,
                                                                        gp=gpar(border=NA,lty='blank',fill='skyblue3',color='skyblue3'))),
            column_names_gp = gpar(fontsize=6),column_title="All-other cancer types")

corr_heatmap <- draw(heatmap_list,
                     column_title_gp = gpar(fontsize = 14,type='bold'),padding = unit(c(2, 8, 2, 2), "mm"),heatmap_legend_side = "bottom",annotation_legend_side='bottom')
decorate_annotation("Pearson's r", {
    grid.lines(c(0, 0), c(0.5, 11.5), gp = gpar(lty = 2),
               default.units = "native")
})


### selecting top 5 amp/del locations
top_locs <- data.frame(chr=c('3q','1q','5p','6p','2p',
                             '16q','8p','9q','22 (22q)','4q'),color=c(rep('firebrick3',5),rep('gold',5)))

anneuploidy.score_sub <- anneuploidy.score %>% select(Sample,top_locs$chr)

### generating matrix for the binary heatmap

colData.heatmap <- data.frame(Cancer.type=as.character(eset.high.corr$twoClass.code),
                              corr.score.bin=eset.high.corr$corr.score.bin,
                              CCS_relative=rescale(eset.high.corr$CCS_change.log2,to=c(0,100)))

rownames(colData.heatmap) <- colnames(eset.high.corr)

Data_col.heatmap = list('Cancer types'=c("ACC"="#C1A72F","BLCA"="#FAD2D9","LGG"="#D49DC7","GBM"="#D49DC7",
                                         "BRCA"="#ED2891","CESC"="#F6B667","COAD"="#9EDDF9","UCEC"="#FBE3C7",
                                         "ESCA"="#007EB5","HNSC"="#97D1A9","KICH"="#ED1C24","KIRC"="#ED1C24","KIRP"="#ED1C24","LIHC"="#CACCDB",
                                         "LUAD"="#D3C3E0","LUSC"="#D3C3E0","OV"="#D97D25","PAAD"="#6E7BA2","PRAD"="#7E1918","SKCM"="#BBD642",
                                         "STAD"="#00AEEF","THCA"="#F9ED32","UCS"="#F89420"),
                        Groups=c("Group 1"="tomato4","Group 2"='#F89420'),
                        'C3 score'=c('Low'='gray29','High'='gray87'),
                        'CCS relative'=CCS_relative.color,
                        CCS=c("Low"="Black","Inter"="gray","High"="gold"))


annu.matrix <- anneuploidy.score_sub;rownames(annu.matrix) <- anneuploidy.score_sub$Sample;annu.matrix <- annu.matrix[,-1]
annu.matrix <- annu.matrix[which(rownames(annu.matrix)%in%rownames(colData.heatmap)),]
colData.heatmap <- colData.heatmap[which(rownames(colData.heatmap)%in%rownames(annu.matrix)),]
annu.matrix <- annu.matrix[rownames(colData.heatmap),]
annu.matrix <- t(annu.matrix);annu.matrix <- as.matrix(annu.matrix)

### remove rows or columns with more than 50% missing
annu.matrix <- annu.matrix[which(rowMeans(!is.na(annu.matrix)) > 0.5), which(colMeans(!is.na(annu.matrix)) > 0.5)]
colData.heatmap <- colData.heatmap[colnames(annu.matrix),]

colData.heatmap.annu.2group <- colData.heatmap[colData.heatmap$Cancer.type %in% c("HNSC", "KICH", "KIRP", "UCEC","CESC", "OV", "UCS"),]
colData.annot.annu.2group = HeatmapAnnotation('Cancer types' = colData.heatmap.annu.2group$Cancer.type,
                                              'CCS relative'=colData.heatmap.annu.2group$CCS_relative,
                                              'C3 score'=ifelse(colData.heatmap.annu.2group$corr.score.bin=='0','Low','High'),
                                              col = Data_col.heatmap,show_legend =T,show_annotation_name = F,na_col='white',
                                              annotation_legend_param = list('CCS relative' = list(
                                                  title = "CCS relative",direction='horizontal',
                                                  at = seq(0,100,by=10),
                                                  labels = seq(0,100,by=10)),
                                                  'Cancer types'=list(ncol=3,title='Cancer types')))
annu.matrix.2group <- annu.matrix[,rownames(colData.heatmap.annu.2group)]


colData.heatmap <- colData.heatmap[rownames(colData.heatmap) %out% rownames(colData.heatmap.annu.2group),]
annu.matrix <- annu.matrix[,colnames(annu.matrix) %out% rownames(colData.heatmap.annu.2group)]

colData.annot = HeatmapAnnotation('Cancer types' = colData.heatmap$Cancer.type,
                                  'CCS relative'=colData.heatmap$CCS_relative,
                                  'C3 score'=ifelse(colData.heatmap$corr.score.bin=='0','Low','High'),
                                  col = Data_col.heatmap,show_legend =T,show_annotation_name = T,na_col='white',
                                  annotation_legend_param = list('CCS relative' = list(
                                      title = "CCS relative",direction='horizontal',
                                      at = seq(0,100,by=10),
                                      labels = seq(0,100,by=10)),
                                      'Cancer types'=list(ncol=3,title='Cancer types')))

row_split = data.frame(rep(c("Amplification", "Deletion"), each = 5))

### generating correlation bars (right annotations)
corr.data.annu.all_others <- merge(t(annu.matrix),pData(eset.high.corr.all_others)['CCS_change.log2'],by='row.names',sort=F)
rownames(corr.data.annu.all_others) <- corr.data.annu.all_others$Row.names;corr.data.annu.all_others <- corr.data.annu.all_others[,-1]

corrbar.annu.all_others <- cor(corr.data.annu.all_others[,-11],corr.data.annu.all_others$CCS_change.log2,use = 'complete.obs',method = 'kendall')





annu.heatmap_list <- Heatmap(matrix = annu.matrix.2group,col = structure(c('firebrick3','grey95','gold3','white'), names =c(1,0,-1,NA)),na_col = 'white',
                             name = 'Status',border = T,cluster_rows = T,heatmap_legend_param = list(labels = c("Gain", "No-change", "Loss","Missing-data")),
                             cluster_columns = T,top_annotation = colData.annot.annu.2group,row_split = row_split,row_names_gp = gpar(fontsize=10),row_dend_reorder = F,show_row_names = T,show_column_names = F,show_column_dend = T,
                             column_names_gp = gpar(fontsize=6),column_title="Two groups") +
    Heatmap(matrix = annu.matrix,col = structure(c('firebrick3','grey95','gold3','white'), names =c(1,0,-1,NA)),na_col = 'white',
            name = 'Status',border = T,cluster_rows = T,heatmap_legend_param = list(labels = c("Gain", "No-change", "Loss","Missing-data")),
            right_annotation = rowAnnotation("Kendall's tau"=anno_barplot(corrbar.annu.all_others,border = T,width = unit(2.5,'cm'),ylim = c(-0.1,0.05),which = 'row',baseline = 0,height = unit(1.5,'cm'),axis = T,bar_width = 0.8, gp=gpar(border=NA,lty='blank',fill='skyblue3',color='skyblue3'))),
            cluster_columns = T,top_annotation = colData.annot,row_split = row_split,row_names_gp = gpar(fontsize=10),row_dend_reorder = F,show_row_names = T,show_column_names = F,show_column_dend = T,
            column_names_gp = gpar(fontsize=6),column_title="All-other cancer types")



### Figure 5B


annu_heatmap <- draw(annu.heatmap_list,
                     column_title_gp = gpar(fontsize = 14,type='bold'),padding = unit(c(2, 8, 2, 2), "mm"),heatmap_legend_side = "bottom",annotation_legend_side='bottom')
decorate_annotation("Kendall's tau", {
    grid.lines(c(0, 0), c(-4.5, 0.45), gp = gpar(lty = 2),
               default.units = "native")
    grid.lines(c(0, 0), c(0.5, 5.5), gp = gpar(lty = 2),
               default.units = "native")
})



