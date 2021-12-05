load(file="~/Box Sync/Karolinska/Projects/CCS_2020/munge/Github/Gene.expression_data.RData")
source(file="~/Box Sync/Karolinska/Projects/CCS_2020/munge/Github/functions.R")
######## finding variable genes

variable.genes.data <- exprs(eSet.Data)

#### choosing the genes that have highest coefficient of variation (most variable genes)
cv <- apply(variable.genes.data, 1, compute_cv)

variable.genes.data <- variable.genes.data[rank(cv) / length(cv) > 1 - 0.25, ]
variable.genes <- rownames(variable.genes.data)

####
# eSet.final = eset of all genes without outliers : assayData: 25453 features, 13117 samples

variable.genes.noOutlier.data <- exprs(eSet.final)[variable.genes,]

########
# script to generate Figure 2
########

# Figure 2A, B - Variable genes
## PCA the data - outliers removed - Variable genes (Figure 2 C,D)

set.seed(1363)
PCA_variable.noOutlier <- gmodels::fast.prcomp(x = t(variable.genes.noOutlier.data))
PCA_variable.noOutlier.data <- as.data.frame(PCA_variable.noOutlier$x)


### calculating variance manually to use in the plot
dev.variable.noOutlier <- PCA_variable.noOutlier$sdev^2/sum(PCA_variable.noOutlier$sdev^2)
plot.study.noOutlier.colors <- eSet.final$PCA.tag
plot.tissue.noOutlier.colors <- eSet.final$primary_site_relevel


PCA.variable.noOutlier.study.plot <- ggplot(PCA_variable.noOutlier.data ,aes(x=PC1,y=PC2)) +
    geom_point(aes(color=(plot.study.noOutlier.colors)),size=0.2,alpha=0.7) + scale_color_jama(name = 'Study origin')+
    theme_few()+ guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.title = element_text(face = "bold")) +labs(title='PCA: Study origin',
                                                            x=paste0("PC1: ",round(dev.variable.noOutlier[1]*100,1),"%"),
                                                            y=paste0("PC2: ",round(dev.variable.noOutlier[2]*100,1),"%"))

PCA.variable.noOutlier.tissue.plot <- ggplot(PCA_variable.noOutlier.data ,aes(x=PC1,y=PC2)) +
    geom_point(aes(color=(plot.tissue.noOutlier.colors)),size=0.2,alpha=0.7) + scale_color_manual(name='Tissue type',values = Tissue.color$tissue.colors)+
    theme_few()+ guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.title = element_text(face = "bold")) +labs(title='PCA: Tissue type',
                                                            x=paste0("PC1: ",round(dev.variable.noOutlier[1]*100,1),"%"),
                                                            y=paste0("PC2: ",round(dev.variable.noOutlier[2]*100,1),"%"))


###  Figure 2A-B - PCA 3D plots

PCA.colors.noOutlier.studies <- pal_jama()(3)[eSet.final$PCA.tag]
PCA.colors.noOutlier.tissues <- eSet.final$tissue.colors


scatterplot3d::scatterplot3d(PCA_variable.noOutlier$x[,1:3],color = PCA.colors.noOutlier.studies,pch = 16, grid=TRUE, box=FALSE,type = 'p',angle = -30,cex.symbols = 0.5,
                             xlab =paste0("PC1: ",round(dev.variable.noOutlier[1]*100,1),"%"),
                             ylab =paste0("PC2: ",round(dev.variable.noOutlier[2]*100,1),"%"),
                             zlab =paste0("PC3: ",round(dev.variable.noOutlier[3]*100,1),"%"))


scatterplot3d::scatterplot3d(PCA_variable.noOutlier$x[,1:3],color = PCA.colors.noOutlier.tissues,pch = 16, grid=TRUE, box=FALSE,type = 'p',angle = -30,cex.symbols = 0.5,
                             xlab =paste0("PC1: ",round(dev.variable.noOutlier[1]*100,1),"%"),
                             ylab =paste0("PC2: ",round(dev.variable.noOutlier[2]*100,1),"%"),
                             zlab =paste0("PC3: ",round(dev.variable.noOutlier[3]*100,1),"%"))



########

# Figure 2C-D, CCS
### eSet.CCS.final = expression set with CCS genes (outliers removed)

CCS.final.genes.data <- exprs(eSet.CCS.final)

set.seed(1363)
PCA_CCS.final <- gmodels::fast.prcomp(x = t(CCS.final.genes.data))
PCA_CCS.final.data <- as.data.frame(PCA_CCS.final$x)

dev.CCS.final <- PCA_CCS.final$sdev^2/sum(PCA_CCS.final$sdev^2)

PCA.CCS.final.study.plot <- ggplot(PCA_CCS.final.data ,aes(x=PC1,y=PC2)) +
    geom_point(aes(color=(eSet.CCS.final$PCA.tag)),size=0.2,alpha=0.7) + scale_color_jama(name = 'Study origin')+
    theme_few()+ guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.title = element_text(face = "bold")) +labs(title='PCA: Study origin',
                                                            x=paste0("PC1: ",round(dev.CCS.final[1]*100,1),"%"),
                                                            y=paste0("PC2: ",round(dev.CCS.final[2]*100,1),"%"))


PCA.CCS.final.tissue.plot <- ggplot(PCA_CCS.final.data ,aes(x=PC1,y=PC2)) +
    geom_point(aes(color=(eSet.CCS.final$primary_site_relevel)),size=0.2,alpha=0.7) + scale_color_manual(name='Tissue type',values = Tissue.color.noTestis$tissue.colors)+
    theme_few()+ guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.title = element_text(face = "bold")) +labs(title='PCA: Tissue type',
                                                            x=paste0("PC1: ",round(dev.CCS.final[1]*100,1),"%"),
                                                            y=paste0("PC2: ",round(dev.CCS.final[2]*100,1),"%"))



########
# UMAPs
########
#UMAP - CCS
# Figure 2E, F, G


Exprs_CCS.noOut <- merge(t(exprs(eSet.CCS.final)),pData(eSet.CCS.final)[c('CCS_ct','primary_site_relevel','PCA.tag')],all.x=T,sort=F,by='row.names')
rownames(Exprs_CCS.noOut) <- Exprs_CCS.noOut$Row.names;Exprs_CCS.noOut <- Exprs_CCS.noOut[,-1]
Exprs_CCS.noOut$CCS_ct <- rescale(Exprs_CCS.noOut$CCS_ct,new.range=c(0,1))

umap.CCS.study.plot <- umap(t(Exprs_CCS.noOut[,-c(448:450)]),labels=Exprs_CCS.noOut$PCA.tag,controlscale=TRUE,scale=3,legendtitle = 'Study origin',dotsize=1,axistextsize=10,
                            seed = 1363,legendtextsize=8)+scale_color_jama(name='Study origin')+
    guides(colour = guide_legend(override.aes = list(size=5),nrow = 2)) +
    theme(legend.title = element_text(face = "bold"))+labs(title='UMAP: Study origin')

umap.CCS.tissue.plot <- umap(t(Exprs_CCS.noOut[,-c(448:450)]),labels=Exprs_CCS.noOut$primary_site_relevel,controlscale=TRUE,scale=3,legendtitle = 'Tissue type',dotsize=1,axistextsize=10,
                             seed = 1363,legendtextsize=8)+scale_color_manual(values = Tissue.color[-18,]$tissue.colors,name='Tissue type')+
    guides(colour = guide_legend(override.aes = list(size=5),nrow = 2)) +
    theme(legend.title = element_text(face = "bold"))+labs(title='UMAP: Tissue type')

umap.CCS.range.plot <- umap(t(Exprs_CCS.noOut[,-c(448:450)]),labels=Exprs_CCS.noOut$CCS_ct,controlscale=TRUE,scale = 1,low = 0,high = 1,legendtitle = 'CCS range',dotsize=1,axistextsize=10,
                            seed = 1363,legendtextsize=8)+scale_color_viridis_c(option="cividis",begin = 0,end = 1)+
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.title = element_text(face = "bold"))+labs(title='UMAP: CCS range')



###############
# Figure 2
# generate figures
###############

plot.legends <- plot_grid(get_legend(PCA.CCS.final.tissue.plot + guides(colour = guide_legend(override.aes = list(size=5),nrow = 7))+
                                         theme(legend.position = 'bottom',legend.title =element_text(size=12,face='bold'),legend.text =element_text(size=10),legend.direction = 'vertical')),
                          get_legend(PCA.CCS.final.study.plot + guides(colour = guide_legend(override.aes = list(size=8),nrow = 7))+
                                         theme(legend.position = 'bottom',legend.title =element_text(size=12,face='bold'),legend.text =element_text(size=10),legend.direction = 'vertical')),
                          get_legend(umap.CCS.range.plot+ guides(colour = guide_legend(override.aes = list(size=8),nrow = 7))+
                                         theme(legend.position = 'bottom',legend.title =element_text(size=12,face='bold'),legend.text =element_text(size=10),legend.direction = 'vertical')),
                          ncol = 1,nrow=3,axis='tblr',rel_widths = c(0.5,0.5,0.5))

#### Figure 2 layout

Figure2 <- plot_grid(PCA.variable.noOutlier.study.plot + theme_few(base_size = 22)+ theme(legend.position ='none',title=element_text(size=16))+
                         labs(title='PCA: Study origin',x=paste0("PC1: ",round(dev.variable.noOutlier[1]*100,1),"%"),y=paste0("Most variable genes\n(gene number = 6364)","\n\nPC2: ",round(dev.variable.noOutlier[2]*100,1),"%")),
                     PCA.variable.noOutlier.tissue.plot + theme_few(base_size = 22)+ theme(legend.position ='none',title=element_text(size=16)),
                     NULL,
                     PCA.CCS.final.study.plot + theme_few(base_size = 22) + theme(legend.position ='none',title=element_text(size=16))+
                         labs(title='PCA: Study origin',x=paste0("PC1: ",round(dev.CCS.final[1]*100,1),"%"),y=paste0("CCS genes\n(gene number = 446)","\n\nPC2: ",round(dev.CCS.final[2]*100,1),"%")),
                     PCA.CCS.final.tissue.plot + theme_few(base_size = 22) + theme(legend.position ='none',title=element_text(size=16)),
                     NULL,
                     umap.CCS.study.plot + theme_few(base_size = 22)+ theme(legend.position ='none',title=element_text(size=16))+
                         labs(title='UMAP: Study origin',y=paste0("CCS genes\n(gene number = 446)","\n\nX2")),
                     umap.CCS.tissue.plot + theme_few(base_size = 22)+ theme(legend.position ='none',title=element_text(size=16)),
                     umap.CCS.range.plot + theme_few(base_size = 22)+ theme(legend.position ='none',title=element_text(size=16)),
                     labels =c('A','B','',
                               'C','D','',
                               'E','F','G'),axis = 'tb',ncol = 3,nrow=3,
                     rel_widths = c(1.13,1,1,1,1,1),rel_heights = c(1,1,1,1))



