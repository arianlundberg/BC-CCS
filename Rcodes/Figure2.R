load(file="~/manuscript_data.RData")
source(file = "~/functions.R")

######## finding variable genes

variable.genes.data <- exprs(eSet.Data)

#### choosing the genes that have highest coefficient of variation (most variable genes)

cv <- apply(variable.genes.data, 1, compute_cv)

variable.genes.data <- variable.genes.data[rank(cv) / length(cv) > 1 - 0.25, ]
variable.genes <- rownames(variable.genes.data)
# length(rownames(variable.genes.data) )
# [1] 6364


########
# script to generate Figure 2 (not Consort here)
########

## PCA the data with outliers (testis samples, majority of them were outliers - all normal testis and two tumours) - Variable genes (Figure 2A)
set.seed(1363)
PCA_variable <- gmodels::fast.prcomp(x = t(variable.genes.data))
PCA_variable.data <- as.data.frame(PCA_variable$x)

#### these outliers selected based on PCA plot

PCA.outliers.idx <- PCA_variable.data$PC1 > 115
PCA.outliers <- PCA_variable.data[PCA.outliers.idx,]
Outlier.data <- pData(eSet.Data)[rownames(PCA.outliers),c('primary_site_relevel','PCA.tag')]
Outlier.data$primary_site_relevel <- factor(Outlier.data$primary_site_relevel)

#table(Outlier.data$primary_site_relevel,Outlier.data$PCA.tag)
#                     GTEX-Normal TCGA-Normal TCGA-Tumour
# Brain                 3           0           0
# Endometrium           0           0           1
# Esophagus             0           1           1
# Lung                  0           0           1
# Stomach               0           0           5
# Testis              156           0           2


### calculating variance manually to use in the plot
dev.variable <- PCA_variable$sdev^2/sum(PCA_variable$sdev^2)
plot.study.colors <- eSet.Data$PCA.tag

########

########
# Figure 2A - Variable genes with outlier

PCA.variable.study.plot <- ggplot(PCA_variable.data ,aes(x=PC1,y=PC2)) +
    geom_point(aes(color=(plot.study.colors)),size=0.2,alpha=0.7) + scale_color_jama(name = 'Study origin')+
    geom_point(data=PCA.outliers,size=0.2,color='red') +
    theme_few()+ guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.title = element_text(face = "bold")) +labs(title='PCA: Study origins*',
                                                            x=paste0("PC1: ",round(dev.variable[1]*100,1),"%"),
                                                            y=paste0("PC2: ",round(dev.variable[2]*100,1),"%"))


### PCA 3D plots

PCA.colors.studies <- pal_jama()(3)[eSet.Data$PCA.tag]
PCA.colors.studies[as.numeric(which(PCA.outliers.idx))] <- 'red'

#Figure 2A - PCA 3D plot

scatterplot3d::scatterplot3d(PCA_variable$x[,1:3],color = PCA.colors.studies,pch = 16, grid=TRUE, box=FALSE,type = 'p',angle = -35,cex.symbols = 0.5,
                             xlab =paste0("PC1: ",round(dev.variable[1]*100,1),"%"),
                             ylab =paste0("PC2: ",round(dev.variable[2]*100,1),"%"),
                             zlab =paste0("PC3: ",round(dev.variable[3]*100,1),"%"))



########

## PCA the data - outliers removed - Variable genes (Figure 2 C,D)

set.seed(1363)
PCA_variable.noOutlier <- gmodels::fast.prcomp(x = t(variable.genes.noOutlier.data))
PCA_variable.noOutlier.data <- as.data.frame(PCA_variable.noOutlier$x)


### calculating variance manually to use in the plot
dev.variable.noOutlier <- PCA_variable.noOutlier$sdev^2/sum(PCA_variable.noOutlier$sdev^2)
plot.study.noOutlier.colors <- eSet.final$PCA.tag
plot.tissue.noOutlier.colors <- eSet.final$primary_site_relevel

########


########
# Figure 2B,C - Variable genes

PCA.variable.noOutlier.tissue.plot <- ggplot(PCA_variable.noOutlier.data ,aes(x=PC1,y=PC2)) +
    geom_point(aes(color=(plot.tissue.noOutlier.colors)),size=0.2,alpha=0.7) + scale_color_manual(name='Tissue type',values = Tissue.color$tissue.colors)+
    theme_few()+ guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.title = element_text(face = "bold")) +labs(title='PCA: Tissue type',
                                                            x=paste0("PC1: ",round(dev.variable.noOutlier[1]*100,1),"%"),
                                                            y=paste0("PC2: ",round(dev.variable.noOutlier[2]*100,1),"%"))

PCA.variable.noOutlier.study.plot <- ggplot(PCA_variable.noOutlier.data ,aes(x=PC1,y=PC2)) +
    geom_point(aes(color=(plot.study.noOutlier.colors)),size=0.2,alpha=0.7) + scale_color_jama(name = 'Study origin')+
    theme_few()+ guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.title = element_text(face = "bold")) +labs(title='PCA: Study origin',
                                                            x=paste0("PC1: ",round(dev.variable.noOutlier[1]*100,1),"%"),
                                                            y=paste0("PC2: ",round(dev.variable.noOutlier[2]*100,1),"%"))

###  Figure 2B-C - PCA 3D plots

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

# Figure 2D, CCS

CCS.genes.data <- exprs(eSet.CCS)

set.seed(1363)
PCA_CCS <- gmodels::fast.prcomp(x = t(CCS.genes.data))
PCA_CCS.data <- as.data.frame(PCA_CCS$x)


PCA.CCS.outlier.idx <- (PCA_CCS.data$PC1 < -50&PCA_CCS.data$PC2 < -50)
PCA.CCS.outlier <- PCA_CCS.data[PCA.CCS.outlier.idx,]

Outlier.CCS.data <- pData(eSet.noTestis)[rownames(PCA.CCS.outlier),c('primary_site_relevel','PCA.tag')]
Outlier.CCS.data$primary_site_relevel <- factor(Outlier.CCS.data$primary_site_relevel)

#table(Outlier.CCS.data$primary_site_relevel,Outlier.CCS.data$PCA.tag)

#                     GTEX-Normal TCGA-Normal TCGA-Tumour
# Adrenal Gland           1           0           0
# Brain                   6           0           0
# Colon                   4           0           0
# Esophagus               4           0           0
# Kidney                  1           0           0
# Lung                    1           0           1
# Ovary                   0           0           1
# Pancreas                2           0           0
# Stomach                 2           0           0
# Thyroid Gland           1           0           0


dev.CCS <- PCA_CCS$sdev^2/sum(PCA_CCS$sdev^2)

PCA.CCS.study.plot <- ggplot(PCA_CCS.data ,aes(x=PC1,y=PC2)) +
    geom_point(aes(color=(eSet.CCS$PCA.tag)),size=0.2,alpha=0.7) + scale_color_jama(name = 'Study origin')+
    geom_point(data=PCA.CCS.outlier,size=0.2,color='red') +

    theme_few()+ guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.title = element_text(face = "bold")) +labs(title='PCA: Study origin**',
                                                            x=paste0("PC1: ",round(dev.CCS[1]*100,1),"%"),
                                                            y=paste0("PC2: ",round(dev.CCS[2]*100,1),"%"))



### PCA 3D plots

PCA.colors.CCS.studies <- pal_jama()(3)[eSet.CCS$PCA.tag]

scatterplot3d::scatterplot3d(PCA_CCS$x[,c(1,3,2)],color = PCA.colors.CCS.studies,pch = 16, grid=TRUE, box=FALSE,type = 'p',angle = -30,cex.symbols = 0.5,cex.axis = 0.5,
                             xlab =paste0("PC1: ",round(dev.CCS[1]*100,1),"%"),
                             zlab =paste0("PC3: ",round(dev.CCS[3]*100,1),"%"),
                             ylab =paste0("PC2: ",round(dev.CCS[2]*100,1),"%"))


######
# removing extra outliers - 24 samples based on PCA - CCS
######

# Figure 2E-1F, CCS

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
# Figure 2G, H, I


Exprs_CCS.noOut <- merge(t(exprs(eSet.CCS.final)),pData(eSet.CCS.final)[c('CCS_ct','primary_site_relevel','PCA.tag')],all.x=T,sort=F,by='row.names')
rownames(Exprs_CCS.noOut) <- Exprs_CCS.noOut$Row.names;Exprs_CCS.noOut <- Exprs_CCS.noOut[,-1]
Exprs_CCS.noOut$CCS_ct <- rescale(Exprs_CCS.noOut$CCS_ct,new.range=c(0,1))

umap.CCS.tissue.plot <- umap(t(Exprs_CCS.noOut[,-c(448:450)]),labels=Exprs_CCS.noOut$primary_site_relevel,controlscale=TRUE,scale=3,legendtitle = 'Tissue type',dotsize=1,axistextsize=10,
                             seed = 1363,legendtextsize=8)+scale_color_manual(values = Tissue.color[-18,]$tissue.colors,name='Tissue type')+
    guides(colour = guide_legend(override.aes = list(size=5),nrow = 2)) +
    theme(legend.title = element_text(face = "bold"))+labs(title='UMAP: Tissue type')

umap.CCS.study.plot <- umap(t(Exprs_CCS.noOut[,-c(448:450)]),labels=Exprs_CCS.noOut$PCA.tag,controlscale=TRUE,scale=3,legendtitle = 'Study origin',dotsize=1,axistextsize=10,
                            seed = 1363,legendtextsize=8)+scale_color_jama(name='Study origin')+
    guides(colour = guide_legend(override.aes = list(size=5),nrow = 2)) +
    theme(legend.title = element_text(face = "bold"))+labs(title='UMAP: Study origin')

umap.CCS.range.plot <- umap(t(Exprs_CCS.noOut[,-c(448:450)]),labels=Exprs_CCS.noOut$CCS_ct,controlscale=TRUE,scale = 1,low = 0,high = 1,legendtitle = 'CCS range',dotsize=1,axistextsize=10,
                            seed = 1363,legendtextsize=8)+scale_color_viridis_c(option="cividis",begin = 0,end = 1)+
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.title = element_text(face = "bold"))+labs(title='UMAP: CCS range')



###############
# Figure 2
# generate figures
###############


plot.legends.Fig2 <- plot_grid(get_legend(umap.CCS.tissue.plot + guides(colour = guide_legend(override.aes = list(size=5),nrow = 7))+
                                              theme(legend.position = 'bottom',legend.title =element_text(size=9,face='bold'),legend.text =element_text(size=7),legend.direction = 'vertical')),
                               get_legend(umap.CCS.study.plot + guides(colour = guide_legend(override.aes = list(size=5),nrow = 7))+
                                              theme(legend.position = 'bottom',legend.title =element_text(size=9,face='bold'),legend.text =element_text(size=7),legend.direction = 'vertical')),
                               get_legend(umap.CCS.range.plot+ guides(colour = guide_legend(override.aes = list(size=5),nrow = 7))+
                                              theme(legend.position = 'bottom',legend.title =element_text(size=9,face='bold'),legend.text =element_text(size=7),legend.direction = 'vertical')),
                               ncol = 3,axis='tblr',rel_widths = c(0.5,0.5,0.5))

#### Figure 2 layout
Figure2 <- plot_grid(
                     PCA.variable.study.plot + theme(legend.position = 'none'), # Fig2A
                     PCA.variable.noOutlier.study.plot + theme(legend.position = 'none'), # Fig2B
                     PCA.variable.noOutlier.tissue.plot + theme(legend.position = 'none'), # Fig2C
                     PCA.CCS.study.plot + theme(legend.position = 'none'), # Fig2D
                     PCA.CCS.final.study.plot + theme(legend.position = 'none'), # Fig2E
                     PCA.CCS.final.tissue.plot + theme(legend.position = 'none'), # Fig2F
                     umap.CCS.study.plot + theme(legend.position = 'none'), # Fig2F
                     umap.CCS.tissue.plot + theme(legend.position = 'none'), # Fig2G
                     umap.CCS.range.plot + theme(legend.position = 'none'), # Fig2H
                     plot.legends.Fig2,                     # legends


                     labels =c('A','B','C','D',
                               'E','F','G',
                               'H','I'),axis = 'tb',ncol = 3,nrow=4,
                     rel_widths = c(1,1,1,1,1),rel_heights = c(1,1,1,0.5))


