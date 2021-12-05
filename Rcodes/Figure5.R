## loading required data
######
load(file="~/Box Sync/Karolinska/Projects/CCS_2020/munge/Github/Gene.expression_data.RData")
load(file="~/Box Sync/Karolinska/Projects/CCS_2020/munge/Github/mutations_data.RData")
source(file="~/Box Sync/Karolinska/Projects/CCS_2020/munge/Github/functions.R")

#######

# Figure5 - A
# Boxplots
#######

MDSIG <- msigdbr(species = 'Homo sapiens',category = 'C5')
ER_module <- subset(MDSIG,subset=MDSIG$gs_name=='GOMF_ESTROGEN_RECEPTOR_ACTIVITY')$gene_symbol
PR_module <- subset(MDSIG,subset=MDSIG$gs_name=='GOBP_PROGESTERONE_RECEPTOR_SIGNALING_PATHWAY')$gene_symbol
AR_module <- subset(MDSIG,subset=MDSIG$gs_name=='GOBP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY')$gene_symbol

### boxplot ER PR and AR comparison
G1.exprs <- exprs(eSet_CC.Group1.complete)
G2.exprs <- exprs(eSet_CC.Group2.complete)
G1.ERAR.PR.AR <- data.frame(t(G1.exprs[c('ESR1','PGR','AR'),]),Group='G1')
G2.ERAR.PR.AR <- data.frame(t(G2.exprs[c('ESR1','PGR','AR'),]),Group='G2')
G1.G2.bp.data <- rbind(G1.ERAR.PR.AR,G2.ERAR.PR.AR)
G1.G2.bp.data <- data.frame(scale(G1.G2.bp.data[,-4]),Group=G1.G2.bp.data$Group)

#### gene modules
module.df.G1 <- data.frame(ER.corr=rowSums(t(G1.exprs[rownames(G1.exprs)%in%ER_module,])),
                           PR.corr=rowSums(t(G1.exprs[rownames(G1.exprs)%in%PR_module,])),
                           AR.corr=rowSums(t(G1.exprs[rownames(G1.exprs)%in%AR_module,])),Group='G1')
module.df.G2 <- data.frame(ER.corr=rowSums(t(G2.exprs[rownames(G2.exprs)%in%ER_module,])),
                           PR.corr=rowSums(t(G2.exprs[rownames(G2.exprs)%in%PR_module,])),
                           AR.corr=rowSums(t(G2.exprs[rownames(G2.exprs)%in%AR_module,])),Group='G2')

module.df <- rbind(module.df.G1,module.df.G2)
module.df <- data.frame(scale(module.df[,-4]),Group=module.df$Group)

new_labels <- c("ER.corr" = "ER module", "PR.corr" = "PR module", "AR.corr" = "AR module")

ER.PR.AR.bp <- ggboxplot(data = melt(G1.G2.bp.data),x ='Group',y='value',facet.by = 'Group',ggtheme = theme_few(),
                         bxp.ERARrorbar = T,bxp.ERARrorbar.width = 0.08,ylab = 'Expression level',size = 0.25,outlier.size=0.7,
                         fill='Group',color='black',palette = c('skyblue1','yellow'))+facet_grid(~variable)+
    stat_compare_means(comparisons = list(c('G1','G2')),label = 'p.signif')+theme(legend.key.size = unit(1, 'cm'))+coord_cartesian(ylim = c(-4,4))


module.bp <- ggboxplot(data = melt(module.df),x ='Group',y='value',facet.by = 'Group',ggtheme = theme_few(),
                       bxp.ERARrorbar = T,bxp.ERARrorbar.width = 0.08,ylab = 'Module score',size = 0.25,outlier.size=0.7,
                       fill='Group',color='black',palette = c('skyblue1','yellow'))+
    facet_grid(~variable,labeller = labeller(variable=new_labels))+
    stat_compare_means(comparisons = list(c('G1','G2')),label = 'p.signif')+theme(legend.key.size = unit(1, 'cm'))+coord_cartesian(ylim = c(-4,4))

#######

Figure5A <- ggarrange(ER.PR.AR.bp+theme(axis.title.x = element_blank()),module.bp,common.legend = T,legend = 'none',ncol = 1)


###########


eSet.all.data <- eSet.final[rownames(eSet.final)%out%cc.Genes$hgnc_symbol,colnames(eSet.final)%in%c(colnames(eSet_CC.Group1.complete),colnames(eSet_CC.Group2.complete))]
eSet.all.data$Group <- factor(ifelse(eSet.all.data$twoClass.code%in%c("CESC","OV","UCS"),'Group2','Group1'),levels=c('Group2','Group1'))
eSet.all.data$twoClass.code <- factor(eSet.all.data$twoClass.code)

exprs.data <- exprs(eSet.all.data)
sample_annot <- pData(eSet.all.data)[c('twoClass.code','Group')]
design.tb <- model.matrix(~sample_annot$Group+sample_annot$twoClass.code)
fit <- lmFit(exprs.data,design.tb)
fit <- eBayes(fit,robust=T,trend=T)
limma.toptable <- topTable(fit = fit,sort.by =  'logFC',coef = 2,adjust.method = 'BH',number = Inf)
ER.exprs <- exprs.data[rownames(exprs.data)=='ESR1',]
AR.exprs <- exprs.data[rownames(exprs.data)=='AR',]

### adjusted for ER continuous = ER.exprs

design.tb.ERAR <- model.matrix(~sample_annot$Group+sample_annot$twoClass.code+ER.exprs+AR.exprs)
fit.ERAR <- lmFit(exprs.data,design.tb.ERAR)
fit.ERAR <- eBayes(fit.ERAR,robust=T,trend=T)
limma.toptable.ERAR <- topTable(fit = fit.ERAR,sort.by =  'logFC',coef = 2,adjust.method = 'BH',number = Inf)

#######
# Figure5D - GSEA - Pathway - MSIGDBR
#######

Hallmarks_db <- msigdbr(species = "Homo sapiens",category = c('H'))
Hallmarks_list = split(x = Hallmarks_db$gene_symbol, f = Hallmarks_db$gs_name)

### all genes + ER continuous

tmp.all.ERAR.cont <- subset(limma.toptable.ERAR,subset=(limma.toptable.ERAR$adj.P.Val < 0.05&abs(limma.toptable.ERAR$logFC)>=2))
tmp.all.ERAR.cont <- data.frame(GeneID=rownames(tmp.all.ERAR.cont),tmp.all.ERAR.cont)
tmp.all.ERAR.cont$fcSign <- sign(tmp.all.ERAR.cont$logFC)
tmp.all.ERAR.cont$logP <- -log10(tmp.all.ERAR.cont$P.Value)
tmp.all.ERAR.cont$gRank <- tmp.all.ERAR.cont$logP/tmp.all.ERAR.cont$fcSign

### handle the tie ranks
set.seed(100)
tmp.all.ERAR.cont <- tmp.all.ERAR.cont %>% arrange(desc(gRank))
genes.all.ERAR.cont <- tmp.all.ERAR.cont$gRank
names(genes.all.ERAR.cont) <- tmp.all.ERAR.cont$GeneID
genes.all.ERAR.cont <- genes.all.ERAR.cont[!is.infinite(genes.all.ERAR.cont)]

HM.all.ERAR.cont <- fgsea::fgseaMultilevel(pathways=Hallmarks_list,genes.all.ERAR.cont,maxSize =500,eps = 0,minSize = 5)

# compare the correlation to CCS-relative score

#######
exprs.eSet.2group <- exprs(eSet.final[,eSet.final$twoClass.code%in%c('HNSC','KIRP','KICH','UCEC','CESC','OV','UCS')])
# assayData: 25453 features, 1831 samples
exprs.eSet.2group <- t(exprs.eSet.2group[rownames(exprs.eSet.2group) %in% tmp.all.ERAR.cont$GeneID,colnames(exprs.eSet.2group)%in%CCS_relative.df$rownames])

cor.data.2group <- merge(CCS_relative.df[c('CCS_change','rownames')],exprs.eSet.2group,by.y='row.names',by.x='rownames')
rownames(cor.data.2group) <- cor.data.2group$rownames;cor.data.2group <- cor.data.2group[,-1]
vars = data.frame(v1=names(cor.data.2group)[1],v2=names(cor.data.2group)[-1])

corrs.2group = do.call(rbind, mapply(corrFunc, vars[,1], vars[,2], MoreArgs=list(data=cor.data.2group),
                                     SIMPLIFY=FALSE))
corrs.2group <- corrs.2group[corrs.2group$p.value <0.05,]
corrs.2group.order <- rbind(corrs.2group[order(corrs.2group$estimate),][1:10,],corrs.2group[order(corrs.2group$estimate,decreasing = T),][1:10,])[c('var2','estimate','p.value')]
corrs.2group.order <- corrs.2group.order[order((corrs.2group.order$estimate)),]
corrs.2group.order$var2 <- factor(corrs.2group.order$var2,levels = corrs.2group.order$var2)
corrs.2group.order$Correlation <- as.factor(ifelse(corrs.2group.order$estimate>0,'Positive','Negative'))

### other cancers

exprs.eSet <- exprs(eSet.final[,eSet.final$twoClass.code%out%c('HNSC','KIRP','KICH','UCEC','CESC','OV','UCS')&eSet.final$twoType=='Tumour'])
#assayData: 25453 features, 6307 samples

exprs.eSet <- t(exprs.eSet[rownames(exprs.eSet) %in% as.character(corrs.2group.order$var2),colnames(exprs.eSet)%in%CCS_relative.df$rownames])

cor.data <- merge(CCS_relative.df[c('CCS_change','rownames')],exprs.eSet,by.y='row.names',by.x='rownames')
rownames(cor.data) <- cor.data$rownames;cor.data <- cor.data[,-1]
vars = data.frame(v1=names(cor.data)[1],v2=names(cor.data)[-1])

corrs = do.call(rbind, mapply(corrFunc, vars[,1], vars[,2], MoreArgs=list(data=cor.data),
                              SIMPLIFY=FALSE))

corrs.order <- rbind(corrs[order(corrs$estimate),][1:10,],corrs[order(corrs$estimate,decreasing = T),][1:10,])[c('var2','estimate','p.value')]
corrs.order <- corrs.order[order((corrs.order$estimate)),]
corrs.order$var2 <- factor(corrs.order$var2,levels = corrs.order$var2)
corrs.order$Correlation <- as.factor(ifelse(corrs.order$estimate>0,'Positive','Negative'))

tmp.pathway.ERAR <- HM.all.ERAR.cont
tmp.pathway.ERAR$pathway <- gsub("HALLMARK_","",tmp.pathway.ERAR$pathway)
tmp.pathway.ERAR$pathway <- gsub("_"," ",tmp.pathway.ERAR$pathway)
tmp.pathway.ERAR$pathway <- str_to_title(tmp.pathway.ERAR$pathway)

#######
# labels on the volcano plot

gene_labels <- unique(c("C19orf57","TCF19","ADRB2","NEFL","DNASE1L3",as.character(corrs.2group.order$var2),
                        rownames(limma.toptable.ERAR[limma.toptable.ERAR$logFC>3&-log10(limma.toptable.ERAR$adj.P.Val)>10,]),
                        rownames(limma.toptable.ERAR[limma.toptable.ERAR$logFC< -6.5&-log10(limma.toptable.ERAR$adj.P.Val)>100,])))
options(ggrepel.max.overlaps = Inf)

# Figure5B - volcano plot - DE analyses

Figure5B <- EnhancedVolcano::EnhancedVolcano(limma.toptable.ERAR,
                                                lab = rownames(limma.toptable.ERAR),
                                                x = 'logFC',selectLab = gene_labels[c(1:10,40:length(gene_labels))],labFace = 'bold',
                                                y = 'adj.P.Val',drawConnectors =T,labCol = 'black',boxedLabels = F,
                                                title = 'Group 1 vs 2',subtitle ='',ylab = expression(-Log[10]~FDR),colConnectors = 'grey20',
                                                pCutoff = 0.05,gridlines.major = F,gridlines.minor = F,widthConnectors = 0.1,
                                                FCcutoff = 2,colAlpha = 0.6,axisLabSize = 13,titleLabSize = 11,
                                                captionLabSize = 0,legendPosition = 'bottom',caption = '',
                                                col=c('gray40','skyblue','palegreen4','firebrick3'),
                                                pointSize = 2.0,legendLabels = c('NS',expression(Log[2]~FC),
                                                                                 'FDR', expression(FDR~and~log[2]~FC)),cutoffLineCol = 'black',
                                                labSize = 3)+ylim(0,250)+
    geom_text(aes(x=6, y=240, label="Up-regulated in\nGroup 2"),size=4, color="black") +
    geom_segment(x = 2, y = 210, xend = 10, yend =210,arrow = arrow(length = unit(0.01, "npc"), ends = 'last',type = 'closed'))+
    geom_text(aes(x=-7, y=240, label="Down-regulated in\nGroup 2"),size=4, color="black") +
    geom_segment(x = -2, y = 210, xend = -12, yend = 210,arrow = arrow(length = unit(0.01, "npc"), ends = "last",type = 'closed'))+
    theme(plot.title = element_text(size = 18,hjust = c(0.5,0.5)))


### Pathway enrichment analyses

Figure5C <- ggplot(tmp.pathway.ERAR[order(tmp.pathway.ERAR$padj),][c(1:13,16,24)], aes(reorder(pathway, NES), -NES)) +
    geom_col(aes(fill=padj<0.1)) +
    coord_flip(xlim = c(-1,15) ) +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways",fill='FDR < 0.1') +
    geom_text(aes(y=1.25, x=-0.8, label="Up-regulated in\nGroup 2"),size=4,color="black") +
    geom_segment(y = 0.02, x = -0.1, yend = 2.5, xend = -0.1,arrow = arrow(length = unit(0.01, "npc"), ends = 'last',type = 'closed'))+
    geom_text(aes(y=-0.75, x=-0.8, label="Down-regulated in\nGroup 2"),size=4, color="black") +
    geom_segment(y= -0.02, x = -0.1, yend = -1.35, xend = -0.1,arrow = arrow(length = unit(0.01, "npc"), ends = "last",type = 'closed'))+
    theme_few() + scale_fill_manual(values = c("slategray3","royalblue3"),breaks = c("FALSE","TRUE"),name="",labels=c("FDR > 0.1","FDR < 0.1"))+theme(legend.position=c(.85,.91),
                                                                                                                                                      axis.text.y=element_text(size = 10,angle = 30,face = 'bold'),
                                                                                                                                                      axis.text.x=element_text(size=12,face = 'bold'))

#### Correlation plots

corrs.2group.order$sig <- corrs.2group.order$p.value < 0.05
colnames(corrs.2group.order)[5] <- 'p < 0.05'


Figure5D <- ggbarplot(corrs.2group.order,x='var2',y='estimate',palette = c('darksalmon','darkorchid4'),ggtheme = theme_few(),title = '2 Groups (n = 1831)',
                    fill='Correlation',xlab = 'Genes',ylab=paste0("Spearman's rank correlation (",'p',")"))+rotate_x_text(angle = 45)+theme(legend.position=c(.2,.9),axis.text.x=element_text(size=8.5,face = 'bold'))+
    scale_color_manual(values = c('black','red'),breaks = c('TRUE','FALSE'))+coord_cartesian(ylim = c(-0.6,0.6))

corrs.order$sig <- corrs.order$p.value < 0.05
colnames(corrs.order)[5] <- 'p < 0.05'

Figure5E <- ggbarplot(corrs.order,x='var2',y='estimate',palette = c('darksalmon','darkorchid4'),ggtheme = theme_few(),title = 'Other cancers (n = 6307)',
                    fill='Correlation',xlab = 'Genes',ylab='')+rotate_x_text(angle = 45)+theme(legend.position=c(.2,.9),axis.text.x=element_text(size=8.5,face = 'bold'))+
    scale_color_manual(values = c('black','red'),breaks = c('TRUE','FALSE'))+ geom_hline(yintercept=0.3, linetype="dashed", color = "black")+coord_cartesian(ylim = c(-0.6,0.6))



##### Final figure 5

Figure5 <- plot_grid(Figure5A+theme(legend.position = 'none'),
                  Figure5B,
                  Figure5C,
                  Figure5D,
                  Figure5E,ncol = 3,nrow=2,labels=c('A','B','C','D','E',''))

