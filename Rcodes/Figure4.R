# loading required data
######
load(file="~/Box Sync/Karolinska/Projects/CCS_2020/munge/Github/manuscript_data.RData")

### functions
corrFunc <- function(var1,var2,data){
    result=cor.test(data[,var1],data[,var2],method='spearman')
    data.frame(var1, var2, result[c("estimate","p.value","statistic","method")],
               stringsAsFactors=FALSE)
}


#######

# Fig4 - A
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

ER.PR.AR.bp <- ggboxplot(data = melt(G1.G2.bp.data),x ='Group',y='value',facet.by = 'Group',ggtheme = theme_few(),
                         bxp.ERARrorbar = T,bxp.ERARrorbar.width = 0.08,ylab = 'Expression level',size = 0.25,outlier.size=0.7,
                         fill='Group',color='black',palette = c('skyblue1','yellow'))+facet_grid(~variable)+
    stat_compare_means(comparisons = list(c('G1','G2')),label = 'p.signif')+theme(legend.key.size = unit(1, 'cm'))+coord_cartesian(ylim = c(-4,4))

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

module.bp <- ggboxplot(data = melt(module.df),x ='Group',y='value',facet.by = 'Group',ggtheme = theme_few(),
                       bxp.ERARrorbar = T,bxp.ERARrorbar.width = 0.08,ylab = 'Module score',size = 0.25,outlier.size=0.7,
                       fill='Group',color='black',palette = c('skyblue1','yellow'))+
    facet_grid(~variable,labeller = labeller(variable=new_labels))+
    stat_compare_means(comparisons = list(c('G1','G2')),label = 'p.signif')+theme(legend.key.size = unit(1, 'cm'))+coord_cartesian(ylim = c(-4,4))

#######

Fig4.A <- ggarrange(ER.PR.AR.bp+theme(axis.title.x = element_blank()),module.bp,common.legend = T,legend = 'none',ncol = 1)


# Figure 4 - B - volcano plot - DE analyses
# generating volcano plot

#####
### All genes - CCS genes excluded

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


### adjusted for ESR1 and AR

design.tb.ERAR <- model.matrix(~sample_annot$Group+sample_annot$twoClass.code+ER.exprs+AR.exprs)
fit.ERAR <- lmFit(exprs.data,design.tb.ERAR)
fit.ERAR <- eBayes(fit.ERAR,robust=T,trend=T)
limma.toptable.ERAR <- topTable(fit = fit.ERAR,sort.by =  'logFC',coef = 2,adjust.method = 'BH',number = Inf)


Fig4.B.ERAR <- EnhancedVolcano::EnhancedVolcano(limma.toptable.ERAR,
                                                lab = rownames(limma.toptable.ERAR),
                                                x = 'logFC',
                                                y = 'adj.P.Val',drawConnectors = F,labCol = 'black',boxedLabels = F,
                                                title = 'Group 1 vs 2',subtitle ='',ylab = expression(-Log[10]~FDR),
                                                pCutoff = 0.05,gridlines.major = F,gridlines.minor = F,
                                                FCcutoff = 2,colAlpha = 1,axisLabSize = 13,titleLabSize = 11,
                                                caption = bquote(~Log[2]~ "fold change cutoff: 2; FDR cutoff: < 0.05"),
                                                captionLabSize = 7,legendPosition = 'right',
                                                col=c('gray40','skyblue','palegreen4','firebrick3'),
                                                pointSize = 2.0,legendLabels = c('NS',expression(Log[2]~FC),
                                                                                 'FDR', expression(FDR~and~log[2]~FC)),cutoffLineCol = 'black',
                                                labSize = 2.5) +ylim(-20,250)+
    geom_text(aes(x=6, y=-17.5, label="Up-regulated in Group 2"),size=2, color="black") +
    geom_segment(x = 2, y = -12, xend = 10, yend = -12,arrow = arrow(length = unit(0.01, "npc"), ends = 'last',type = 'closed'))+
    geom_text(aes(x=-7, y=-17.5, label="Down-regulated in Group 2"),size=2, color="black") +
    geom_segment(x = -2, y = -12, xend = -12, yend = -12,arrow = arrow(length = unit(0.01, "npc"), ends = "last",type = 'closed'))


# Figure 4D - GSEA - Pathway - MSIGDBR
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
exprs.eSet <- exprs(eSet.final)
## 1501 genes - 8138 samples
exprs.eSet <- t(exprs.eSet[rownames(exprs.eSet) %in% tmp.all.ERAR.cont$GeneID,colnames(exprs.eSet)%in%CCS_relative.df$rownames])

cor.data <- merge(CCS_relative.df[c('CCS_change','rownames')],exprs.eSet,by.y='row.names',by.x='rownames')
rownames(cor.data) <- cor.data$rownames;cor.data <- cor.data[,-1]
vars = data.frame(v1=names(cor.data)[1],v2=names(cor.data)[-1])

corrs = do.call(rbind, mapply(corrFunc, vars[,1], vars[,2], MoreArgs=list(data=cor.data),
                              SIMPLIFY=FALSE))
corrs <- corrs[corrs$p.value <0.05,]
corrs.order <- rbind(corrs[order(corrs$estimate),][1:10,],corrs[order(corrs$estimate,decreasing = T),][1:10,])[c('var2','estimate','p.value')]
corrs.order <- corrs.order[order((corrs.order$estimate)),]
corrs.order$var2 <- factor(corrs.order$var2,levels = corrs.order$var2)
corrs.order$Correlation <- as.factor(ifelse(corrs.order$estimate>0,'Positive','Negative'))

#######

Fig4.C <- ggplot(HM.all.ERAR.cont[order(HM.all.ERAR.cont$padj),][c(1:13,16,24)], aes(reorder(pathway, NES), -NES)) +
    geom_col(aes(fill=padj<0.1)) +
    coord_flip(xlim = c(-1,15) ) +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways",fill='FDR < 0.1') +
    geom_text(aes(y=1.25, x=-0.6, label="Up-regulated in Group 2"),size=2.5,color="black") +
    geom_segment(y = 0.02, x = -0.2, yend = 2.5, xend = -0.2,arrow = arrow(length = unit(0.01, "npc"), ends = 'last',type = 'closed'))+
    geom_text(aes(y=-0.75, x=-0.6, label="Down-regulated in Group 2"),size=2.5, color="black") +
    geom_segment(y= -0.02, x = -0.2, yend = -1.3, xend = -0.2,arrow = arrow(length = unit(0.01, "npc"), ends = "last",type = 'closed'))+
    theme_few() + scale_fill_manual(values = c("slategray3","royalblue3"),breaks = c("FALSE","TRUE"))+theme(axis.text.y=element_text(size = 7,angle = 20,face = 'bold'),axis.text.x=element_text(size=12,face = 'bold'))

Fig4.D <- ggbarplot(corrs.order,x='var2',y='estimate',palette = c('darksalmon','darkorchid4'),ggtheme = theme_few(),
                    fill='Correlation',color='Correlation',xlab = 'Genes',ylab=paste0("Spearman's rank correlation (",'p',")"))+rotate_x_text(angle = 45)+theme(axis.text.x=element_text(size=8.5,face = 'bold'))




#generating heatmap - most mutated genes

#######

### Group1
matrix.adj.G1 <- melt(mutated.number.Group1.adjust[c('cancer',names(sort(colSums(mutated.number.Group1.adjust[,-1]),decreasing = T)))])

### Group2
matrix.adj.G2 <-  melt(mutated.number.Group2.adjust[c('cancer',names(sort(colSums(mutated.number.Group2.adjust[,-1]),decreasing = T)))])

### select genes that have higher mutational load in G1 vs G2

tmp.G1.G2 <- merge(data.frame(G1=colSums(mutated.number.Group1.adjust[,-1])),data.frame(G2=colSums(mutated.number.Group2.adjust[,-1])),by='row.names')
mut.G1.G2.genes <- tmp.G1.G2[tmp.G1.G2$G1 < tmp.G1.G2$G2,]$Row.names


mut.G1.temp <- cbind(mutated.number.Group1[1],mutated.number.Group1[,colnames(mutated.number.Group1)%in%mut.G1.G2.genes],total=as.numeric(table(eSet_CC.Group1.complete$twoClass.code)))
mut.G2.temp <- cbind(mutated.number.Group2[1],mutated.number.Group2[,colnames(mutated.number.Group2)%in%mut.G1.G2.genes],total=as.numeric(table(eSet_CC.Group2.complete$twoClass.code)))

mut.data <- rbind(cbind(cancer=mutated.number.Group1$cancer,round(sweep(mut.G1.temp[-c(1,ncol(mut.G1.temp))],1,mut.G1.temp$total,FUN = "/")*100,digits = 2),Group='Group 1'),
                  cbind(cancer=mutated.number.Group2$cancer,round(sweep(mut.G2.temp[-c(1,ncol(mut.G2.temp))],1,mut.G2.temp$total,FUN = "/")*100,digits = 2),Group='Group 2'))

mut.data.melt <- melt(mut.data)

#reoder columns based on higher mean in Group 2
mut.data.order.idx <- names(sort(colMeans(mut.data[4:7,-c(1,ncol(mut.data))])))
### in case we want to not show 0
mut.data.melt$value[mut.data.melt$value==0] <- NA
mut.data.melt$variable <- factor(mut.data.melt$variable,levels=mut.data.order.idx)

### AMPLIFICATION

tmp.amp.G1.G2 <- merge(data.frame(G1=colSums(Gain.Group1.adjust[,-1])),data.frame(G2=colSums(Gain.Group2.adjust[,-1])),by='row.names')
amp.G1.G2.chr <- tmp.amp.G1.G2[tmp.amp.G1.G2$G1 < tmp.amp.G1.G2$G2,]$Row.names

amp.G1.temp <- cbind(Gain.Group1[1],Gain.Group1[,colnames(Gain.Group1)%in%amp.G1.G2.chr],total=as.numeric(table(tmp.Group1$Type)))
amp.G2.temp <- cbind(Gain.Group2[1],Gain.Group2[,colnames(Gain.Group2)%in%amp.G1.G2.chr],total=as.numeric(table(tmp.Group2$Type)))
amp.data <- rbind(cbind(cancer=Gain.Group1$Type,round(sweep(amp.G1.temp[-c(1,ncol(amp.G1.temp))],1,amp.G1.temp$total,FUN = "/")*100,digits = 2),Group='Group 1'),
                  cbind(cancer=Gain.Group2$Type,round(sweep(amp.G2.temp[-c(1,ncol(amp.G2.temp))],1,amp.G2.temp$total,FUN = "/")*100,digits = 2),Group='Group 2'))

amp.data.melt <- melt(amp.data)


#reoder columns based on higher mean in Group 2
amp.data.order.idx <- names(sort(colMeans(amp.data[4:7,-c(1,ncol(amp.data))])))
### in case we want to not show 0
amp.data.melt$value[amp.data.melt$value==0] <- NA
amp.data.melt$variable <- factor(amp.data.melt$variable,levels=amp.data.order.idx)


### in case we want to not show 0
amp.data.melt$value[amp.data.melt$value==0] <- NA


##### DELETION

tmp.del.G1.G2 <- merge(data.frame(G1=colSums(Loss.Group1.adjust[,-1])),data.frame(G2=colSums(Loss.Group2.adjust[,-1])),by='row.names')
del.G1.G2.chr <- tmp.del.G1.G2[tmp.del.G1.G2$G1 < tmp.del.G1.G2$G2,]$Row.names

del.G1.temp <- cbind(Loss.Group1[1],Loss.Group1[,colnames(Loss.Group1)%in%del.G1.G2.chr],total=as.numeric(table(tmp.Group1$Type)))
del.G2.temp <- cbind(Loss.Group2[1],Loss.Group2[,colnames(Loss.Group2)%in%del.G1.G2.chr],total=as.numeric(table(tmp.Group2$Type)))
del.data <- rbind(cbind(cancer=Loss.Group1$Type,round(sweep(del.G1.temp[-c(1,ncol(del.G1.temp))],1,del.G1.temp$total,FUN = "/")*100,digits = 2),Group='Group 1'),
                  cbind(cancer=Loss.Group2$Type,round(sweep(del.G2.temp[-c(1,ncol(del.G2.temp))],1,del.G2.temp$total,FUN = "/")*100,digits = 2),Group='Group 2'))

del.data.melt <- melt(del.data)


#reoder columns based on higher mean in Group 2
del.data.order.idx <- names(sort(colMeans(del.data[4:7,-c(1,ncol(del.data))])))
del.data.melt$variable <- factor(del.data.melt$variable,levels=del.data.order.idx)

### in case we want to not show 0
del.data.melt$value[del.data.melt$value==0] <- NA

mut.data.melt$round.value <- ifelse(mut.data.melt$value < 1, '< 1%',paste0(round(mut.data.melt$value,digits = 0),'%'))
amp.data.melt$round.value <- ifelse(amp.data.melt$value < 1, '< 1%',paste0(round(amp.data.melt$value,digits = 0),'%'))
del.data.melt$round.value <- ifelse(del.data.melt$value < 1, '< 1%',paste0(round(del.data.melt$value,digits = 0),'%'))

mut.data.HM <- ggplot(mut.data.melt,aes(x=cancer,y=(variable),fill=value))+geom_tile(color='black',width=0.9, height=0.8,size=0.1)+
    geom_text(aes(label=round.value),color='black',size=4) +
    labs(x='',y='')+scale_fill_gradient(low="white", high="skyblue4",guide = 'colorbar',na.value = NA,limits=c(0,75),breaks=c(0,25,50,75)
                                        ,name='Adj-mutation load (%)',labels = function(x) paste0(x, '%'))+
    facet_grid(~Group,scales = 'free_x',space='free',shrink = T)+
    theme_few()+theme(legend.text = element_text(size = 9,hjust = 1),
                      legend.position = 'right',strip.text = element_text(size=12,face = 'bold'),
                      axis.text.x  = element_text(size=10,face = 'bold'),axis.text.y  = element_text(size=10,face = 'bold'))


Gain.data.HM <- ggplot(amp.data.melt,aes(x=cancer,y=(variable),fill=value))+geom_tile(color='black',width=0.9, height=0.8,size=0.1)+
    geom_text(aes(label=round.value),color='black',size=4) +
    labs(x='',y='')+scale_fill_gradient(low="white", high="firebrick3",guide = 'colorbar',na.value = NA,limits=c(0,75),breaks=c(0,25,50,75)
                                        ,name='Adj-Amplification (%)',labels = function(x) paste0(x, '%'))+
    facet_grid(~Group,scales = 'free_x',space='free',shrink = T)+
    theme_few()+theme(legend.text = element_text(size = 9,hjust = 1),
                      legend.position = 'right',strip.text = element_text(size=12,face = 'bold'),
                      axis.text.x  = element_text(size=10,face = 'bold'),axis.text.y  = element_text(size=10,face = 'bold'))


Loss.data.HM <- ggplot(del.data.melt,aes(x=cancer,y=(variable),fill=value))+geom_tile(color='black',width=0.9, height=0.8,size=0.1)+
    geom_text(aes(label=round.value),color='black',size=4) +
    labs(x='',y='')+scale_fill_gradient(low="white", high="gold",guide = 'colorbar',na.value = NA,limits=c(0,75),breaks=c(0,25,50,75)
                                        ,name='Adj-Deletion (%)',labels = function(x) paste0(x, '%'))+
    facet_grid(~Group,scales = 'free_x',space='free',shrink = T)+
    theme_few()+theme(legend.text = element_text(size = 9,hjust = 1),
                      legend.position = 'right',strip.text = element_text(size=12,face = 'bold'),
                      axis.text.x  = element_text(size=10,face = 'bold'),axis.text.y  = element_text(size=10,face = 'bold'))


#######



# Figure 4.F
#### selecting top 10 mutated genes - FDR < 0.05 , Fisher exact

######

### select unique type of mutation for each gene
tmp.G1.mut.type <- unique(as.data.frame(CC.Group1_maf@data[,c('Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification')]))
tmp.G1.mut.type <- merge(tmp.G1.mut.type,pData(eSet_CC.Group1.complete)['twoClass.code'],by.x='Tumor_Sample_Barcode',by.y='row.names')
tmp.G1.mut.type$cancer.sample <- as.numeric(as.character(factor(tmp.G1.mut.type$twoClass.code,labels = as.numeric(table(eSet_CC.Group1.complete$twoClass.code)))))
tmp.G1.mut.type$mut.number <- (1/tmp.G1.mut.type$cancer.sample)*dim(eSet_CC.Group1.complete)[[2]]
tmp.G1.mut.type <- tmp.G1.mut.type[order(tmp.G1.mut.type$Hugo_Symbol),]
tmp.G1.mut.type$mut.counts <- 1

G1.mut.type <- aggregate(tmp.G1.mut.type$mut.number,by=list(Gene=tmp.G1.mut.type$Hugo_Symbol,Variant_Classification=tmp.G1.mut.type$Variant_Classification),FUN=sum)


tmp.G2.mut.type <- unique(as.data.frame(CC.Group2_maf@data[,c('Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification')]))
tmp.G2.mut.type <- merge(tmp.G2.mut.type,pData(eSet_CC.Group2.complete)['twoClass.code'],by.x='Tumor_Sample_Barcode',by.y='row.names')
tmp.G2.mut.type$cancer.sample <- as.numeric(as.character(factor(tmp.G2.mut.type$twoClass.code,labels = as.numeric(table(eSet_CC.Group2.complete$twoClass.code)))))
tmp.G2.mut.type$mut.number <- (1/tmp.G2.mut.type$cancer.sample)*dim(eSet_CC.Group2.complete)[[2]]
tmp.G2.mut.type <- tmp.G2.mut.type[order(tmp.G2.mut.type$Hugo_Symbol),]
tmp.G2.mut.type$mut.counts <- 1

G2.mut.type <- aggregate(tmp.G2.mut.type$mut.number,by=list(Gene=tmp.G2.mut.type$Hugo_Symbol,Variant_Classification=tmp.G2.mut.type$Variant_Classification),FUN=sum)

mut.type.counts.G1.G2 <- merge(G1.mut.type,G2.mut.type,by=c('Gene','Variant_Classification'),all=T)
# exclude genes with NA in either of groups
mut.type.counts.G1.G2 <- na.omit(mut.type.counts.G1.G2);colnames(mut.type.counts.G1.G2)[3:4] <- c('Group1','Group2')

### subset genes that have more mutations in G1 vs G2
mut.more.G1.G2 <- mut.type.counts.G1.G2[mut.type.counts.G1.G2$Group1 < mut.type.counts.G1.G2$Group2,]
mut.more.G1.G2 <- cbind(mut.more.G1.G2, t(apply(mut.more.G1.G2[3:4], 1, function(x) {
    ch <- chisq.test(x)
    c(unname(ch$statistic), ch$p.value)})))
colnames(mut.more.G1.G2)[5:6] <- c('Chi.stat','p')
mut.more.G1.G2$signif <- mut.more.G1.G2$p < 0.05
mut.more.G1.G2.list <- split(mut.more.G1.G2,mut.more.G1.G2$Variant_Classification)


genes.sig.mut.type <- unique(mut.more.G1.G2[mut.more.G1.G2$signif==T,]$Gene)
#[1] "BRCA2" "FBXW7" "KLF5"  "NIPBL" "PTCH1"

# mut.G1.G2.genes
# [1] "BCL2L11" "CEBPA"   "FBXW7"   "FOXQ1"   "KLF5"    "MYC"     "PTMA"    "RHOB"    "ZBTB7B"

Group1.mut.top <- tmp.G1.mut.type[tmp.G1.mut.type$Hugo_Symbol%in%genes.sig.mut.type,-c(1,5)]
Group1.mut.top <- aggregate(Group1.mut.top$mut.counts,by=list(Gene=Group1.mut.top$Hugo_Symbol,Variant_Classification=factor(Group1.mut.top$Variant_Classification),twoClass.code=Group1.mut.top$twoClass.code),FUN=sum)
Group1.mut.top$total <- as.numeric(as.character(factor(Group1.mut.top$twoClass.code,labels = as.numeric(table(eSet_CC.Group1.complete$twoClass.code)))))

# table(factor(Group1.mut.top$Variant_Classification))
# Frame_Shift_Del   Frame_Shift_Ins      In_Frame_Del      In_Frame_Ins Missense_Mutation Nonsense_Mutation       Splice_Site
# 9                 7                 2                 1                14                 8                 4
Group1.mut.top$mut.type <- paste0(Group1.mut.top$Gene,":",factor(Group1.mut.top$Variant_Classification,labels =c('Frame shift deletion','Frame shift insertion','In frame deletion',
                                                                                                                 'In frame insertion','Missense mutation','Nonsense mutation','Splice site')))


Group2.mut.top <- tmp.G2.mut.type[tmp.G2.mut.type$Hugo_Symbol%in%genes.sig.mut.type,-c(1,5)]
Group2.mut.top <- aggregate(Group2.mut.top$mut.counts,by=list(Gene=Group2.mut.top$Hugo_Symbol,Variant_Classification=factor(Group2.mut.top$Variant_Classification),twoClass.code=Group2.mut.top$twoClass.code),FUN=sum)
Group2.mut.top$total <- as.numeric(as.character(factor(Group2.mut.top$twoClass.code,labels = as.numeric(table(eSet_CC.Group2.complete$twoClass.code)))))

# table(factor(Group2.mut.top$Variant_Classification))
# Frame_Shift_Del   Frame_Shift_Ins      In_Frame_Del Missense_Mutation Nonsense_Mutation       Splice_Site
# 5                 2                 1                14                 9                 2
Group2.mut.top$mut.type <- paste0(Group2.mut.top$Gene,":",factor(Group2.mut.top$Variant_Classification,labels =c('Frame shift deletion','Frame shift insertion','In frame deletion',
                                                                                                                 'Missense mutation','Nonsense mutation','Splice site')))

mut.type.data <- rbind.data.frame(cbind(mutation.type=Group1.mut.top$mut.type,cancer=as.character(Group1.mut.top$twoClass.code),round((Group1.mut.top$x/Group1.mut.top$total)*100,digits = 2),Group='Group 1'),
                                  cbind(mutation.type=Group2.mut.top$mut.type,cancer=as.character(Group2.mut.top$twoClass.code),round((Group2.mut.top$x/Group2.mut.top$total)*100,digits = 2),Group='Group 2'))

mut.type.data$V3 <- as.numeric(mut.type.data$V3)

### in case we want to not show 0
mut.type.data$V3[mut.type.data$V3==0] <- NA
mut.type.data$round.value <- ifelse(mut.type.data$V3 < 1, '< 1%',paste0(round(mut.type.data$V3,digits = 0),'%'))
mut.type.data$mutation.type <- factor(mut.type.data$mutation.type,rev(names(table(mut.type.data$mutation.type))))
######
mut.type.HM <- ggplot(mut.type.data,aes(x=cancer,y=(mutation.type),fill=V3))+geom_tile(color='black',width=0.9, height=0.8,size=0.1)+
    geom_text(aes(label=round.value),color='black',size=3) +
    labs(x='',y='')+scale_fill_gradient(low="white", high="darkturquoise",guide = 'colorbar',na.value = NA,limits=c(0,45),breaks=c(0,15,30,45)
                                        ,name='Adj-mutational type (%)',labels = function(x) paste0(x, '%'))+
    facet_grid(~Group,scales = 'free_x',space='free',shrink = T)+
    theme_few()+theme(legend.text = element_text(size = 9,hjust = 1),
                      legend.position = 'right',strip.text = element_text(size=12,face = 'bold'),
                      axis.text.x  = element_text(size=10,face = 'bold'),axis.text.y  = element_text(size=8,face = 'bold'))



###### COSMIC signature plot - heatmap
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library("NMF")

Group1.tnm = trinucleotideMatrix(maf = CC.Group1_maf, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
Group2.tnm = trinucleotideMatrix(maf = CC.Group2_maf, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

Group1.sig = extractSignatures(mat = Group1.tnm,n = 4)
Group2.sig = extractSignatures(mat = Group2.tnm, n = 4)

Group1.sig.cosm = compareSignatures(nmfRes = Group1.sig, sig_db = "legacy")
Group2.sig.cosm = compareSignatures(nmfRes = Group2.sig, sig_db = "legacy")

Group1.heatmap <- Heatmap(mat=Group1.sig.cosm$cosine_similarities,cluster_rows = F,row_title="Group 1",row_labels = c(paste('Sig',1:4)),name = 'Cosine Similarity',column_names_gp=gpar(fontsize=7),column_names_rot = 90,heatmap_legend_param = list(direction = "horizontal"),column_labels = paste0(1:30," (",Group1.sig.cosm$aetiology_db$aetiology,")"))
Group2.heatmap <- Heatmap(mat=Group2.sig.cosm$cosine_similarities,cluster_rows = F,row_title="Group 2",row_labels = c(paste('Sig',1:4)),name = 'Cosine Similarity',column_names_gp=gpar(fontsize=7),column_names_rot = 90,column_labels = paste0(1:30," (",Group2.sig.cosm$aetiology_db$aetiology,")"))
ht_list <- Group1.heatmap + Group2.heatmap
COSMIC_heatmap <- draw(ht_list,row_title_side='right',heatmap_legend_side="bottom",
                       #column_title = "Cosine similarity against validated signatures",
                       column_title_gp = gpar(fontsize = 14,type='bold'))


##### left space for H panel - added manually to the figure
Fig4.EH <- plot_grid(Gain.data.HM+theme(legend.position = 'none'),
                     Loss.data.HM+theme(legend.position = 'none'),
                     mut.data.HM+theme(legend.position = 'none'),
                     NULL,NULL,
                     theme(legend.position = 'none'),ncol = 5,nrow=1,align = 'v',rel_widths = c(1,1,1,0.2,1.5),labels=c('E','F','G','','H'))


Fig4.leg <-  plot_grid(get_legend(ER.PR.AR.bp+theme(legend.position = 'bottom')),get_legend(Fig4.B.ERAR),get_legend(Fig4.C),get_legend(Fig4.D),
                       get_legend(Gain.data.HM),get_legend(Loss.data.HM),get_legend(mut.data.HM),align = 'v',nrow=1,rel_widths = c(0.1,0.1,0.05,0.05,
                                                                                                                                   0.05,0.05,0.05),rel_heights = c(0.2))

Fig4CD <- plot_grid(Fig4.C+theme(legend.position = 'none'),Fig4.D+theme(legend.position = 'none'),align = 'h',labels=c('C','D'),rel_widths = c(1,0.8))

Fig4.AD <- plot_grid(Fig4.A,Fig4.B.ERAR+theme(legend.position = 'none'),Fig4CD,
                     rel_widths = c(0.75,0.75,1.5),scale = c(0.9,0.9,1),labels=c('A','B'),ncol=3,axis='tb')

Fig4.AH <- plot_grid(Fig4.AD,NULL,Fig4.EH,nrow=3,rel_heights = c(1,0.1,1.2))


#### Final Figure 4

Fig4 <- plot_grid(Fig4.AH,
                  Fig4.leg,align = 'hv',
                  byrow = F,nrow = 2,ncol=1,rel_heights = c(1.7,0.3),scale=c(1,0.8))

