load(file="~/Rcodes/Gene.expression_data.RData")
source(file="~/Rcodes/functions.R")

#######
# Figure3A

#######
### add PAM50 subtypes
#######

BRCA.pData$PAM50 <- factor(BRCA.pData$Subtype,levels = c('BRCA_Normal','BRCA_LumA','BRCA_LumB','BRCA_Her2','BRCA_Basal'),labels = c('Normal-Like','Luminal-A','Luminal-B','Her2-Enriched','Basal-like'))

tmp <- merge(pData(eSet.final),BRCA.pData[c('Sample ID','PAM50')],by.x='row.names',by.y='Sample ID',sort=F,all.x=T)
rownames(tmp) <- tmp$Row.names; tmp <- tmp[,-1]
pData(eSet.final) <- tmp[colnames(eSet.final),]

######
set.seed(1363)
CCS_tissue.df <- subset(pData(eSet.final), select = c(primary_site_relevel,PAM50,twoType,twoClass,twoClass.code,tissue.colors,CCS_ct))
CCS_tissue.df$CCS_100 <- rescale(eSet.final$CCS_ct,newrange = c(0,1))
tissue_numbers <- data.frame(type=as.factor(names(table(CCS_tissue.df$twoClass.code))),
                             sum=as.numeric(table(CCS_tissue.df$twoClass.code)))

tissue_numbers$combine <- paste(tissue_numbers$type," (",tissue_numbers$sum,")",sep="")


## relevel boxplots (normals first)

CCS_tissue.df$twoClass.code <- factor(CCS_tissue.df$twoClass.code,levels = c(levels(CCS_tissue.df$twoClass.code)[levels(CCS_tissue.df$twoClass.code) %in% levels(CCS_tissue.df$primary_site_relevel)],levels(CCS_tissue.df$twoClass.code)[levels(CCS_tissue.df$twoClass.code) %out% levels(CCS_tissue.df$primary_site_relevel)]))
highlight_df <- CCS_tissue.df %>%
    filter(twoType=='Normal')

stat.test <- CCS_tissue.df%>%
    group_by(primary_site_relevel) %>%
    rstatix::t_test(CCS_100 ~ twoClass.code,p.adjust.method = 'fdr',conf.level = .95)
stat.test$p.format <- rstatix::p_format(add.p = T,space = T,
                                        stat.test$p.adj, accuracy = 0.001,
                                        leading.zero = T
)


### Figure3A

Figure3A <-  ggplot(CCS_tissue.df, aes(x = twoClass.code, y = CCS_100, fill = primary_site_relevel)) +
    geom_boxplot(outlier.fill = NULL,outlier.size = NULL,show.legend = T,alpha=1,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL,outlier.shape = NULL) +
    geom_jitter(show.legend = F,size=0.8,aes(color=primary_site_relevel,fill=primary_site_relevel)) +
    stat_boxplot(geom='errorbar',width=0.2)+
    geom_boxplot(outlier.fill = NULL,outlier.size = 3,show.legend = F,alpha=1,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL,outlier.shape = NULL)  +
    scale_fill_manual(values = Tissue.color.noTestis$tissue.colors,name='Tissue type') +
    scale_color_manual(values = Tissue.color.noTestis$tissue.colors,name='Tissue type') +
    labs(x='Tissue type',title = 'CCS - Tumor vs normal tissues',y='Cell Cycle Score (CCS)') +
    scale_y_continuous(name = expression(atop(bold("Cell-cycle score (CCS)"))))+
    facet_grid(~primary_site_relevel,scales = 'free_x')+theme_few()+
    scale_x_discrete(breaks=tissue_numbers$type,labels=tissue_numbers$combine)+
    geom_boxplot(data=highlight_df,
                 aes(x=twoClass.code,y=CCS_100),
                 fill='white',alpha=0.9,outlier.fill = NULL,outlier.size = NULL,color='black',outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL,outlier.shape = NULL)+
    theme(legend.key = element_rect(fill='white'),
          strip.text.x =element_text(size = 10,face='bold'),
          axis.title.x = element_text(size = 16,face='bold'),
          axis.title.y = element_text(size = 16,face='bold'),
          plot.title = element_text(size = 17, face = "bold",hjust = c(0.5,0.5)),
          plot.subtitle = element_text(size = 12,hjust = c(0.5,0.5)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.y.right =  element_blank(),
          axis.text.y.left=element_text(size=12,face = 'bold',colour = 'black'),
          axis.text.y.right=element_text(size=14,angle = 90,face = 'bold',colour = c('black','darkgray','gold4'),hjust = c(0,0.4,0.5)),
          axis.title.y.left=element_text(size=16,face = 'bold',colour = 'black'),
          axis.text.x=element_text(size=10,angle = 70,hjust = c(1,1)),axis.ticks.y = element_blank(),
          plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),legend.position = 'top',
          legend.key.size = unit(8, "mm"),legend.title = element_text(size = 10,face='bold')) +  stat_pvalue_manual(stat.test, label = "p.adj.signif",label.size = 3,
                                                                                                                    y.position = 1.05,step.group.by = 'primary_site_relevel',tip.length = 0,step.increase = -0.03)



########

eSet.CCS.final$id_order <- 1:nrow(pData(eSet.CCS.final))
pData(eSet.CCS.final) <-merge(pData(eSet.CCS.final),tumour_purity[c('array','purity')],by.x='row.names',by.y='array',all.x=T,sort=FALSE)
rownames(pData(eSet.CCS.final)) <- eSet.CCS.final$Row.names;pData(eSet.CCS.final) <- pData(eSet.CCS.final)[,-1]
pData(eSet.CCS.final) <- pData(eSet.CCS.final)[order(eSet.CCS.final$id_order),]
eSet.CCS.final$purity <- ifelse(is.na(eSet.CCS.final$purity)==T,1,eSet.CCS.final$purity)

BC_CCS.df <- pData(eSet.CCS.final)[,c('primary_site_relevel','twoClass.code','twoType','tissue.colors','CCS','CCS_ct','purity')]
BC_CCS.df$rownames <- rownames(BC_CCS.df)

### calculate the median value from normal tumours
normal.median <- summaryBy(CCS_ct ~ primary_site_relevel,data = BC_CCS.df[BC_CCS.df$twoType=='Normal',],FUN=c(median))
tumor.median <- summaryBy(CCS_ct ~ primary_site_relevel,data = BC_CCS.df[BC_CCS.df$twoType=='Tumour',],FUN=c(median))

BC_CCS.df <- merge(BC_CCS.df[BC_CCS.df$twoType=='Tumour',],normal.median,by='primary_site_relevel',sort.x=F)

### adjusting for tumour purity

BC_CCS.df$CCS_change <- (BC_CCS.df$CCS_ct-BC_CCS.df$CCS_ct.median)
BC_CCS.df$CCS_change <- BC_CCS.df$CCS_change * BC_CCS.df$purity
BC_CCS.df$twoClass.code <- factor(BC_CCS.df$twoClass.code, levels=new.tumour.color$twoClass.code)

Figure3B <- ggplot(BC_CCS.df, aes(y = rescale(CCS_change,new.range=c(0,1)), x = reorder(twoClass.code, CCS_change, FUN = median), fill = twoClass.code)) +
    geom_boxplot(width=0.5,outlier.fill = NULL,outlier.size = NULL,show.legend = T,alpha=1,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL,outlier.shape = NULL) +
    geom_jitter(width=0.1, show.legend = F,size=0.5,aes(color=twoClass.code,fill=twoClass.code)) +
    stat_boxplot(geom='errorbar',width=0.2)+
    geom_boxplot(width=0.5,outlier.fill = NULL,outlier.size = 3,show.legend = F,alpha=1,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL,outlier.shape = NULL)  +
    geom_half_violin(side = "r",position = position_nudge(x=0.3))+
    geom_segment(x = 1, y = 1.1, xend = 1, yend = 1.03,color='royalblue',linejoin='mitre',size=2,
                 arrow = arrow(length = unit(0.025, "npc"), ends = "last",type = 'closed'))+
    geom_segment(x = 2, y = 1.1, xend = 2, yend = 1.03,color='royalblue',linejoin='mitre',size=2,
                 arrow = arrow(length = unit(0.025, "npc"), ends = "last",type = 'closed'))+
    geom_segment(x = 3, y = 1.1, xend = 3, yend = 1.03,color='royalblue',linejoin='mitre',size=2,
                 arrow = arrow(length = unit(0.025, "npc"), ends = "last",type = 'closed'))+
    geom_segment(x = 5, y = 1.1, xend = 5, yend = 1.03,color='royalblue',linejoin='mitre',size=2,
                 arrow = arrow(length = unit(0.025, "npc"), ends = "last",type = 'closed'))+
    geom_segment(x = 21, y = 1.1, xend = 21, yend = 1.03,color='red',linejoin='mitre',size=2,
                 arrow = arrow(length = unit(0.025, "npc"), ends = "last",type = 'closed'))+
    geom_segment(x = 22, y = 1.1, xend = 22, yend = 1.03,color='red',linejoin='mitre',size=2,
                 arrow = arrow(length = unit(0.025, "npc"), ends = "last",type = 'closed'))+
    geom_segment(x = 23, y = 1.1, xend = 23, yend = 1.03,color='red',linejoin='mitre',size=2,
                 arrow = arrow(length = unit(0.025, "npc"), ends = "last",type = 'closed'))+

    scale_fill_manual(values = new.tumour.color$tissue.colors,name='Cancer types',breaks = new.tumour.color$twoClass.code) +
    scale_color_manual(values = new.tumour.color$tissue.colors,name='Cancer types',breaks = new.tumour.color$twoClass.code)+
    labs(y='BC-CCS',title = 'Relative change of CCS',x='Cancer types') +
    coord_cartesian(ylim=c(0,1.1)) +
    theme_cowplot()+
    theme(legend.key = element_rect(fill='white'),
          strip.text.y =element_text(size = 10,face='bold'),
          axis.title.y = element_text(size = 22,face='bold'),
          axis.title.x = element_text(size = 22,face='bold'),
          plot.title = element_text(size = 24, face = "bold",hjust = c(0.5,0.5)),
          plot.subtitle = element_text(size = 12,hjust = c(0.5,0.5)),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),axis.text.x=element_text(size = 20,angle = 45,hjust = c(1,1)),
          axis.text.y=element_text(size=20),
          legend.key.size = unit(8, "mm"),legend.title = element_text(size = 10,face='bold'))


median.df <- data.frame(Type=as.character(tumor.median$primary_site_relevel),Normal=rescale(normal.median$CCS_ct.median,c(0,1)),
                        Tumour=rescale(tumor.median$CCS_ct.median,c(0,1)))

Figure3C <- ggplot(median.df, aes(x = Normal, y = Tumour, label = Type,fill=Type,color=Type)) +
    stat_cor(size=6,label.x = 1,label.y=1,method = 'spearman',label.sep = '\t',label.y.npc = 'top',p.accuracy = 0.001)+
    geom_vline(xintercept = 0.5,linetype = "dashed") +
    geom_hline(yintercept = 0.5,linetype = "dashed") +
    geom_point(aes(color=Type), size = 7,shape=21) +
    geom_text_repel(size = 5,nudge_y = +0.04,colour='black',fontface='bold',max.overlaps = 2,segment.color=NA)+
    scale_fill_manual(values = Tissue.color$tissue.colors,name='Tissue type',breaks = Tissue.color$tissue, expand = c(0, 0))+
    scale_color_manual(values = rep('black',20),name='Tissue type',breaks = Tissue.color$tissue, expand = c(0, 0)) +

    scale_x_continuous(sec.axis = dup_axis(breaks = c(0.25,0.75),name = "",
                                           labels=c("0.25" = "Low", "0.75" = "High"))) +
    scale_y_continuous(sec.axis = dup_axis(breaks = c(0.25,0.75),name = "",
                                           labels=c("2" = "Low", "8" = "High"))) +
    theme_few(base_size = 22)+
    labs(y='Median CCS - Tumours',x='Median CCS - Normals') +theme(legend.position = 'none',axis.text.y.right = element_text(angle = -90),
                                                                   axis.ticks.y.right = element_blank(),axis.ticks.x.top = element_blank(),
                                                                   axis.title.x = element_text(size = 22,face='bold'),
                                                                   axis.title.y = element_text(size = 22,face='bold'),
                                                                   axis.text= element_text(size = 20,face='bold'))



Figure3.legends <- plot_grid(get_legend(Figure3A+guides(fill = guide_legend(nrow = 2))+
                                            theme(legend.position = 'bottom',legend.title =element_text(size = 20,face='bold'),
                                                  legend.key.size = unit(2,'cm'),legend.text =element_text(size = 15))),nrow = 1,rel_widths = c(1.3,0.7))

Figure3.BC <- plot_grid(Figure3B+theme(legend.position = 'none'),NULL,
                        Figure3C+theme(legend.position = 'none'),ncol=3,align = 'hv',rel_widths = c(1.5,0.01,0.8),labels=c('B','','C'),label_size = 20,scale=c(1,0.8))


#### FIGURE 3.

Figure3 <- plot_grid(Figure3A+theme(legend.position = 'none'),
                     Figure3.BC,
                     Figure3.legends,label_size = 20,labels=c('A','',''),
                     ncol = 1,nrow=3,align = 'h',axis='lr',rel_heights = c(1,1,0.2))

