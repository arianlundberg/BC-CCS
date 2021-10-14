load(file="~/Box Sync/Karolinska/Projects/CCS_2020/munge/Github/manuscript_data.RData")

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
          legend.key.size = unit(8, "mm"),legend.title = element_text(size = 10,face='bold')) +  stat_pvalue_manual(stat.test, label = "p.format",label.size = 3,
                                                                                                                    y.position = 1.05,step.group.by = 'primary_site_relevel',tip.length = 0,step.increase = -0.03)



########

eSet.CCS.final$id_order <- 1:nrow(pData(eSet.CCS.final))
pData(eSet.CCS.final) <-merge(pData(eSet.CCS.final),tumour_purity[c('array','purity')],by.x='row.names',by.y='array',all.x=T,sort=FALSE)
rownames(pData(eSet.CCS.final)) <- eSet.CCS.final$Row.names;pData(eSet.CCS.final) <- pData(eSet.CCS.final)[,-1]
pData(eSet.CCS.final) <- pData(eSet.CCS.final)[order(eSet.CCS.final$id_order),]
eSet.CCS.final$purity <- ifelse(is.na(eSet.CCS.final$purity)==T,1,eSet.CCS.final$purity)

CCS_relative.df <- pData(eSet.CCS.final)[,c('primary_site_relevel','twoClass.code','twoType','tissue.colors','CCS','CCS_ct','purity')]
CCS_relative.df$rownames <- rownames(CCS_relative.df)

### calculate the median value from normal tumours
normal.median <- summaryBy(CCS_ct ~ primary_site_relevel,data = CCS_relative.df[CCS_relative.df$twoType=='Normal',],FUN=c(median))
CCS_relative.df <- merge(CCS_relative.df[CCS_relative.df$twoType=='Tumour',],normal.median,by='primary_site_relevel',sort.x=F)
CCS_relative.df$CCS_change <- (CCS_relative.df$CCS_ct-CCS_relative.df$CCS_ct.median)

### adjust for tumour purity

CCS_relative.df$CCS_change <- CCS_relative.df$CCS_change * CCS_relative.df$purity

new.tumour.color <- merge(Tissue.color,unique(pData(eSet.CCS.final)[eSet.CCS.final$twoType=='Tumour',c('twoClass.code','tissue.colors')]),by='tissue.colors')
new.tumour.color <- new.tumour.color[order(new.tumour.color$tissue),c('tissue','twoClass.code','tissue.colors')]
new.tumour.color$twoClass.code <- factor(new.tumour.color$twoClass.code)

CCS_relative.df$twoClass.code <- factor(CCS_relative.df$twoClass.code, levels=new.tumour.color$twoClass.code)

anatomy.male <- merge(Tissue.color,hgMale_key[hgMale_key$type %in% c('digestion','respiratory','reproductive'),],by='organ',sort=F)
anatomy.male <- anatomy.male[anatomy.male$organ %out% c('testis','kidney'),]
anatomy.male <- anatomy.male[,c("organ","type","tissue.colors","value")];colnames(anatomy.male)[3] <- 'colour'

anatomy.female <- merge(Tissue.color,rbind(hgFemale_key[hgFemale_key$type %out% c('digestion','respiratory'),],hgFemale_key[hgFemale_key$organ=='kidney',]),by='organ',sort=F)
anatomy.female <- anatomy.female[c("organ","type","tissue.colors","value")];colnames(anatomy.female)[3] <- 'colour'

CCS_change.median <- summaryBy(CCS_change~primary_site_relevel,CCS_relative.df,FUN=c(median))
CCS_change.median <- CCS_change.median[order(CCS_change.median$primary_site_relevel),]
CCS_change.median$organ <- Tissue.color.noTestis$organ
anatomy.female <- merge(anatomy.female,CCS_change.median[c('organ','CCS_change.median')],by='organ')
anatomy.male <- merge(anatomy.male,CCS_change.median[c('organ','CCS_change.median')],by='organ')
anatomy.male$value <- rescale(log2(anatomy.male$CCS_change.median),new.range=c(0,1))
anatomy.female$value <- rescale(log2(anatomy.female$CCS_change.median),new.range=c(0,1))

Figure3C1 <- gganatogram(data = anatomy.male[-18,],organism = 'human',sex='male',fillOutline = 'snow1',fill = 'value')+theme_void()+scale_fill_viridis(name='CCS change',option="cividis")
Figure3C2 <- gganatogram(data = anatomy.female,organism = 'human',sex='female',fillOutline = 'snow1',fill = 'value')+theme_void()+scale_fill_viridis(name='CCS change',option="cividis")

CCS_relative.df$CCS_change.log <- log2(CCS_relative.df$CCS_change)+1

# PAN-kidney 3: "KICH", "KIRP", "KIRC"
# PAN-Glioma 2: "GBM", "LGG"
# PAN-Squaumous 8: "SKCM", "LUAD", "LUSC", "THCA", "HNSC","PRAD","BLCA","ACC"
# PAN-GI 5: "ESCA","COAD", "STAD", "PAAD", "LIHC"
# PAN-GYN 5: "BRCA","OV", "UCS", "CESC", "UCEC"

Figure3B <- ggplot(CCS_relative.df, aes(y = (CCS_change.log), x = fct_reorder(twoClass.code,CCS_change,.fun =median,.desc = F), fill = twoClass.code)) +
    geom_boxplot(width=0.5,outlier.fill = NULL,outlier.size = NULL,show.legend = T,alpha=1,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL,outlier.shape = NULL) +
    geom_jitter(width=0.1, show.legend = F,size=0.5,aes(color=twoClass.code,fill=twoClass.code)) +
    stat_boxplot(geom='errorbar',width=0.2)+
    geom_boxplot(width=0.5,outlier.fill = NULL,outlier.size = 3,show.legend = F,alpha=1,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL,outlier.shape = NULL)  +
    geom_half_violin(side = "r",position = position_nudge(x=0.3))+
    scale_fill_manual(values = new.tumour.color$tissue.colors,name='Cancer types',breaks = new.tumour.color$twoClass.code) +
    scale_color_manual(values = new.tumour.color$tissue.colors,name='Cancer types',breaks = new.tumour.color$twoClass.code)+
    labs(y='Adjusted CCS change\n(log2)',title = 'Relative change of CCS',x='Cancer types') +
    coord_cartesian(ylim=c(0,12)) +
    theme_cowplot()+
    theme(legend.key = element_rect(fill='white'),
          strip.text.y =element_text(size = 10,face='bold'),
          axis.title.y = element_text(size = 16,face='bold'),
          axis.title.x = element_text(size = 16,face='bold'),
          plot.title = element_text(size = 17, face = "bold",hjust = c(0.5,0.5)),
          plot.subtitle = element_text(size = 12,hjust = c(0.5,0.5)),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),axis.text.x=element_text(angle = 45,hjust = c(1,1)),
          axis.text.y=element_text(size=20),
          legend.key.size = unit(8, "mm"),legend.title = element_text(size = 10,face='bold'))


Figure3C <- plot_grid(Figure3C1+theme(legend.position = 'none'),Figure3C2+theme(legend.position = 'none'),ncol = 2)
Figure3B.C <- plot_grid(Figure3B+theme(legend.position = 'none'),Figure3C,ncol=2,rel_widths = c(1.5,0.5),align = 'hv',label_size = 20,labels = c('B','C'))
Figure3.legends <- plot_grid(get_legend(Figure3A+guides(fill = guide_legend(nrow = 2))+
                                            theme(legend.position = 'bottom',legend.title =element_text(size = 20,face='bold'),
                                                  legend.key.size = unit(1.5,'cm'),legend.text =element_text(size = 12))),
                             NULL,
                             get_legend(Figure3C1+
                                            theme(legend.position = 'bottom',legend.title =element_text(size = 20,face='bold'),
                                                  legend.key.size = unit(1,'cm'),legend.text =element_text(size = 12))),nrow=1,rel_widths = c(1,0.3,0.5))


### Final Figure3

Figure3 <- plot_grid(Figure3A+theme(legend.position = 'none'),
                     Figure3.legends,
                     Figure3B.C,
                     label_size = 20,labels=c('A','',''),
                     ncol = 1,nrow=3,align = 'h',axis='lr',rel_heights = c(1,0.3,1))

