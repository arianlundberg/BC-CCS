
# loading required data
######
load(file="~/Box Sync/Karolinska/Projects/CCS_2020/munge/Github/Gene.expression_data.RData")
load(file="~/Box Sync/Karolinska/Projects/CCS_2020/munge/Github/mutations_data.RData")
source(file="~/Box Sync/Karolinska/Projects/CCS_2020/munge/Github/functions.R")

#######
### Barplots of mutational data

### Group1
matrix.adj.G1.top20 <- melt(mutated.number.Group1.adjust[c('cancer',names(sort(colSums(mutated.number.Group1.adjust[,-1]),decreasing = T)))][1:21])
matrix.adj.G1.top20$value <- matrix.adj.G1.top20$value/100;matrix.adj.G1.top20 <- merge(matrix.adj.G1.top20,new.tumour.color[c('twoClass.code',"tissue.colors")],by.x='cancer',by.y='twoClass.code',sort=F)

### Group2
matrix.adj.G2.top20 <-  melt(mutated.number.Group2.adjust[c('cancer',names(sort(colSums(mutated.number.Group2.adjust[,-1]),decreasing = T)))][1:21])
matrix.adj.G2.top20$value <- matrix.adj.G2.top20$value/100;matrix.adj.G2.top20 <- merge(matrix.adj.G2.top20,new.tumour.color[c('twoClass.code',"tissue.colors")],by.x='cancer',by.y='twoClass.code',sort=F)


Figure4A <- ggbarplot(matrix.adj.G1.top20,x = 'variable',y = 'value',fill = 'cancer',color = 'cancer',ggtheme = theme_few())+
    scale_fill_manual(values = matrix.adj.G1.top20$tissue.colors,name='Cancer types',breaks = matrix.adj.G1.top20$cancer) +
    scale_color_manual(values = matrix.adj.G1.top20$tissue.colors,name='Cancer types',breaks = matrix.adj.G1.top20$cancer) +
    scale_y_continuous(name=expression(bold(paste("Adjusted mutational load (%)"))),labels = function(x) paste0(x, "%")) + scale_x_discrete(name="") +
    ggtitle(label = 'Group 1')+theme(axis.text.x=element_text(size = 8,angle = 45,hjust = c(1,1)),axis.title.y = element_text(size = 8,face='bold'),
                                     strip.text.x = element_text(size = 10, face = "bold"),strip.text.y = element_text(
                                         size = 7, face = "bold"),legend.position = 'bottom')+coord_cartesian(ylim = c(0,25))


Figure4B <- ggbarplot(matrix.adj.G2.top20,x = 'variable',y = 'value',fill = 'cancer',color = 'cancer',ggtheme = theme_few())+
    scale_fill_manual(values = matrix.adj.G2.top20$tissue.colors,name='Cancer types',breaks = matrix.adj.G2.top20$cancer) +
    scale_color_manual(values = matrix.adj.G2.top20$tissue.colors,name='Cancer types',breaks = matrix.adj.G2.top20$cancer) +
    scale_y_continuous(name=expression(bold(paste("Adjusted mutational load (%)"))),labels = function(x) paste0(x, "%")) + scale_x_discrete(name="") +
    ggtitle(label = 'Group 2')+theme(axis.text.x=element_text(size = 8,angle = 45,hjust = c(1,1)),axis.title.y = element_text(size = 8,face='bold'),
                                     strip.text.x = element_text(size = 10, face = "bold"),strip.text.y = element_text(
                                         size = 7, face = "bold"),legend.position = 'bottom')+coord_cartesian(ylim = c(0,25))

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

### not showing 0% from the plot
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

Figure4C <- ggplot(mut.data.melt,aes(x=cancer,y=(variable),fill=value))+geom_tile(color='black',width=0.9, height=0.8,size=0.1)+
    geom_text(aes(label=round.value),color='black',size=4) +
    labs(x='',y='')+scale_fill_gradient(low="white", high="skyblue4",guide = 'colorbar',na.value = NA,limits=c(0,75),breaks=c(0,25,50,75)
                                        ,name='Adj-mutation load (%)',labels = function(x) paste0(x, '%'))+
    facet_grid(~Group,scales = 'free_x',space='free',shrink = T)+
    theme_few()+theme(legend.text = element_text(size = 9,hjust = 1),
                      legend.position = 'right',strip.text = element_text(size=12,face = 'bold'),
                      axis.text.x  = element_text(size=10,face = 'bold'),axis.text.y  = element_text(size=10,face = 'bold'))


Figure4D <- ggplot(del.data.melt,aes(x=cancer,y=(variable),fill=value))+geom_tile(color='black',width=0.9, height=0.8,size=0.1)+
    geom_text(aes(label=round.value),color='black',size=4) +
    labs(x='',y='')+scale_fill_gradient(low="white", high="gold",guide = 'colorbar',na.value = NA,limits=c(0,75),breaks=c(0,25,50,75)
                                        ,name='Adj-Deletion (%)',labels = function(x) paste0(x, '%'))+
    facet_grid(~Group,scales = 'free_x',space='free',shrink = T)+
    theme_few()+theme(legend.text = element_text(size = 9,hjust = 1),
                      legend.position = 'right',strip.text = element_text(size=12,face = 'bold'),
                      axis.text.x  = element_text(size=10,face = 'bold'),axis.text.y  = element_text(size=10,face = 'bold'))


Figure4E <- ggplot(amp.data.melt,aes(x=cancer,y=(variable),fill=value))+geom_tile(color='black',width=0.9, height=0.8,size=0.1)+
    geom_text(aes(label=round.value),color='black',size=4) +
    labs(x='',y='')+scale_fill_gradient(low="white", high="firebrick3",guide = 'colorbar',na.value = NA,limits=c(0,75),breaks=c(0,25,50,75)
                                        ,name='Adj-Amplification (%)',labels = function(x) paste0(x, '%'))+
    facet_grid(~Group,scales = 'free_x',space='free',shrink = T)+
    theme_few()+theme(legend.text = element_text(size = 9,hjust = 1),
                      legend.position = 'right',strip.text = element_text(size=12,face = 'bold'),
                      axis.text.x  = element_text(size=10,face = 'bold'),axis.text.y  = element_text(size=10,face = 'bold'))

#########
# Generating the final figure4.

Fig4AB <- plot_grid(Figure4A+theme(legend.position=c(.85,.81),axis.title = element_text(face = 'bold'),plot.title = element_text(face = 'bold'),legend.key.size = unit(0.5, 'cm'),legend.text = element_text(size = 15,hjust = 1),legend.title =element_text(size = 17,hjust = 1,face = 'bold')),
                    Figure4B+theme(legend.position=c(.85,.82),axis.title = element_text(face = 'bold'),plot.title = element_text(face = 'bold'),axis.title.y=element_blank(),legend.key.size = unit(0.5, 'cm'),legend.text = element_text(size = 15,hjust = 1),legend.title =element_text(size = 17,hjust = 1,face = 'bold')),ncol = 2,align = 'vh',labels = c('A','B'))

Fig4CE <- plot_grid(Figure4C+theme(legend.position = 'bottom'),
                    Figure4D+theme(legend.position = 'bottom'),
                    Figure4E+theme(legend.position = 'bottom'),ncol = 3,align = 'vh',labels = c('C','D','E'))

Figure4 <- plot_grid(Fig4AB,
                     Fig4CE,
                     ncol = 1,nrow=2,align = 'hv',rel_heights = c(0.9,1.1))

