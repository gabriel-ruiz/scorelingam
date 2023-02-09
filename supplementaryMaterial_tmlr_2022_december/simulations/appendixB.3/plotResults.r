#
path <- NA
if(is.na(path)){
  print("Set path to directory of supplementary material!")
}
setwd(paste(path,'supplementaryMaterial_tmlr2022/',sep='/'))
#
load('simulations/section3.1/bnlearn architectures/adjmats.rdata')
fname <- 'simulations/appendixB.3/results.rdata'
#
load(fname)
dim(results)
#
results$ratio <- round(results$n/results$p,digits=1)
colnames(results)
##
library(reshape2)
res <- melt(data=results,id.vars=c('bn','n','p','rep','famSpec','famTrue','ratio'),
                measure.vars = c('err.ours','err.dlingam','err.highD'),value.name='error')
res2 <- melt(data=results,id.vars=c('bn','n','p','rep','famSpec','famTrue','ratio'),
                measure.vars = c('tm.ours','tm.dlingam','tm.highD'),value.name='time')
colnames(res)
res$algorithm <- ifelse(res$variable=='err.ours','ScoreLiNGAM',
                        ifelse(res$variable=='err.dlingam','DirectLiNGAM','HighDimLiNGAM'))
res$time <- res2$time
res$famSpec <- ifelse(res$famSpec=='laplace','Laplace',ifelse(res$famSpec=='logistic','Logistic',
               'Scaled-t (df=10)' ))
res$famTrue <- ifelse(res$famTrue=='laplace','Laplace',ifelse(res$famTrue=='logistic','Logistic',
               'Scaled-t (df=10)' ))
res$fam <- paste('Specified: ',res$famSpec,'\nTrue:',res$famTrue)
res <- dplyr::filter(res,!(bn %in% c('asia','sachs')))
res$bn <- paste(res$bn,'\n(p=',res$p,')',sep='')
res$network <- factor(res$bn,levels=unique(res$bn)[c(1,3,2)])
##
library(ggplot2)
# library(ggpubr)
res$famIndic <- ifelse(res$famSpec==res$famTrue,TRUE,FALSE)
#

p1a <- ggplot(dplyr::filter(res,famIndic==TRUE),aes(x=factor(ratio),y=error,fill=algorithm))+
  xlab('ratio: n/p')+ylab('proportion of parents sorted after a child')+
  geom_boxplot()+
  facet_grid(rows=vars(fam),cols=vars(network))+
  theme_classic()+
  ggtitle('Results without standardizing data matrix, known Markov blankets')+
  theme(
    plot.title = element_text( face='bold',size = 14,hjust=0.5),
    legend.title = element_text(face='bold', size = 12),legend.text=element_text( size = 12),
    axis.title = element_text(face = "bold", size = 14), axis.text = element_text( size = 13),
    strip.text = element_text(face='bold',size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = 'bottom',
    legend.position = 'bottom',
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank()
  )
p1a
ggsave(filename = '~/UCLA/nonGauss_scms/simulations/14 November 2022/sortError_correctSpec_unscaled.png',
       plot=p1a,units='in',width=6.75*2,height=4.35*2)
#
p1b <- ggplot(dplyr::filter(res,famIndic==FALSE),aes(x=factor(ratio),y=error,fill=algorithm))+
  xlab('ratio: n/p')+ylab('proportion of parents sorted after a child')+
  geom_boxplot()+
  facet_grid(rows=vars(fam),cols=vars(network))+
  theme_classic()+
  theme(
    # plot.title = element_text( size = 14),
    legend.title = element_text(face='bold', size = 12),legend.text=element_text( size = 12),
    axis.title = element_text(face = "bold", size = 14), axis.text = element_text( size = 13),
    strip.text = element_text(face='bold',size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = 'bottom',
    legend.position = 'bottom',
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank()
  )
p1b
####################################
########## Sorting Time ############
####################################
#
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
p2 <- ggplot(res,aes(x=factor(ratio),y=log10(time),fill=algorithm))+
  xlab('ratio: n/p')+ylab('log10(sorting time in seconds)')+
  geom_boxplot()+
  stat_summary(
    fun.y = median,
    geom = 'line',
    aes(group = algorithm, colour = algorithm),
    position = position_dodge(width = 0.9) #this has to be added
  )+ 
  facet_grid(cols=vars(fam),rows=vars(network))+
  theme_classic()+
  theme(
    # plot.title = element_text( size = 14),
    legend.title = element_text(face='bold', size = 12),legend.text=element_text( size = 12),
    axis.title = element_text(face = "bold", size = 14), axis.text = element_text( size = 13),
    strip.text = element_text(face='bold',size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = 'bottom',
    legend.position = 'bottom',
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank()
  )
p2
#

