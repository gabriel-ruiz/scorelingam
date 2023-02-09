path <- NA
if(is.na(path)){
  print("Set path to directory of supplementary material!")
}
setwd(paste(path,'supplementaryMaterial_tmlr2022/',sep='/'))
#
########################################################
  ######## ScoreLingam: TMLR Submission 2022 ########
########################################################
library(Rcpp)
Rcpp::sourceCpp('source/scorelingam/source.cpp')
source('source/scorelingam/helperFunctions.r')
################
require(data.table)
#
X <- as.matrix(fread('simulations/section3.3/data/datTrain.csv')[,-1])
Xtest <- as.matrix(fread('simulations/section3.3/data/datTest.csv')[,-1])
X <- rbind(X,Xtest)
colnames(X) <- NULL
rm(Xtest)
n <- nrow(X)
p <- ncol(X)
# example
inds.tr <- sample(n,ceiling(0.5*n))
inds.ts <- setdiff(1:n,inds.tr)
#####
mbHatPos <- readRDS(file='simulations/section3.3/section3.3.2-allGenes/mbHatPos.rds')$mbHatPos
#####
#
namesCol <- c('r','method','rsq.ts','num.par','tmSort','Index')
res.ScoreLiNGAM <- data.frame(matrix(NA,nrow=0,ncol=length(namesCol)))
res.Random <- data.frame(matrix(NA,nrow=0,ncol=length(namesCol)))
colnames(res.ScoreLiNGAM) <- colnames(res.Random) <- namesCol
#####
nreps <- 3e1
fname.ScoreLiNGAM <- 'simulations/section3.3/section3.3.2-allGenes/resultsScoreLiNGAM_comparisonToRandom.rds'
fname.Random <- 'simulations/section3.3/section3.3.2-allGenes/resultsRandom_comparisonToRandom.rds'
#
for(r in 1:nreps){
  #
  inds.tr <- sample(n,ceiling(0.5*n))
  inds.ts <- setdiff(1:n,inds.tr)
  #
  start <- Sys.time()
  perm.r <- sort_llrmbCPP(Xmat=X[inds.tr,],mb = mbHatPos,numUpdates = 100,family = 'laplace')
  end <- Sys.time()
  tm.sort <- difftime(end,start,units='secs')
  #
  X.ts <- scale(X[inds.ts,])
  #
  perm.rand.r <- sample(1:p,size=p,replace=F)
  #
  pa.r <- getParents(mb=mbHatPos,ordering=perm.r)
  numPar.r <- unlist(lapply(pa.r,length))
  B.r <- getWeights(X=X.ts,pa = pa.r)
  #
  pa.rand.r <- getParents(mb=mbHatPos,ordering=perm.rand.r)
  numPar.rand.r <- unlist(lapply(pa.rand.r,length))
  B.rand.r <- getWeights(X=X.ts,pa = pa.rand.r)
  #
  rsq.ScoreLiNGAM <- getRsq(X.ts,B.r)
  rsq.Random <- getRsq(X.ts,B.rand.r)
  #
  rm(X.ts)
  #
  res.r.ScoreLiNGAM <- data.frame(r,'ScoreLiNGAM',rsq.ScoreLiNGAM,numPar.r,tm.sort=as.numeric(tm.sort), Matrix::invPerm(perm.r) )
  res.r.Random <- data.frame(r,'Random',rsq.Random,numPar.rand.r,0, Matrix::invPerm(perm.rand.r) )
  colnames(res.r.ScoreLiNGAM) <- colnames(res.r.Random) <- namesCol
  #
  res.ScoreLiNGAM <- rbind(res.ScoreLiNGAM,res.r.ScoreLiNGAM)
  res.Random <- rbind(res.Random,res.r.Random)
  if(r %% 5 == 0){
    saveRDS(res.ScoreLiNGAM,file=fname.ScoreLiNGAM,compress=F)
    saveRDS(res.Random,file=fname.Random,compress = F)
  }
  ####
  end <- Sys.time()
  cat(
    paste(r,': ',difftime(end,start,units='mins'),'\n.\n.\n.',sep='')
  )
}

###

library(ggplot2)
res.Random <- readRDS(file=fname.Random)
res.ScoreLiNGAM <- readRDS(file=fname.ScoreLiNGAM)
#
results <- rbind(res.Random,res.ScoreLiNGAM)
results.adj <- results; results.adj$rsq.ts <- 1-(1-results$rsq.ts)*(length(inds.ts)-1)/(length(inds.ts)-results$num.par)
results <- rbind(results,results.adj)
results$type <- rep(c('Unadjusted','Adjusted'),each=nrow(results.adj))
#
colnames(results)
res.med <- aggregate(rsq.ts~method+r+type,FUN=median,data=results)
res.q80 <- aggregate(rsq.ts~method+r+type,FUN=function(v){quantile(v,probs=0.8)},data=results)
res.q90 <- aggregate(rsq.ts~method+r+type,FUN=function(v){quantile(v,probs=0.9)},data=results)
res.q95 <- aggregate(rsq.ts~method+r+type,FUN=function(v){quantile(v,probs=0.95)},data=results)
res.qs <- rbind(res.med,res.q80,res.q90,res.q99)
res.qs$percentile <- rep(c('50th Percentile','80th Percentile','90th Percentile','95th Percentile'),each=nrow(res.med))
#
(g1 <- ggplot(dplyr::filter(res.qs,percentile!='99th Percentile'))+
  geom_boxplot(aes(y=rsq.ts,x=method))+
  facet_grid(cols=vars(type),rows=vars(percentile),scales = 'free_y')+
  ylab('Coefficient of Determination\'s Quantiles\nAcross all 10,012 Genes')+
  xlab('\nPermutation Method')+
  theme_classic()+
  scale_y_continuous(n.breaks = 10)+
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
  )+
  theme(axis.text.x = element_text(angle = -10))
)



