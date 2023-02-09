# 50 top correlated variables as nbrs
# 1000 random genes 
path <- NA
if(is.na(path)){
  print("Set path to directory of supplementary material!")
}
setwd(paste(path,'supplementaryMaterial_tmlr2022/',sep='/'))
# 
########################################################
  ######### sparsebn: Aragam et. al (2019)  ##########
########################################################
require(sparsebn) # devtools::install_github(repo='itsrainingdata/sparsebn')

########################################################
  ######## ScoreLingam: TMLR Submission 2022 ########
########################################################
require(Rcpp)
Rcpp::sourceCpp('source/scorelingam/source.cpp')
source('source/scorelingam/helperFunctions.r')
#
source('simulations/section3.3/section3.3.1-only1000genes/helper.r') # helper function for making comparisons
################################
require(data.table)
X <- fread('simulations/section3.3/data/datTrain.csv')
Xtst <- fread('simulations/section3.3/data/datTest.csv')
X <- rbind(X,Xtst); rm(Xtst)
dim(X)
X <- as.data.frame(X)
X <- X[,-1]
################################
set.seed(seed=12292021)
numGenes <- 1000
propZero <- colMeans(X==0)
colInds <- sample(ncol(X),numGenes,F)
#
nreps <- 50
#
cnames <- c('split','density','perm','neighborhood','method','cross.entropy','rep')
results <- data.frame(matrix(nrow=0,ncol=length(cnames)))
colnames(results) <- cnames
#
for(r in 1:nreps){
  start <- Sys.time()
  sz <- 2000
  # sz <- nrow(X)
  inds <- sample(nrow(X),size=sz)
  inds.train <- inds[1:floor(sz/2)]
  inds.test <- inds[(floor(sz/2)+1):(sz)]
  ###
  Xtrain <- scale(X[inds.train,colInds],scale = T)
  Xtesting <- scale(X[inds.test,colInds],scale = T)
  ###
  tmp <- capture.output( res.r <- evalEntropy(Xtest=Xtesting,Xtrain=Xtrain,numNbrs=50,prop.cor=0.2) )
  res.r$rep <- r
  ###
  results <- rbind(results,res.r) 
  print(dplyr::arrange(dplyr::filter(res.r,split=='test'),cross.entropy))
  print(difftime(Sys.time(),start,units='mins'))
  ###
  save(results,file='simulations/section3.3/section3.3.1-only1000genes/results.rdata')
}

##############################################3

#
####For TMLR paper####
load('simulations/section3.3/section3.3.1-only1000genes/results.rdata')
results$split <- ifelse(results$split=='test','Test Set','Training Set')
results$neighborhood <- as.character(results$neighborhood)
results$neighborhood[results$neighborhood=='Own'] <- 'Sparsebn'
results$neighborhood <- as.factor(results$neighborhood)

library(ggplot2)
#
#
results$method <- paste('(',results$perm,',',results$neighborhood,')',sep='')
p <- ggplot(dplyr::filter(results,perm!='Random'&split!='Training Set'&method!='(Sparsebn,Correlation)'),      
       aes(x=density,y=-cross.entropy,fill=perm) )+
  geom_boxplot()+ylab('Mean Log-Likelihood (Test Set)')+xlab('Residual Density Specification')+
  theme_bw()+labs(fill='Permutation Method')+
  theme(axis.text.x = element_text(angle= -00),legend.position = 'bottom')+
  theme(
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

p
