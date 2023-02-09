path <- NA
if(is.na(path)){
  print("Set path to directory of supplementary material!")
}
setwd(paste(path,'supplementaryMaterial_tmlr2022/',sep='/'))
#
########################################################
  ######## ScoreLingam: TMLR Submission 2022 ########
########################################################
require(Rcpp)
Rcpp::sourceCpp('source/scorelingam/source.cpp')
source('source/scorelingam/helperFunctions.r')
#
##############################
fname <- 'simulations/section3.2/results.rdata'
p2try <- c(5000,10000)
nfac2try <- c(0.1,0.25,0.5)
nreps <- 50 
prop.roots <- 0.05
#
pa.min <- 1; pa.max <- 2
pa.wt.max <- 0.9; pa.wt.min <- 0.4
b.min <- 0.25; b.max <- 0.9
##
results <- data.frame(matrix(NA,nrow=0,ncol=6))
colnames(results) <- c('rep','p','n','mbhat','error','time')
##############################
iter <- 1
lastTime <- Sys.time()
for(p in rev(p2try)){
  for(nfac in rev(nfac2try)){
    n <- floor(nfac*p)
    for(rep in 1:nreps){
      #
      nroots <- floor(prop.roots*p)
      perm <- sample(1:p,size=p,F)
      b <- runif(n=p,min=b.min,b.max)
      B <- rand.wtd.adj.mat(p=p,num.roots=nroots,pa.min=pa.min,pa.max=pa.max,
              pa.wt.min=pa.wt.min,pa.wt.max=pa.wt.max,prob.pos=0.5,perm=perm)
      perm <- B$order
      B <- B$B
      mb <- moralize(B) # markov blankets
      # generate data
      X <- genSCM.data(B=B,perm=perm,shape=b,n=n,family='laplace')
      X <- scale(X) # standardize columns
      ############################################
      start <- Sys.time()
      est.perm <- sort_llrmbCPP(Xmat=X,mb=mb,numUpdates=0)
      end <- Sys.time()
      ##
      tmTrueMb <- as.numeric(difftime(end,start,units='secs'))
      errTrueMb <- check.valid.sort(est.perm,t(B))
      ############################################
      inds.cor <- sample(n,size=floor(0.2*n),F)
      inds.sort <- setdiff(1:n,inds.cor)
      corMat <- corCpp(Xmat=X[inds.cor,])
      numNbrs <- 10
      mbhat <- lapply(1:ncol(X),function(j){
        # most correlated variables 
        inds <- order(abs(corMat[j,-j]),decreasing=T)[1:numNbrs]
        return( ((1:p)[-j])[inds] )
      })
      start <- Sys.time()
      est.perm <- sort_llrmbCPP(Xmat=X[inds.sort,],mb=mbhat,numUpdates=0)
      end <- Sys.time()
      ##
      tmMbhat <- as.numeric(difftime(end,start,units='secs'))
      errMbhat <- check.valid.sort(est.perm,t(B))
      ############################################
      ### update results ###
      results <- rbind(results,
                  data.frame(rep=rep,p=p,n=n,mbhat=c('true','estimated'),
                    error=c(errTrueMb,errMbhat),time=c(tmTrueMb,tmMbhat)) )
      save(results,file=fname)
      ############################################
      cat(paste( 'Time: ',difftime(Sys.time(),lastTime,units='mins'),'...n=',n,'..p=',p,'..rep: ',rep,'/',nreps,'..\n',sep=''))
      #
      iter <- iter + 2
      lastTime <- Sys.time()
    }
    ### 
  }
}

###########################################
  ########### plot results ##############
load(fname)
results$neighborhood <- ifelse(results$mbhat=='estimated','10 Most Correlated','True Markov Blanket')
results$p2 <- factor(paste('p=',results$p,sep=''),levels=paste('p=',p2try,sep=''))
results$ratio <- results$n/results$p

library(ggplot2)
g1 <- ggplot(results,aes(x=factor(ratio),y=error,fill=neighborhood))+
  geom_boxplot()+facet_grid(rows=vars(p2),scales = 'free_x')+
  geom_boxplot()+
  xlab('ratio: n/p')+ylab('proportion of parents sorted after a child')+
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
g1
#
g2 <- ggplot(results,aes(x=factor(ratio),y=time,fill=neighborhood))+
  geom_boxplot()+facet_grid(cols=vars(p2),scales = 'free_x')+
  stat_summary(
    fun.y = median,
    geom = 'line',
    aes(group = neighborhood, colour = neighborhood),
    position = position_dodge(width = 0.9) #this has to be added
  )+
  xlab('ratio: n/p')+ylab('sorting time in seconds')+
  theme_classic()+
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
g2
#
