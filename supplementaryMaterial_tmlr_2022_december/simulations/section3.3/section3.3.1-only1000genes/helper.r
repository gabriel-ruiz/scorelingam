crossEntropy <- function(X,B,s,laplace=T,parents){
  assertthat::assert_that(ncol(B)==ncol(X) && nrow(B)==ncol(B) && ncol(X)==length(s))
  p <- ncol(X)
  #
  # return joint cross entropy for laplace noise specification
  if(laplace){
    return(-sum(  
      sapply(1:p,function(j){
        paj <- parents[[j]]
        Rj <- X[,j]
        if(length(paj)!=0){
          Rj <- X[,j]-X[,paj]%*%matrix(B[paj,j],ncol=1)
        }
        mean( dlap(x = Rj,center=0,shape = s[j],log = T) )
      })
    ))
  }
  # return joint cross entropy for normal noise specification
  return(-sum(  
    sapply(1:p,function(j){
      paj <- parents[[j]]
      Rj <- X[,j]
      if(length(paj)!=0){
        Rj <- X[,j]-X[,paj]%*%matrix(B[paj,j],ncol=1)
      }
      mean( dnorm(x = Rj,mean=0,sd = s[j],log = T) )
    })
  ))
}

evalEntropy <- function(Xtest,Xtrain,numNbrs=50,prop.cor=0.2){
  #data splitting
  inds.cor <- sample(1:nrow(Xtrain),size=floor(prop.cor*nrow(Xtrain))) 
  inds.train <- setdiff(1:nrow(Xtrain),inds.cor)
  Xcor <- Xtrain[inds.cor,]
  Xtrain2 <- Xtrain[inds.train,]
  #
  cormat <- corCpp(Xmat=Xcor)
  mbhat <- lapply(1:ncol(Xcor),function(j){
    v <- (1:ncol(Xcor))[-j]
    v[ order(abs(cormat[j,-j]),decreasing=T)[1:numNbrs] ]
  })
  #results for our method using correlated neighbors
  start <- Sys.time()
  perm <- sort_llrmbCPP(Xmat=Xtrain2,mb=mbhat)
  end <- Sys.time();
  cat('\n')
  cat(paste('ScoreLiNGAM: ',difftime(end,start,units='secs'),'\n',sep=''))
  paHat <- getParents(mbhat,ordering=c(perm))
  B <- getWeights(X=Xtrain2,pa=paHat)
  R <- Xtrain2-Xtrain2%*%B
  b.train <- colMeans(abs(R))
  sd.train <- sqrt(colMeans(R^2))
  #results using a random permutation
  perm.rand <- sample(1:ncol(Xtrain2),size=ncol(Xtrain2),F)
  pa.rand <- getParents(mbhat,perm.rand)
  B.rand <- getWeights(Xtrain2,pa.rand)
  R.rand <- Xtrain2-Xtrain2%*%B.rand
  b.rand <- colMeans((abs(R.rand)))
  sd.rand <- sqrt(colMeans(R.rand^2))
  #results with sparsebn
  start <- Sys.time()
  dat <- sparsebnData(Xtrain,type='continuous')
  blacklist <- NULL
  #
  mod.sbn <- estimate.dag(data=dat,blacklist = blacklist)
  mod.selected <- select.parameter(x = mod.sbn,data = dat,type = 'profile',alpha = 0.1)
  end <- Sys.time();cat(paste('Sparsebn: ',difftime(end,start,units='secs'),'\n',sep=''))
  params.sbn <- estimate.parameters(mod.sbn,data=dat)[[mod.selected]]
  B.sbn <- as.matrix(params.sbn$coefs)
  sd.sbn <- sqrt(Matrix::diag(params.sbn$vars))
  perm.sbn <- topoSort(B.sbn)
  pa.sbn <- lapply(1:ncol(B.sbn),function(j){which(B.sbn[,j]!=0)})
  R.sbn <- Xtrain - Xtrain%*%B.sbn
  b.sbn <- colMeans(abs(R.sbn))
  ##results for our method using neighbors given by sparsebn
  paHat2 <- getParents(mbhat,ordering=c(perm.sbn))
  B2 <- getWeights(X=Xtrain2,pa=paHat2)
  R2 <- Xtrain2 - Xtrain2%*%B2
  b.train2 <- colMeans(abs(R2))
  sd.train2 <- sqrt(colMeans(R2^2))
  ###checking empirical cross-entropy
  return(
    data.frame(
      split = rep(c('train','test'),each=8),
      density = rep(c('Laplace','Gaussian','Laplace','Gaussian','Laplace','Gaussian','Laplace','Gaussian'),2),
      perm = rep(c('ScoreLiNGAM','ScoreLiNGAM','Random','Random','Sparsebn','Sparsebn','Sparsebn','Sparsebn'),2),
      neighborhood = rep(c('Correlation','Correlation','Correlation','Correlation','Own','Own','Correlation','Correlation'),2),
      method = rep(
        c('ScoreLiNGAM (Laplace)',
          'ScoreLiNGAM (Gaussian)',
          'Random (Laplace)',
          'Random (Gaussian)',
          'Sparsebn (Laplace)',
          'Sparsebn (Gaussian)',
          'Sparsebn (Laplace)',
          'Sparsebn (Gaussian)'),2),
      cross.entropy = c(
        # training data results
        crossEntropy(X=Xtrain,B=B,s=b.train,laplace=T,parents=paHat),
        crossEntropy(X=Xtrain,B=B,s=sd.train,laplace=F,parents=paHat),
        crossEntropy(X=Xtrain,B=B.rand,s=b.rand,laplace=T,parents=pa.rand),
        crossEntropy(X=Xtrain,B=B.rand,s=sd.rand,laplace=F,parents=pa.rand),
        crossEntropy(X=Xtrain,B=B.sbn,s=b.sbn,laplace=T,parents=pa.sbn),#
        crossEntropy(X=Xtrain,B=B.sbn,s=sd.sbn,laplace=F,parents=pa.sbn),
        crossEntropy(X=Xtrain,B=B2,s=b.train2,laplace=T,parents=paHat2),
        crossEntropy(X=Xtrain,B=B2,s=sd.train2,laplace=F,parents=paHat2),
        # testing data results
        crossEntropy(X=Xtest,B=B,s=b.train,laplace=T,parents=paHat),
        crossEntropy(X=Xtest,B=B,s=sd.train,laplace=F,parents=paHat),
        crossEntropy(X=Xtest,B=B.rand,s=b.rand,laplace=T,parents=pa.rand),
        crossEntropy(X=Xtest,B=B.rand,s=sd.rand,laplace=F,parents=pa.rand),
        crossEntropy(X=Xtest,B=B.sbn,s=b.sbn,laplace=T,parents=pa.sbn),
        crossEntropy(X=Xtest,B=B.sbn,s=sd.sbn,laplace=F,parents=pa.sbn),
        crossEntropy(X=Xtest,B=B2,s=b.train2,laplace=T,parents=paHat2),
        crossEntropy(X=Xtest,B=B2,s=sd.train2,laplace=F,parents=paHat2)
      )
    )
    
  )
}

###########

topoSort <- function(B){
  p <- ncol(B)
  assertthat::assert_that(p==nrow(B))
  perm <- rep(NA,p)
  B2use <- B
  unordered <- 1:p
  #
  ind <- 1
  while(ind <= p){
    roots.inds <- which(colSums(B2use!=0)==0)
    perm[ind:(ind+length(roots.inds)-1)] <- unordered[roots.inds]
    unordered <- setdiff(unordered,unordered[roots.inds])
    B2use <- as.matrix(B2use[-roots.inds,-roots.inds])
    #####
    ind <- ind+length(roots.inds)
  }
  return(perm)
}

###
getParents <- function(mb,ordering){
  lapply(1:length(ordering),function(j){
    place <- which(ordering==j)
    if(place==1){return(numeric(0))}
    intersect(mb[[j]],ordering[1:(place-1)])
  })
}
