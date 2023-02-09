#generate a random wtd. adj. matrix with num.roots number of root nodes
#and each non-root with b/w pa.min and pa.max number of parents
#parent weight (abs value) is b/w pa.wt.min and pa.wt.max
#the parent weight is positive with probability prob.pos
rand.wtd.adj.mat <- function( p,num.roots=p-1, pa.min=1, pa.max=1,pa.wt.min=1,pa.wt.max=1,
              prob.pos=0.5,perm=sample(x = p,size = p,replace = F)){
  ######checks######
  assertthat::assert_that(prob.pos >= 0);assertthat::assert_that(prob.pos <= 1)
  assertthat::assert_that(length(perm)==p)
  assertthat::assert_that(num.roots < p,msg = "Error: too many root nodes (>= p).")
  assertthat::assert_that(num.roots >= 1,msg = "Error: too little root nodes (< 1).")
  ######
  #Upper triangular wtd. adj. matrix (unpermuted)
  B <- matrix(0,ncol=p,nrow=p)
  #
  for(j in (num.roots+1):p){
    num.pa <- sample(x = pa.min:pa.max,size=1)
    if(num.pa > j-1){num.pa <- j-1}
    #indices of the parents
    inds.pa <- sample(1:(j-1),size=num.pa,replace=F)
    #parent weights
    wts.pa <- runif(n = num.pa,min = pa.wt.min,max = pa.wt.max)
    #sign of ancestor weights
    sign.pa <- ifelse(rbinom(n=num.pa,size=1,prob=prob.pos)==1,1,-1)
    #assigning to mixing matrix
    B[inds.pa,j] <- sign.pa * wts.pa
    #
  }
  #
  #
  invPerm <- 1:p
  invPerm[perm] <- 1:p
  #
  return(list('order'=perm,'B'=B[invPerm,invPerm]))
  # return(list('order'=perm,'M'=M[invPerm,invPerm],'B'=B[invPerm,invPerm]))
}



# generate Linear (laplace error) Structural Causal Model Data
# M := p by p mixing matrix, upto permutation lower triangular
# n := desired sample size
# noise.density := length p vector s.t. noise_unscaled[j] ~ noise.density[j]
## currently allow 't' and 'laplace'
# scale := length p vector s.t. var[noise[j]] = scale[j]*var[noise_unscaled[j]]
# df := degrees of freedom for e.g. t-dist noise
genSCM.data <- function(B,shape=rep(1,ncol(B)),perm=topoSort(as.matrix(B)),n=1e3,family='laplace',df=5){
  #
  p <- ncol(B)
  assertthat::assert_that(ncol(B)==p && p == length(shape),msg=
                'Error B not square matrix, or length(shape)!=ncol(B).')
  X <- NULL
  # add noise
  if(family=='laplace'){
    X <- matrix(rlap(n=n*p,mu=0,b=1),nrow=n,ncol=p)%*%diag(shape)
  }else if(family=='logistic'){
    X <- matrix(rlogis(n=n*p,location=0,scale=1),nrow=n,ncol=p)%*%diag(shape)
  }else{ # t-distributed random variables
    X <- matrix(rt(n=n*p,df=df,ncp=0),nrow=n,ncol=p)%*%diag(shape)
  }
  ####
  notDone <- T
  BCurr <- B
  nPars <- colSums(BCurr != 0)
  currRoots <- which(nPars==0)
  currNonRoots <- which(nPars!=0 )
  lastRoots <- c()
  while(notDone){
    #
    X[,currNonRoots] <- X[,currNonRoots]+X[,currRoots]%*%matrix(BCurr[currRoots,currNonRoots],ncol=length(currNonRoots))
    ####
    BCurr[currRoots,currNonRoots] <- 0
    ####
    nPars[currNonRoots] <- colSums( matrix(BCurr[,currNonRoots] != 0,ncol=length(currNonRoots)))
    currRoots <- currNonRoots[which(nPars[currNonRoots]==0)]
    currNonRoots <- currNonRoots[which(nPars[currNonRoots] != 0)]
    #
    if(length(currNonRoots)==0){notDone <- F}
  }
  return(X)
}

 


#laplace r.v.
rlap <- function(n=1,mu=0,b=1){
  u <- runif(n = n,min=0,max=1)
  v <- runif(n = n,min=0,max=1)
  x <- mu + b*log(u/v)
  return(x)
}
#
#laplace density of x
dlap <- function(x,center=mean(x),shape=mean(abs(x-center)),log=F){
  d <- -abs(x-center)/shape - log(2*shape) #log-density
  if(log==F){return(exp(d))}
  return(d)
}



#############



## Returns proportion of ancestors sorted after descendants
## if argument M is a DAG adjancency matrix, then returns 
## proportion of parents sorted after children
#M: the mixing matrix
check.valid.sort <- function(estOrder,M){
  d <- length(estOrder)
  p <- ncol(M)
  err.vec <- rep(0,d)
  err.possible <- rep(0,d)
  #
  unsorted <- setdiff(estOrder,1:p)
  #
  for(j in 1:(length(estOrder)-1) ){
    #what are the ancestors
    anc.j <- which(M[estOrder[j],]!=0)
    if(length(anc.j)==0){next()}#no ancestors
    #count number of ancestors that come after it in the estimated sort
    err.vec[j] <- sum( anc.j %in% c(estOrder[(j+1):d] ,unsorted) )
    # err.vec[j] <- sum( anc.j %in% estOrder[(j+1):d]  )+ sum( !(anc.j %in% estOrder) )
    #number of possible errors = number of ancestors
    err.possible[j] <- length(anc.j)
  }
  #metric between 0 and 1
  # return(sum(err.vec)/sum(err.possible+sum(M[est.order[d],]!=0)))
  return(sum(err.vec)/sum(M!=0))
  # err.possible[err.possible==0] <- 1
  # return(err.vec/err.possible)
}



# obtain markov blankets for each node based on moralization of DAG
# returns a list, where mb[[k]] is markov blanket of node k
moralize <- function(B){ #markov blanket for each node
  p <- ncol(B)
  #
  mb <- list()
  #
  for(j in 1:p){
    pa <- which(B[,j] !=0 ) #parents
    ch <- which(B[j,] !=0 ) #children
    co_pa <- c()
    if(length(ch)>1){ #co-parents
      co_pa <- unlist( sapply(X=ch,function(k){which(B[,k]!=0)}) ) 
    }
    #
    mb[[j]] <- setdiff(sort( unique(c(pa,ch,co_pa) ) ),j)
  }
  #
  return(mb)
}

#
#prior skeleton based on mb#
#to input into directLingam#
#0 indicates to edge in dag
#1 indicates directed path
#-1 indicates unkown
skeleton <- function(mb){
  p <- length(mb)
  skel <- matrix(data = 0,nrow = p,ncol = p)
  #
  for(j in 1:p){
    skel[j,mb[[j]]] <- 1
    skel[mb[[j]],j] <- 1
  }
  diag(skel) <- 0
  return(skel)
}


# get parent sets given markov blanket and topological ordering
getParents <- function(mb,ordering){
  lapply(1:length(ordering),function(j){
    place <- which(ordering==j)
    if(place==1){return(numeric(0))}
    intersect(mb[[j]],ordering[1:(place-1)])
  })
}


# get weighted adjancency matrix and laplace scale parameters
#####
getParams <- function(X,mbhat,ordering,laplace=F){
  p <- ncol(X)
  n <- nrow(X)
  #
  Bhat <- matrix(0,nrow=p,ncol=p)
  B0hat <- rep(0,p)
  lapScale <- rep(NA,p)
  sdResid <- rep(NA,p)
  rsq <- rep(NA,p); rsq[ordering[1]] <- 0
  rsq.adj <- rep(NA,p); rsq.adj[ordering[1]] <- 0
  #
  lapScale[ordering[1]] <- mean(abs(X[,ordering[1]]-median(X[,ordering[1]])))
  sdResid[ordering[1]] <- sd(X[,ordering[1]])
  B0hat[ordering[1]] <- median(X[,ordering[1]])
  #
  if(laplace==F){ 
    for(j in 2:p){
      mbjt <- intersect(ordering[1:(j-1)],mbhat[[ordering[j]]])
      if(length(mbjt)==0){
        lapScale[ordering[j]] <- mean(abs(X[,ordering[j]]-median(X[,ordering[j]])) )
        sdResid[ordering[j]] <- sd(X[,ordering[j]])
        #
        B0hat[ordering[j]] <- median(X[,ordering[j]])
        #
        rsq[ordering[j]] <- 0
        rsq.adj[ordering[j]] <- 0
        #
        next()
      }
      #
      mod <- lm(X[,ordering[j]] ~ X[,mbjt])
      smry <- summary(mod)
      #
      Bhat[mbjt,ordering[j]] <- coef(mod)[-1]
      B0hat[ordering[j]] <- coef(mod)[1]
      #
      resids <- X[,ordering[j]]- mod$fitted.values
      lapScale[ordering[j]] <- mean(abs(resids-median(resids)))
      sdResid[ordering[j]] <- sd(resids)
      rsq[ordering[j]] <- smry$r.squared
      rsq.adj[ordering[j]] <- smry$adj.r.squared
    }
  }else{
      for(j in 2:p){
        require(quantreg)
        mbjt <- intersect(ordering[1:(j-1)],mbhat[[ordering[j]]])
        #
        if(length(mbjt)==0){
          lapScale[ordering[j]] <- mean(abs(X[,ordering[j]]-median(X[,ordering[j]])) )
          sdResid[ordering[j]] <- sd(X[,ordering[j]])
          #
          rsq[ordering[j]] <- 0
          rsq.adj[ordering[j]] <- 0
          #
          next()
        }
        #
        rqfit <- rq(X[,ordering[j]] ~ X[,mbjt],tau=0.5)
        #
        Bhat[mbjt,ordering[j]] <- coef(rqfit)[-1]
        B0hat[ordering[j]] <- coef(rqfit)[1]
        #
        resids <- rqfit$residuals
        lapScale[ordering[j]] <- mean(abs(resids-mean(resids)))
        sdResid[ordering[j]] <- sd(resids)
        rsq[ordering[j]] <- 1-var(resids)/var(X[,ordering[j]])
        rsq.adj[ordering[j]] <- 1-(1-rsq[ordering[j]])*(n-1)/(n-1-length(mbjt))
      }
  }
  #
  return(list( 'Bhat'=Bhat,'B0hat'=B0hat,'lapScale'=lapScale,'sdResid'=sdResid,
               'rsq'=rsq,'rsq.adjusted'=rsq.adj))
}

###

