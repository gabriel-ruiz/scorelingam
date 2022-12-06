# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' The ScoreLiNGAM sorting procedure.
#' @param X An n by p data matrix which is used to estimate the permutation. 
#' @param mb A length p list object whose j-th entry gives the Markov blanket of node j (or its superset).
#' @param numUpdates The number of updates to give while sorting. Default: numUpdates=5.
#' @param family Currently allows 't', 'laplace' (default), or 'logistic'.
#' @param df Degrees of freedom for scaled-t-distributed noise. Default: df=10.
#' @return A length p vector with the specified topological ordering of the nodes in the underlying DAG. Its unique elements correspond to the indices 1,2,...,p. 
#' @export
scorelingam <- function(X, mb, numUpdates = 5L, family = "laplace", df = 10) {
  scorelingam_( X, mb, numUpdates, family, df)
}

#' Get weighted adjancency matrix for linear SEM.
#' Assumes that X is zero-centered, which can be done via X=scale(Xoriginal,center=T,scale=F).
#' Estimation is done via Ordinary Least Squares (OLS) regression.
#' @param X An n by p data matrix which is used to estimate the weighted adjacency matrix. 
#' @param pa A length p list object whose j-th entry gives the support for node j's parent set.
#' @return A p by p weighted adjacency matrix for the linear SEM. 
#' @export
getWeights <- function(X, pa) {
  getWeights_(X,pa) 
}

#' Quickly estimate a Pearson correlation matrix using RcppArmadillo.
#' @param X An n by p data matrix which is used to estimate the correlation matrix. 
#' @return A p by p estimate of the correlation matrix.
#' @export
corrMat <- function(X) {
  corrMat_(X) 
}







#' Generate a random weighted adjacency matrix with num.roots number of root nodes.
#' Each non-root node will have between pa.min and pa.max number of parents.
#' The parent weight (in absolute value) is between pa.wt.min and pa.wt.max.
#' The parent weight is positive with probability prob.pos (else, negative).
#' @param p Number of nodes in underlying DAG.
#' @param num.roots Number of root nodes in underlying DAG. Default: num.roots=p-1.
#' @param pa.min Minimum number of parents a node can have. Default: pa.min=1.  
#' @param pa.max Maximum number of parents a node can have. Default: pa.max=1.  
#' @param pa.wt.min The minimum (in absolute value) coefficient value for the parent to a child in the linear SEM. Default: pa.wt.min=1.
#' @param pa.wt.max The maximum (in absolute value) coefficient value for the parent to a child in the linear SEM. Default: pa.wt.max=1.
#' @param prob.pos The probability of a positive coefficient value for the parent to a child in the the linear SEM. Default: prob.pos=0.5.
#' @param perm A topological ordering for underlying DAG. Default is random.
#' @return A list with two elements: perm (the topological ordering used) and B (the weighted adjacency matrix).
#' @export
randomWeightedAdjacencyMatrix <- function( p,num.roots=p-1, pa.min=1, pa.max=1,pa.wt.min=1,pa.wt.max=1,
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



#' Generate Linear (laplace, logistic, or scaled-t error family) Structural Causal Model Data
#' @param B A p by p weighted adjacency matrix for an acyclic linear SEM. 
#' @param n The desired sample size. Default: n=1000.
#' @param family Currently allow 't', 'laplace' (default), or 'logistic'.
#' @param scale Length p vector s.t. \eqn{{\rm var}[\epsilon_j] = }scale[j]\eqn{\times {\rm var}[\epsilon_0]} with \eqn{\epsilon_0\sim {\rm family}(1) } (the chosen family with scale parameter 1). Default: scale=rep(1,ncol(B)).
#' @param df Degrees of freedom for scaled-t-distributed noise. Default: df=5.
#' @return An n by p matrix object.
#' @export
generateData <- function(B,n=1e3,family='laplace',scale=rep(1,ncol(B)),perm=topoSort(as.matrix(B)),df=5){
  #
  p <- ncol(B)
  assertthat::assert_that(ncol(B)==p && p == length(scale),msg=
                            'Error B not square matrix, or length(shape)!=ncol(B).')
  X <- NULL
  # add noise
  if(family=='laplace'){
    X <- matrix(rlap(n=n*p,mu=0,b=1),nrow=n,ncol=p)%*%diag(scale)
  }else if(family=='logistic'){
    X <- matrix(rlogis(n=n*p,location=0,scale=1),nrow=n,ncol=p)%*%diag(scale)
  }else{ # t-distributed random variables
    X <- matrix(rt(n=n*p,df=df,ncp=0),nrow=n,ncol=p)%*%diag(scale)
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




#' laplace r.v.
rlap <- function(n=1,mu=0,b=1){
  u <- runif(n = n,min=0,max=1)
  v <- runif(n = n,min=0,max=1)
  x <- mu + b*log(u/v)
  return(x)
}
#' 
#' laplace density of x
dlap <- function(x,center=mean(x),shape=mean(abs(x-center)),log=F){
  d <- -abs(x-center)/shape - log(2*shape) #log-density
  if(log==F){return(exp(d))}
  return(d)
}






#' Sorting Error. 
#' @param estOrder: A length p vector with the estimated topological ordering. Its unique elements correspond to the indices 1,2,...,p. 
#' @param A: The true DAG adjacency matrix. 
#' @return Returns proportion of parents sorted after children. See Equation (11) in the paper, which defines the output of this function: \eqn{{\rm ERR}(\hat{\pi})}.
#' @export
checkSortingErrors <- function(estOrder,A){
  M <- A
  d <- length(estOrder)
  p <- ncol(M)
  err.vec <- rep(0,d)
  err.possible <- rep(0,d)
  #
  unsorted <- setdiff(estOrder,1:p)
  #
  for(j in 1:(length(estOrder)-1) ){
    #what are the ancestors
    anc.j <- which(M[,estOrder[j]]!=0)
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



#' Moralize a DAG. 
#' Obtain Markov Blanket, the union of all children, parents, co-parents.
#' @param B The p by p weighted adjacency matrix. Assumes B[j,k]!=0 if \eqn{j\in PA_k}.
#' @return A list of length p. Each entry j corresponds to the Markov Blanket of node j according to adjacency matrix B.
#' @export
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
#' Skeleton for an undirected graph based on Markov Blanket (or neighborhood estimates) of each node.
#' @param mb The Markov blanket or neighborhood estiamtes. A list of length p.
#' @return A p by p adjacency matrix for an undirected graph. 0 indicates no undirected edge in the graph, while 1 indicates an directed edge in the graph.
#' @export
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


#' Get parent sets given Markov blanket and topological ordering
#' @param mb A length p list object whose j-th entry gives the Markov blanket of node j. The j-th entry can also be some other set that defines the possible support for the parent set.
#' @param ordering: A length p vector with the specified topological ordering of the nodes in the underlying DAG. Its unique elements correspond to the indices 1,2,...,p. 
#' @return A list object of length p; entry j denotes the parent set of node j in the underlying DAG.
#' @export
getParents <- function(mb,ordering){
  lapply(1:length(ordering),function(j){
    place <- which(ordering==j)
    if(place==1){return(numeric(0))}
    intersect(mb[[j]],ordering[1:(place-1)])
  })
}


#' #' Get weighted adjancency matrix and laplace scale parameters
#' #' @param X An n by p data matrix which is used to estimate the weighted adjacency matrix. 
#' #' @param mb A length p list object whose j-th entry gives the Markov blanket of node j. The j-th entry can also be some other set that defines the possible support for the parent set.
#' #' @param ordering A length p vector which gives the specified topological ordering of the nodes in the underlying DAG.
#' #' @return A list object of length 2. Entry Bhat gives the estimated p by p adjacency matrix, while B0hat is a length p vector containing the intercept terms, if any, in the linear SEM. 
#' getParams <- function(X,mbhat,ordering){
#'   p <- ncol(X)
#'   n <- nrow(X)
#'   #
#'   Bhat <- matrix(0,nrow=p,ncol=p)
#'   B0hat <- rep(0,p)
#'   lapScale <- rep(NA,p)
#'   sdResid <- rep(NA,p)
#'   rsq <- rep(NA,p); rsq[ordering[1]] <- 0
#'   rsq.adj <- rep(NA,p); rsq.adj[ordering[1]] <- 0
#'   #
#'   lapScale[ordering[1]] <- mean(abs(X[,ordering[1]]-median(X[,ordering[1]])))
#'   sdResid[ordering[1]] <- sd(X[,ordering[1]])
#'   B0hat[ordering[1]] <- median(X[,ordering[1]])
#'   #
#'   for(j in 2:p){
#'     mbjt <- intersect(ordering[1:(j-1)],mbhat[[ordering[j]]])
#'     if(length(mbjt)==0){
#'       lapScale[ordering[j]] <- mean(abs(X[,ordering[j]]-median(X[,ordering[j]])) )
#'       sdResid[ordering[j]] <- sd(X[,ordering[j]])
#'       #
#'       B0hat[ordering[j]] <- median(X[,ordering[j]])
#'       #
#'       rsq[ordering[j]] <- 0
#'       rsq.adj[ordering[j]] <- 0
#'       #
#'       next()
#'     }
#'     #
#'     mod <- lm(X[,ordering[j]] ~ X[,mbjt])
#'     smry <- summary(mod)
#'     #
#'     Bhat[mbjt,ordering[j]] <- coef(mod)[-1]
#'     B0hat[ordering[j]] <- coef(mod)[1]
#'     #
#'     resids <- X[,ordering[j]]- mod$fitted.values
#'     lapScale[ordering[j]] <- mean(abs(resids-median(resids)))
#'     sdResid[ordering[j]] <- sd(resids)
#'     rsq[ordering[j]] <- smry$r.squared
#'     rsq.adj[ordering[j]] <- smry$adj.r.squared
#'   }
#'   
#'   #
#'   return(list( 'Bhat'=Bhat,'B0hat'=B0hat))
#' }

###

