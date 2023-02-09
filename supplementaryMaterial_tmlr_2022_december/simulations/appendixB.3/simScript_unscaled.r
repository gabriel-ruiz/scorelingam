path <- NA
if(is.na(path)){
  print("Set path to directory of supplementary material!")
}
setwd(paste(path,'supplementaryMaterial_tmlr2022/',sep='/'))
#
#########################################################
  ######## DirectLiNGAM: Shimizu et. al (2011) ########
#########################################################
# for running python in R
require(reticulate) # devtools::install_cran('reticulate')
use_python("/usr/bin/python3")
# source a wrapper function for DirectLiNGAM
source_python(file='source/directlingam/wrapper.py') # will give error code if need to install python modules (note: many are required by DirectLingam's lingam module)
########################################################
  ######## HighDimLingam: Wang & Drton (2020) ########
########################################################
require(highDLingam) # devtools::install_github(repo='ysamwang/highDNG/highDlingam/')
require(doParallel)
ncores <- detectCores()-1 # for highDLingam parallelization
########################################################
  ######## ScoreLingam: TMLR Submission 2022 ########
########################################################
require(Rcpp)
Rcpp::sourceCpp('source/scorelingam/source.cpp')
source('source/scorelingam/helperFunctions.r')
#
########################################################
  ########## Architectures from bnlearn.com ##########
########################################################
path <-'simulations/section3.1/bnlearn architectures/' 
names <- list.files(path)
adjmats <- list() # will store adjmats to simulate data later
# names(adjmats) <- gsub(pattern='.rda',replacement='',x=names)
exclude <- c(1,
             which(names %in% paste(c('hailfinder','hepar2','win95pts','andes','sachs','asia'),'.rda',sep='') ))
for(nm in names[-exclude]){
  require(bnlearn)
  nm2use <- gsub(pattern='.rda',replacement='',x=nm)
  load(paste(path,nm,sep=''))
  adjmats[[nm2use]] <- as.matrix(igraph::as_adj(bnlearn::as.igraph(bn)) )
  rm(bn)
  print(nm2use)
  print(dim(adjmats[[nm2use]]))
}

########################
# example #
X <- genSCM.data(B=adjmats[[1]],n=1e5)
mb <- moralize(adjmats[[1]])
res <- sort_llrmbCPP(X,mb,numUpdates = 20)
check.valid.sort(res,t(adjmats[[1]]))
########################################################
  ############# Simulation Parameters ################
########################################################
pa.min <- 0.4
pa.max <- 0.9
b.min  <- 0.4
b.max  <- 0.7
#
famSpecified = c('laplace','laplace','laplace','logistic','t')
famTrue = c('laplace','logistic','t','logistic','t')
#
# FIXME: set to 1 but change to 30 if would like to replicate paper
nreps <- 30 
nfac2try <- c(0.5,1,2,10,50) # factor to multiply to ncol(adjmat) = p
#
names2use <- c('bn','n','p','rep','famSpec','famTrue','err.ours','err.dlingam','err.highD','tm.ours','tm.dlingam','tm.highD')
results <- data.frame( matrix(NA,nrow=0,ncol=length(names2use)) )
#
#
fname <- './simulations/appendixB.3/results.rdata'
###########################
#### simulation loops  ####
###########################
for( nm in rev(names(adjmats)) ){
  # if(nm == 'andes'){next()}
  B <- adjmats[[nm]]
  p <- ncol(B)
  nEdges <- sum(c(B!=0))
  B[B != 0] <- runif(n=nEdges,min=pa.min,max=pa.max) * rbinom(n=nEdges,size=1,prob=0.5)
  sc <- runif(n=p,min=b.min,max=b.max)
  # markov blanket
  mb <- moralize(B) # markov blankets
  skel <- skeleton(mb) # moral graph adjacency matrix
  #
  for(nfac in rev(nfac2try)){
    n <- ceiling(nfac*p)
    s <- paste('n=',n,sep='')
    # will store errors
    #### save work  #####
    save(results,file=fname)
    #
    for(f in 1:length(famSpecified)){
      for(r in 1:nreps){
        cat(paste('',nm,' (',famTrue[f],',',famSpecified[f],'): ','p=',p,'..n=',n,'..rep: ',r,'/',nreps,'..\n',sep=''))
        X <- genSCM.data(B=B,shape=sc,n=n,family=famTrue[f],df=10)
        # X <- scale(X) # standardize columns
        #####################
        ####high dim sort####
        cat('..highD..')
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        start <- Sys.time()
        perm.highd <- highDLingam::findGraphMulti(Y=X,verbose = F,B = skel)$topOrder
        end <- Sys.time()
        stopCluster(cl)
        tm.highD <- as.numeric( difftime(end,start,units='secs') )
        err.highD <- check.valid.sort(perm.highd,t(B))
        #####################
        ##DirectLiNGAM Sort##
        tm.dlingam <- NA
        err.dlingam <- NA
        if(n>p & nm != 'andes'){ #NOTE: DirectLiNGAM too slow on andes network
          cat('..DirectLingam..')
          start <- Sys.time()
          mod_lingam <- lingam_fit(X=X)
          perm.dlingam <- unlist(mod_lingam$causal_order_)+1
          end <- Sys.time(); 
          tm.dlingam <- as.numeric(difftime(end,start,units='secs'))
          err.dlingam <- check.valid.sort(perm.dlingam,t(B))
        }
        #####################
        ##ScoreLiNGAM Sort##
        cat('..ScoreLingam..\n')
        start <- Sys.time()
        perm.score <- sort_llrmbCPP(Xmat=X,mb=mb,numUpdates = 0,family=famSpecified[f],df=10)
        end <- Sys.time()
        tm.ours <- as.numeric(difftime(end,start,units='secs'))
        err.ours <- check.valid.sort(perm.score,t(B))
        #####################
        #### save work  #####
        res.curr <- data.frame(bn=nm,n=n,p=ncol(B),rep=r,famSpec=famSpecified[f],famTrue=famTrue[f],
                      err.ours=err.ours,err.dlingam=err.dlingam,err.highD=err.highD,tm.ours=tm.ours,tm.dlingam=tm.dlingam,tm.highD=tm.highD)
        results <- rbind(results,res.curr)
        save(results,file=fname)
        #####################
      }
    }
  }
}



