#####
# get gene neighborhoods for sorting procedure based on non-negative least squares
#####
path <- NA
if(is.na(path)){
  print("Set path to directory of supplementary material!")
}
setwd(paste(path,'supplementaryMaterial_tmlr2022/',sep='/'))
####
require(data.table)
Xcor <- fread('simulations/section3.3/data/datCor.csv')
Xcor <- scale(as.matrix(Xcor[,-1]))

# neighborhood non-negative least squares regression
mbPos <- function(X,numUpdates=100,parallel=F){
  require(glmnet)
  mbHat <- list()
  p <- ncol(X); n <- nrow(X)
  if(!parallel){
    for(j in 1:p){
      mod <- glmnet(y=X[,j],x=X[,-j],lambda=c(0),lower.limits = rep(0,p-1) )
      mbHat[[j]] <- ((1:p)[-j])[which(mod$beta[,1]!=0)]
      if(j==1||j==p||j%%numUpdates==0){
        cat(paste(j,'..',sep=''))
      }
    };cat('\n')
  }else{
    require(future.apply)
    mbHat <- future.apply::future_sapply(1:p,function(j){
      mod <- glmnet(y=X[,j],x=X[,-j],lambda=c(0),lower.limits = rep(0,p-1) )
      return( 
        ((1:p)[-j])[which(mod$beta[,1]!=0)]
      ) 
    })
  }
  return(mbHat)
}
#
start <- Sys.time()
mbHatPos <- mbPos(Xcor,numUpdates = 1e3,parallel = F)
end <- Sys.time()
difftime(end,start)
saveRDS(mbHatPos,file='simulations/section3.3/section3.3.2-allGenes/mbHatPos.rds',compress=F)
#