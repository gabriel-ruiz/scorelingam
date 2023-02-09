path <- NA
if(is.na(path)){
  print("Set path to directory of supplementary material!")
}
setwd(paste(path,'supplementaryMaterial_tmlr2022/',sep='/'))
#
require(data.table) # for fast reading in of data
meta <- fread(file='meta.tsv') # need to download from http://cells.ucsc.edu/?ds=allen-celltypes+mouse-cortex&meta=regionlabel
start <- Sys.time()
expr <- fread(file='exprMatrix.tsv') # need to download from http://cells.ucsc.edu/?ds=allen-celltypes+mouse-cortex&meta=regionlabel
end <- Sys.time(); difftime(end,start,units='mins')
#
genes <- expr[,'gene']$gene
# largest combo
toKeep <- meta$cellId[meta$injection_materials_label==''&meta$region_label=='VISp'&meta$class_label=='Glutamatergic']
length(toKeep)
# subset
expr <- expr[,..toKeep]
# make matrix
expr <- as.matrix(expr)
rownames(expr) <- genes
dim(expr)
# prop zeros
propZeros <- c()#rowMeans(expr==0)
for(j in 1:nrow(expr)){
  expr.j <- expr[j,]
  propZeros[j] <- mean(expr.j==0)
  if(j %% 1e3 == 0){
    cat(paste(j,'..',sep=''))
  }
}
cat('\n')
expr <- expr[propZeros<=0.5,]
dim(expr)
# transpose
expr <- t(expr)
# data splitting
set.seed(07192022)
inds.tr <- sample(nrow(expr),size=floor(nrow(expr)/2))
inds.ts <- setdiff(1:nrow(expr),inds.tr) 
inds.cor <- sample(inds.ts,size=floor(0.1*length(inds.ts)))
inds.ts <- setdiff(inds.ts,inds.cor)
#
X.tr <- expr[inds.tr,]
X.cor <- expr[inds.cor,]
X.ts <- expr[inds.ts,]
#
write.csv(x = X.tr,file = 'datTrain.csv')
write.csv(x = X.cor,file = 'datCor.csv')
write.csv(x = X.ts,file = 'datTest.csv')
#