-   <a href="#overview" id="toc-overview">Overview</a>
-   <a href="#install-package-from-github."
    id="toc-install-package-from-github.">Install package from Github.</a>
    -   <a href="#load-scorelingam-once-installed"
        id="toc-load-scorelingam-once-installed">Load scorelingam once
        installed</a>
-   <a href="#an-example-for-our-sorting-procedure-p5."
    id="toc-an-example-for-our-sorting-procedure-p5.">An example for our
    sorting procedure: <span class="math inline"><em>p</em> = 5</span>.</a>
    -   <a href="#generate-data." id="toc-generate-data.">Generate data.</a>
        -   <a href="#linear-sem-parameters." id="toc-linear-sem-parameters.">Linear
            SEM Parameters.</a>
        -   <a href="#data-matrix-with-n5000." id="toc-data-matrix-with-n5000.">Data
            matrix with <span class="math inline"><em>n</em> = 5000</span>.</a>
    -   <a
        href="#obtain-a-topological-ordering-estimate-and-check-proportion-of-parents-sorted-after-a-child."
        id="toc-obtain-a-topological-ordering-estimate-and-check-proportion-of-parents-sorted-after-a-child.">Obtain
        a topological ordering estimate and check proportion of parents sorted
        after a child.</a>
        -   <a href="#specify-the-neighborhood-sets."
            id="toc-specify-the-neighborhood-sets.">Specify the neighborhood
            sets.</a>
        -   <a href="#obtain-the-ordering-estimate"
            id="toc-obtain-the-ordering-estimate">Obtain the ordering estimate:</a>
        -   <a href="#check-the-accuracy-of-the-ordering"
            id="toc-check-the-accuracy-of-the-ordering">Check the accuracy of the
            ordering:</a>
        -   <a href="#check-accuracy-of-estimated-weighted-adjacency-matrix"
            id="toc-check-accuracy-of-estimated-weighted-adjacency-matrix">Check
            accuracy of estimated weighted adjacency matrix:</a>
-   <a href="#a-higher-dimensional-example-p10000."
    id="toc-a-higher-dimensional-example-p10000.">A higher dimensional
    example: <span class="math inline"><em>p</em> = 10, 000</span>.</a>
    -   <a href="#generate-data.-1" id="toc-generate-data.-1">Generate data.</a>
        -   <a href="#linear-sem-parameters.-1"
            id="toc-linear-sem-parameters.-1">Linear SEM Parameters.</a>
        -   <a href="#data-matrix-with-n5000.-1"
            id="toc-data-matrix-with-n5000.-1">Data matrix with <span
            class="math inline"><em>n</em> = 5000</span>.</a>
    -   <a
        href="#obtain-a-topological-ordering-estimate-and-time-the-algorithm."
        id="toc-obtain-a-topological-ordering-estimate-and-time-the-algorithm.">Obtain
        a topological ordering estimate and time the algorithm.</a>
        -   <a href="#specify-the-neighborhood-sets"
            id="toc-specify-the-neighborhood-sets">Specify the neighborhood
            sets:</a>
        -   <a href="#obtain-the-ordering-estimate-1"
            id="toc-obtain-the-ordering-estimate-1">Obtain the ordering
            estimate:</a>
    -   <a href="#check-accuracy-of-estimated-ordering"
        id="toc-check-accuracy-of-estimated-ordering">Check accuracy of
        estimated ordering:</a>

# Overview

This package implements the estimation of a topological ordering for a
Linear Structural Equation Model (SEM) with non-Gaussian errors, as
outlined in [Ruiz et. al
(2022+)](https://openreview.net/forum?id=4pCjIGIjrt):

> G Ruiz, OH Madrid-Padilla, Q Zhou. “Sequentially learning the
> topological ordering of directed acyclic graphs with likelihood ratio
> scores.” Transactions of Machine Learning Research, 2022+.
> \[[Link](https://openreview.net/forum?id=4pCjIGIjrt)\]

See below for some examples.

# Install package from Github.

    # install devtools, if needed. 
    require(devtools)
    # install the scorelingam package. 
    devtools::install_github(repo='gabriel-ruiz/scorelingam')

### Load scorelingam once installed

    library(scorelingam)

# An example for our sorting procedure: *p* = 5.

## Generate data.

### Linear SEM Parameters.

    p = 5;numRoots = 1; numParentsMin = 1; numParentsMax = 2
    scaleParam = runif(n=p,min=0.5,max=1.2)
    causalOrder = 5:1
    lingamParams = randomWeightedAdjacencyMatrix(p=p,num.roots=numRoots,pa.min=numParentsMin,pa.max=numParentsMax,
                        pa.wt.min = 0.25,pa.wt.max = 0.9,prob.pos = 0.5,perm = causalOrder)

    # weighted adjacency matrix
    print(lingamParams$B)

    ##           [,1]      [,2]      [,3]       [,4] [,5]
    ## [1,] 0.0000000 0.0000000 0.0000000  0.0000000    0
    ## [2,] 0.0000000 0.0000000 0.0000000  0.0000000    0
    ## [3,] 0.0000000 0.5048810 0.0000000  0.0000000    0
    ## [4,] 0.3853393 0.5883873 0.4772225  0.0000000    0
    ## [5,] 0.2603125 0.0000000 0.3124259 -0.8937634    0

### Data matrix with *n* = 5000.

When generating the data matrix, the possible options for the family
argument are ‘laplace’, ‘logistic’, and ‘t’, in which case the
additional argument df is needed (df=10 is the default)

    n = 5000
    X = generateData(B=lingamParams$B,scale = scaleParam,perm = causalOrder,n = n,family = 'laplace')
    dim(X)

    ## [1] 5000    5

## Obtain a topological ordering estimate and check proportion of parents sorted after a child.

When estimating a topological ordering for the linear SEM, the possible
options for the family argument are ‘laplace’, ‘logistic’, and ‘t’, in
which case the additional argument df is needed (df=10 is the default).

### Specify the neighborhood sets.

    # neighborhoods specified to be all other nodes
    mbhat = lapply(1:ncol(X),function(j){(1:ncol(X))[-j]}) 

### Obtain the ordering estimate:

    estOrder = scorelingam(X=X,mb=mbhat,numUpdates=ncol(X),family='laplace')
    print(estOrder)

    ##      [,1]
    ## [1,]    5
    ## [2,]    4
    ## [3,]    3
    ## [4,]    2
    ## [5,]    1

### Check the accuracy of the ordering:

    # check errors (ideally close to zero)
    checkSortingErrors(estOrder=estOrder,A=lingamParams$B)

    ## [1] 0

### Check accuracy of estimated weighted adjacency matrix:

    paHat = getParents(mb=mbhat,ordering=estOrder)
    (Bhat = getWeights(X=scale(X,scale=F,center=T),pa=paHat))

    ##              [,1]        [,2]      [,3]       [,4] [,5]
    ## [1,]  0.000000000  0.00000000 0.0000000  0.0000000    0
    ## [2,] -0.010345117  0.00000000 0.0000000  0.0000000    0
    ## [3,]  0.002221982  0.49206986 0.0000000  0.0000000    0
    ## [4,]  0.406513432  0.60945725 0.4780302  0.0000000    0
    ## [5,]  0.272617690 -0.01703144 0.3077930 -0.8893524    0

    # maximum entry-wise difference in absolute value
    norm(x=lingamParams$B-Bhat,type = 'i')

    ## [1] 0.0430518

# A higher dimensional example: *p* = 10, 000.

## Generate data.

### Linear SEM Parameters.

    p = 10000;numRoots = 100; numParentsMin = 1; numParentsMax = 2
    scaleParam = runif(n=p,min=0.5,max=1.2)
    causalOrder = c(seq(2,p,by=2),seq(1,p-1,by=2))
    lingamParams = randomWeightedAdjacencyMatrix(p=p,num.roots=numRoots,pa.min=numParentsMin,
                                                 pa.max=numParentsMax, pa.wt.min = 0.25,
                                                 pa.wt.max = 0.9,prob.pos = 0.5,perm = causalOrder)

### Data matrix with *n* = 5000.

When generating the data matrix, the possible options for the family
argument are ‘laplace’, ‘logistic’, and ‘t’, in which case the
additional argument df is needed (df=10 is the default)

    n = 5000
    X = generateData(B=lingamParams$B,scale = scaleParam,perm = causalOrder,n = n,family = 'laplace')
    dim(X)

    ## [1]  5000 10000

## Obtain a topological ordering estimate and time the algorithm.

When estimating a topological ordering for the Linear SEM, the possible
options for the family argument are ‘laplace’, ‘logistic’, and ‘t’, in
which case the additional argument df is needed (df=10 is the default).

### Specify the neighborhood sets:

    # neighborhoods specified to be true markov blankets
    mbhat = moralize(lingamParams$B) 
    #

### Obtain the ordering estimate:

    start = Sys.time()
    estOrder = scorelingam(X=X,mb=mbhat,numUpdates=10,family='laplace')
    end = Sys.time()
    difftime(end,start,units='mins')

    ## Time difference of 1.45492 mins

## Check accuracy of estimated ordering:

    # check sorting errors (ideally close to zero)
    # proportion of parents sorted after child
    checkSortingErrors(estOrder=estOrder,A=lingamParams$B)

    ## [1] 0.07580136
