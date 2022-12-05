#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma; 
using namespace std;

arma::mat fastLm(const arma::vec & y, const arma::mat & X) {
  //
  int n = X.n_rows, k = X.n_cols;
  
  arma::mat coef = arma::solve(X, y); 
  arma::mat resid = y - X*coef; 
  return(resid);
  
}


// laplace density 
arma::vec dlap(const arma::vec & y,const double & mu, const double & b, const bool & log){
    arma::vec d = -arma::abs( y-mu )/b - std::log(2*b);
    if(log==false){
        return(arma::exp(d));
    }
    return(d);
}


// normal density
arma::vec dnorm(const arma::vec & y,const double & mu, const double & sigma, const bool & log){
    double sqrtpitimes2 = std::pow(2*M_PI,0.5); 
    double c = std::log(sigma* sqrtpitimes2  );
    arma::vec d = -arma::square(( y-mu )/(sigma))/2 - c;
    if(log==false){
        return(arma::exp(d));
    }
    return(d);
}
// logistic density
arma::vec dlogis(const arma::vec & y,const double & mu, const double & sigma, const bool & log){
    double s = std::pow( 3 * sigma * sigma / (M_PI * M_PI), 0.5 );
    arma::vec c = arma::log( s * arma::pow( 1 + arma::exp( -1*(y-mu)/s ),2 ) );
    arma::vec d = -1*(y-mu)/s - c;
    if(log==false){
        return(arma::exp(d));
    }
    return(d);
}

// scaled t density
arma::vec dt(const arma::vec & y,const double & mu, const double & sigma, const bool & log,const double & df=10){
    double s = std::pow( (df-2) * sigma * sigma / df ,0.5 );
    double c = std::log( tgamma( (df+1)/2 )/ (std::pow(df*M_PI * s * s,0.5) * tgamma( df/2 ) ) );
    arma::vec d = -1*(df+1)/2 * arma::log( 1+ arma::pow((y-mu)/s,2)/df ) + c;
    if(log==false){
        return(arma::exp(d));
    }
    return(d);
}

// likelihood ratio used in the algorithm
double llr(const arma::vec & y,const std::string & family="laplace",const double & df=10){ //family = "laplace", "logistic", "t"
    double mu = 0.0;
    double b = arma::mean( arma::abs( y ) );
    double sigma = std::pow( arma::mean( arma::square( y ) ), 0.5 );
    double val = 0;
    //
    if(family=="logistic"){//Logistic case
      val = arma::mean( dlogis(y,mu,sigma,true) - dnorm(y,mu,sigma,true) );
    }
    else if(family=="t"){// Scaled-t case
      val = arma::mean( dt(y,mu,sigma,true,df) - dnorm(y,mu,sigma,true) );
    }
    else{// Otherwise assume Laplace case
        val = arma::mean( dlap(y,mu,b,true)-dnorm(y,mu,sigma,true) )  ;
    }
    return(val);
}



// the sorting procedure used in practice
// [[Rcpp::export()]]
arma::uvec sortllr(const arma::mat & Xmat, const std::vector <std::vector<int>> & mb, const int & numUpdates=5,
                         const std::string & family="laplace", const double & df = 10){//family = "laplace", "logistic", "t" (df is for this case)
    int p = Xmat.n_cols;
    int n = Xmat.n_rows;
    
    // mixing matrix s.t. Xmat = Emat*(M^T)
    arma::mat M = zeros<arma::mat>(p,p);
    M.diag().ones(); // populate diagonal with ones

    // where the ordering will be stored
    arma::uvec perm(p);

    //
    std::vector<unsigned int> unordered_nodes(p);//keep track of unordered nodes
    std::vector<std::vector<unsigned int>> ancestors(p);//ancestor set for each node
    // arma::vec unordered_nodes(p);
    // insert some values: include itself for the sake of this code
    for(int i=0; i<p; ++i){
        unordered_nodes[i] = i;  // 0 1 2 3 4 ... p-1
        ancestors[i].push_back(i);  // 0 1 2 3 4 ... p-1
    }
    //standardizing X
    arma::mat X = Xmat;
    arma::mat mu(1,p);
    arma::mat sd(1,p);
    for(int i=0;i<p;++i){
      arma::vec Xi = X.col(i);
      mu(0,i) = arma::mean(Xi);
      sd(0,i) = arma::stddev(Xi);
    }
    arma::mat unos = ones<mat>(n,1);
    X -= unos*mu;
    X *= arma::diagmat(1/sd);
    //keep track of ordered neighbors
    std::vector <std::vector<int>> mbt;
    for(int k=0; k < p; ++k ){
        mbt.push_back(std::vector<int>()); //start off as empty
    }
    //reisduals used to compute likelihood ratio
    arma::mat resids = X;
    arma::vec llr_k(p);
    for(unsigned int k=0; k < p; ++k){
        arma::vec residsk = resids.col(k);
        llr_k(k) = llr(residsk,family,df);
    }
    //
    double negInf = std::numeric_limits<double>::min(); // value to set llr_k[k] when node is ordered
    //
    for(int t=0; t < p; ++t){
        //have to convert to use functionality of arma::uvec
        arma::uvec unordered_now = arma::conv_to<arma::uvec>::from(unordered_nodes);
        int k_pr_ind = llr_k(unordered_now).index_max();
        //the next node to continue the ordering
        perm(t) = unordered_nodes[k_pr_ind];
        unordered_nodes.erase(unordered_nodes.begin()+k_pr_ind);
        llr_k(perm(t)) = negInf;//so that this node not selected again
        //
        if(unordered_nodes.size()>1){
          std::vector<int> neigh_t = mb[perm(t)] ;// markov blanket of node perm(t)
          arma::uvec curr_ord = perm( span(0,t) );// partial ordering after perm(t) is appended
          // the ancestors of node perm(t), including itself for the sake of the below update to mixing matrix
          std::vector<unsigned int> anc_t = ancestors[perm(t)];
          //update residuals and llr for perm(t)'s neighbors
          for(unsigned int j=0;j < neigh_t.size();j++){
            unsigned int k = neigh_t[j]- 1;// subtract 1 b/c mb[ perm(t) ] is 1-indexed from R
            if(mb[k].size()==0 || llr_k(k)==negInf){
              continue;
            }
            // update mixing matrix entries in row k and residual.col(k)
            for(unsigned int l=0;l<anc_t.size();l++){
              unsigned int anc = anc_t[l];
              if( M(perm(t),anc) !=0 & M(k,anc) == 0 ){
                arma::mat coef = solve( resids.col(anc), resids.col(k) );// regress X[,k] on resids[,perm(t)]
                M(k,anc) = coef(0,0); // signal X_k gets from epsilon_k in LiNGAM
                resids.col(k) -= resids.col(anc)*M(k,anc); // update residual
                //
                ancestors[k].push_back(anc);// ancestor of perm(t) must also be an ancestor of node k
              }
            }
            // update likelihood ratio now
            llr_k(k) = llr(resids.col(k),family,df);
          }
          //
        }
        //
        if(numUpdates > 0){
          if(t % (p/numUpdates) == 0){
            //sleep(10);
            
            cout <<  t << ".." << flush;
            //
          }
        }
    }
    cout << p << endl; 
    return  perm + 1;// add 1 because R 1-indexes vectors
}

// weighted adjancency matrix based on known/estimated parent sets
// [[Rcpp::export()]]
arma::mat getWeights(const arma::mat & X,const std::vector <std::vector<unsigned int>> & pa){
  unsigned int p = X.n_cols;
  arma::mat B(p,p);
  //
  for(unsigned int i=0;i<p;++i){
    arma::uvec node(1);
    node(0) = i;
    if(pa[i].size()==0){
      continue;
    }
    //
    arma::uvec pa2use = conv_to<arma::uvec>::from(pa[i])-1;//-1 because parent indices are 1-indexed 
    B(pa2use,node) = arma::solve(X.cols(pa2use), X.col(i)); 
  }
  return(B);
}

// correlation matrix calculation
// [[Rcpp::export()]]
arma::mat corrMat(const arma::mat & Xmat){
  int p = Xmat.n_cols;
  int n = Xmat.n_rows;
  //
  arma::mat X = Xmat;
  arma::mat muvec(1,p);
  arma::mat sdvec(1,p);
  //
  for(int i=0;i<p;++i){
    arma::vec Xi = X.col(i);
    muvec(0,i) = arma::mean(Xi);
    sdvec(0,i) = arma::stddev(Xi);
  }
  arma::mat unos = ones<mat>(n,1);
  X -= unos*muvec;
  X *= arma::diagmat(1/sdvec);
  //
  arma::mat corMat = X.t()*X/n;
  return(corMat);
}

// Residual calculation given parent weights in B (wtd. adj. matrix)
// [[Rcpp::export()]]
arma::mat getResids(const arma::mat & X,const arma::mat & B){
  return( X-X*B );
}

// R-squared calculation given parent weights in B (wtd. adj. matrix)
// [[Rcpp::export()]]
arma::vec getRsq(const arma::mat & X,const arma::mat &B){
  unsigned int p = X.n_cols;
  unsigned int n = X.n_rows;
  arma::mat Resids = getResids(X,B);
  arma::vec rsq(p);
  for(unsigned int i =0; i<p; ++i){
    arma::vec resid = Resids.col(i);
    arma::vec residsq = arma::pow(resid,2.0);
    double sd = arma::stddev(X.col(i));
    rsq(i) = 1-arma::mean(residsq *(n)/(n-1) )/ (sd*sd);
  }
  return(rsq);
}