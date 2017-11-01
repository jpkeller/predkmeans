//
//  mlogitcpp.cpp
//  
//
//  Created by Joshua Keller on 3/26/15.
//
//

#include <stdio.h>
//#include <Rcpp.h>
#include "R.h"
#include "Rmath.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


/* Original Version
// [[Rcpp::export]]
double loglikeCpp(arma::mat X, arma::mat b, arma::mat y)
{
    arma::mat eXb;
    arma::colvec d;

    X *= b;
    eXb=exp(X);
    double res;
    d = log(1+ sum(eXb,1));
    y %= X;
    res = accu(y) - accu(d);
    return res;
} */



// [[Rcpp::export]]
double loglikeCpp(arma::mat X, arma::mat b, arma::mat y, int n)
{
    arma::mat eta;
    arma::vec denomMXmax;
    arma::mat denomMX;
    arma::vec denom;
    double res;
    arma::mat zerovec = arma::zeros<arma::mat>(n, 1);
    
    eta = X*b;
    denomMX = join_rows(zerovec, eta);
    denomMXmax = max(denomMX, 1);

    denomMX = denomMX.each_col() - denomMXmax;
    denom = denomMXmax + log(sum(exp(denomMX), 1));
    y %= eta;
    res = accu(y) - accu(denom);
    return res;
}

// [[Rcpp::export]]
arma::mat gradientMultinomialCpp(arma::mat X, arma::mat b, arma::mat y, int k)
{
    arma::mat eta;
    arma::mat etaTemp;
    arma::mat prob;
    arma::mat res;
    
    /*
    eta = exp(X*b);
    prob = 1+ sum(eta, 1);
    eta.each_col() /=prob;
    eta = y - eta;
    res=X.t()*eta;
     */
    
    eta = X*b;
    prob = eta; //overwrite this later
    for (int j=0; j<(k-1); j++){
        etaTemp = exp(eta.each_col() - eta.col(j));
        prob.col(j) = exp(-eta.col(j)) + sum(etaTemp, 1);
    }
    prob = pow(prob, -1);
    prob = y - prob;
    res=X.t()*prob;
    return res;
}


// [[Rcpp::export]]
arma::mat hessianMultinomialCpp(arma::mat X, arma::mat b, arma::mat y, int p, int k)
{
    //arma::mat probA;
    //arma::mat prob;
    arma::mat eta;
    arma::mat etaTemp;
    arma::mat prob;

    arma::mat H(p*(k-1),p*(k-1));
    arma::mat Hn(p, p);
    arma::mat Xtemp;

    
    /*
    prob = exp(X*b);
    probA = 1+ sum(prob, 1);
    prob.each_col() /=probA;
    */
    eta = X*b;
    prob = eta; //overwrite this later
    for (int j=0; j<(k-1); j++){
        etaTemp = exp(eta.each_col() - eta.col(j));
        prob.col(j) = exp(-eta.col(j)) + sum(etaTemp, 1);
    }
    prob = pow(prob, -1);
    
    for (int i=0; i<(k-1); i++){
        for (int j=0; j<(k-1); j++){
            Xtemp=X;
            if(i==j){
                Xtemp.each_col() %=prob.col(i)%(1-prob.col(j));
                Hn= X.t()*Xtemp;
                //Hn = X.t()*(prob.col(i)%(1-prob.col(j))%X);
            } else {
                Xtemp.each_col() %=-prob.col(i)%prob.col(j);
                Hn=X.t()*Xtemp;
               // Hn=X.t()*(-prob.col(i)%prob.col(j)%X);
            }
            H( arma::span(i*p, (i+1)*p-1), arma::span(j*p, (j+1)*p-1) )=Hn;
        }
    }
    return H;
}




// [[Rcpp::export]]
arma::mat getUproxy(arma::mat X, arma::mat b, arma::mat y)
{
    arma::mat eta;
    arma::mat prob;
    arma::mat res;
    
    eta = exp(X*b);
    prob = 1+ sum(eta, 1);
    eta.each_col() /=prob;
    eta = y - eta;
    res=X.t()*eta;
    return res;
}


// For each row Xi of X, this computes
//      exp(-1/(2*sigma2) * (Xi - mu)^T(Xi - mu))
//
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec getExpMahalRcpp(arma::mat x, arma::rowvec center, double sigma2) {
    arma::vec out;
    x.each_row() -= center;
    out=sum(x % x,1);
    out=exp(-0.5/sigma2*out);
    return out;
}



// For each row Xi of X, this computes
//      exp(-1/(2*sigma2) * (Xi - mu)^T(Xi - mu))
//
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat getHRcpp(arma::mat x, arma::mat R, arma::mat gamma, arma::mat mu, double sigma2) {
    arma::mat Rg;
    arma::colvec workingsum;
    int K=mu.n_rows;
    int n=x.n_rows;
    arma::mat out(n, K);

    for (int k=0; k<K; k++){
        out.col(k)= getExpMahalRcpp(x, mu.row(k), sigma2);
    }
    Rg = exp(R*gamma);
    workingsum = sum(Rg, 1);
    Rg.each_col() /=workingsum;
    out=out%Rg;
    workingsum=sum(out, 1);
    out.each_col()/=workingsum;
    return out;
}

