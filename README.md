# IMP
Importance sampling for aSPU


# IMPaSPU: R package for IMP without using Rcpp
*IMPaSPU_0.13.tar.gz*
See package manual for more description







# R code for IMP with Rcpp
*aimp.cpp*

Usage:
library('Rcpp')
library('inline')
library('RcppArmadillo')
sourceCpp("aimp.cpp")
aimp(stZ,D,B,Ga,wga,ddd,RR)
aimpW(stZ,D,B,Ga,wga,ddd,RR,div,enlarge)
Arguments:
item{stZ}{a numeric vector containing the Z-scores of p SNPs.}
item{D}{a p*p estimated correlation matrix of the p SNPs.}
item{B}{the number of iterations.}  
item{Ga}{a vector containing the power indexes for the SPU and aSPU tests (e.g. Ga=c(1,2,4,8,999)).}   
item{wga}{a vector containing the (initial) weights for the indexes (e.g. wga=c(0.2,0.2,0.2,0.2,0.2)).}    
item{ddd}{an abandoned parameter that has no effects on the results. Simply use ddd=1.}    
item{RR}{a B*p matrix of random variables (following i.i.d. N(0,1)). Will be used to simulate test statistics. Just generate RR=matrix(rnorm(B*p),nrow=B) beforehand, where p is the number of SNPs.}  
item{div}{the number of parts the iterations will be divided into (for aimpW). 5 is recommended.}  
item{enlarge}{how much the weight of the best index will be increased each time (for aimpW). 2 is recommended.}    
