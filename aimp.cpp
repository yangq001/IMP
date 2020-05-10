#define ARMA_NO_DEBUG

#include <armadillo>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <vector>

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec islarger( arma::vec x, double e) {
  int ncs = x.size();
  arma::vec out(ncs);
  for(int i = 0; i < ncs; i++){
    if(x[i]>e){
      out(i)=1;
    }else{
      out(i)=0;
    }
  }
  return(out);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec rep( double x, int e) {
  arma::vec x0(e);
  x0 = x0 * 0 + x;
  return x0;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec colSums(arma::mat X){
  int ncs = X.n_cols;
  arma::vec out(ncs);
  for(int i = 0; i < ncs; i++){
    out(i) = sum(X.col(i));
  }
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec rowSums(arma::mat X){
  int ncs = X.n_rows;
  arma::vec out(ncs);
  for(int i = 0; i < ncs; i++){
    out(i) = sum(X.row(i));
  }
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec rowMax(arma::mat X){
  int ncs = X.n_rows;
  arma::vec out(ncs);
  for(int i = 0; i < ncs; i++){
    out(i) = max(X.row(i));
  }
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sub1( arma::mat x, arma::uword e) {
  x.shed_col(e);
  x.shed_row(e);
  return x;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sub2( arma::mat x, arma::uword e) {
  //int nn=x.n_rows;
  //if (e==0){
  //  arma::mat shit3 = x.submat(1,0,nn-1,nn-1);
  //}else{
  //  if (e==x.n_rows){
  //    arma::mat shit3 = x.submat(0,0,nn-2,nn-1);
  //  }else{
  //    arma::mat shit3a = x.submat(0,0,e-1,nn-1);
  //    arma::mat shit3b = x.submat(e+1,0,nn-1,nn-1);
  //    arma::mat shit3 = join_vert(shit3a,shit3b);
  //  }
  //}
  //x.shed_col(e);
  x.shed_row(e);
  return x;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sub3( arma::mat x, arma::uword e) {
  x.shed_col(e);
  //x.shed_row(e);
  return x;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List aimp (arma::vec stZ, arma::mat& D,int B, arma::vec Ga, arma::vec wga, int ddd, arma::mat& RR) {

  const int n = D.n_rows;
  const int nga = Ga.size();
  
  arma::mat CovSsqrt = sqrtmat_sympd(D);

  //get miu's
  arma::vec mius(nga);
  
  for (int i = 0; i < nga; i++) {
    int ga = Ga[i];
  //    mius[i] = ga;
    if (ga < 100) {
      if(ga <= 2){
        mius[i] = pow(accu(pow(stZ,ga))/n,1.0/ga); 
      }else{
        mius[i] = pow(accu(pow(stZ,ga)),1.0/ga); 
      }
      //mius[i] = pow(abs(accu(pow(stZ,ga))),1.0/ga);    WRONG
    } else {
      arma::vec abstZ = abs(stZ);
      mius[i] = abstZ.max();
    }
  }
  
  
  //mius = mius * 0; //experiment 
  
  
  //conditional sqrt(cov)
  Rcpp::List conM;
  Rcpp::List conDsqrt;  

  for (int t = 0; t < n; t++) {
    arma::mat DD2 = sub2(D,t);
    
    conM[std::to_string(t)]=DD2.submat(0,t,n-2,t);
    
    arma::mat LL = sub1(D,t) - DD2.submat(0,t,n-2,t) * DD2.submat(0,t,n-2,t).t();
    conDsqrt[std::to_string(t)] = sqrtmat_sympd(LL);
  }  
  
  arma::vec wga0 = wga;
  
  //sample gamma   Ga   0:(nga-1)
  arma::vec aa = RcppArmadillo::sample(linspace(0,nga-1,nga), B, "T",wga0);    
  //arma::vec aa = RcppArmadillo::sample(Ga, B, "T");
  //sample index 0:(n-1)
  arma::vec tt = RcppArmadillo::sample(linspace(0,n-1,n), B, "T");  
  //sample sign 0:1
  arma::vec pm = RcppArmadillo::sample(linspace(0,1,2), B, "T");  
  pm = 2*pm-1;
  
  //simulate Z-scores
  arma::mat Z(B,n);
  
  for (int i = 0; i < B; i++) {
    int a = aa[i];
    int ga = Ga[a];
    
    if (ga <= 2) {
      //arma::vec dieplz = rnorm(n);
      //Z.submat(i,0,i,n-1) = rnorm(n);
      //arma::mat yo0 = RR.submat(i,0,i,n-1);
      Z.submat(i,0,i,n-1) = CovSsqrt * (RR.submat(i,0,i,n-1).t()) + mius[a] * pm[i];
      //Z.submat(i,0,i,n-1) = 0 * RR.submat(i,0,i,n-1);
    } else {
      Z.submat(i,0,i,n-1) = 0 * RR.submat(i,0,i,n-1);
      int t = tt[i];
      double yo1 = RR(i,t);
      //arma::mat yo2 = sub3(RR.submat(i,0,i,n-1),t); 
      
      //arma::mat yo0 = RR.submat(i,0,i,n-1);
      //arma::mat yo2 = sub3(yo0,t); 
      
      Z(i,t) = yo1 + mius[a] * pm[i];
      
      arma::mat shit3 = conM[std::to_string(t)];
      
      arma::mat miunew = Z.submat(i,0,i,n-2) * 0  + shit3.t() * Z(i,t);  //row
      
      arma::mat shit4 = conDsqrt[std::to_string(t)];      
      
      arma::mat yo3 = shit4 * (sub3(RR.submat(i,0,i,n-1),t).t()) + miunew.t();   //col
      //arma::mat yo3 = shit4 * yo2.t() + miunew.t();   //col
 
      //insert
      if(t==0){
        Z.submat(i,1,i,n-1) = yo3;
      }else{
        if(t==n-1){
          Z.submat(i,0,i,n-2) = yo3;
        }else{
          //arma::mat yo3a = yo3.submat(0,0,t-1,0);
          //arma::mat yo3b = yo3.submat(t,0,n-2,0);
          
          //Z.submat(i,0,i,t-1) = yo3a;
          //Z.submat(i,t+1,i,n-1) = yo3b;
          
          Z.submat(i,0,i,t-1) = yo3.submat(0,0,t-1,0).t();
          Z.submat(i,t+1,i,n-1) = yo3.submat(t,0,n-2,0).t();
        }
      }
    }
  }  
  
  //weights
  arma::mat soD = inv(D);
  //arma::mat x = Z.t();
  arma::mat xx = pow(Z,2);
  
  arma::mat ww(B,nga);
  ww = ww * 0;
  
  for (int j = 0; j < nga; j++) {
    int ga = Ga[j];    
    if(ga<=2){
      arma::vec miu = rep(mius[j],n);
      
      arma::mat mRm = miu.t() * soD * miu;  
      arma::mat mRZ = miu.t() * soD * Z.t();  
      arma::mat w1 = -0.5*(mRm(0,0)-2*mRZ);
      w1 = 0.5 * exp(w1);
      arma::mat w2 = -0.5*(mRm(0,0)+2*mRZ);
      w2 = 0.5 * exp(w2);
      arma::mat w3 = w1 + w2;
      
      ww.col(j) = w3;
      
    }else{
      double br = mius[j];
      
      ww.col(j) = 0.5*(colSums(exp( (xx.t() - pow(Z-br,2).t()) /2))+colSums(exp( (xx.t() - pow(Z+br,2).t()) /2)))/n;
    }
    
  }
  
  //arma::mat ww2=ww;
  
  arma::mat ww2=ww;
  
  ww=1/(ww * wga);
  
                  

  //p SPU
  arma::mat spu(B+1,nga);
  spu = spu*0;
  
  for (int j = 0; j < nga; j++) {
    int ga = Ga[j];    
    if(ga < 100){
      spu(0,j)=accu(pow(stZ,ga));
      spu.submat(1,j,B,j)=rowSums(pow(Z,ga));
    }else{
      spu(0,j)=max(abs(stZ));
      spu.submat(1,j,B,j)=rowMax(abs(Z));      
    }
  }
  
  arma::mat pspu(B+1,nga);  
  pspu = pspu*0;  
  
  //arma::mat oo(B+1,nga);
  for (int j = 0; j < nga; j++) {
    pspu(0,j)=accu(ww % islarger(abs(spu.submat(1,j,B,j)),abs(spu(0,j))))/B;
    
    arma::vec dieplz = -abs(spu.submat(1,j,B,j));
    
    //arma::vec hi(B);
    
    uvec hi = sort_index(dieplz);
    //add weights
    arma::vec dw(B);
    dw = dw * 0;
    pspu(hi[0]+1,j)=0;
    for (int k = 1; k < B; k++){
      dw[k] = dw[k-1] + ww[hi[k-1]];
      pspu(hi[k]+1,j)=dw[k]/(B-1);
    }
  }
  

  //p aSPU
  double aspu0 = min(pspu.row(0));
  arma::vec aspu(B);
  for (int i = 0; i < B; i++) {
    aspu[i] = min(pspu.row(i+1));
  }  
  
  double Paspu = accu(ww % islarger(-aspu,-aspu0))/B;
  
  
  //se
  arma::vec se(nga+1);
  for (int j = 0; j < nga; j++) {
    se[j]=accu( pow( ww % islarger(abs(spu.submat(1,j,B,j)),abs(spu(0,j))) - pspu(0,j), 2) )/B;
  }
  
  se[nga]=accu( pow( ww % islarger(-aspu,-aspu0) - Paspu, 2) )/B;
  
  se=pow(se/B,1.0/2);
  
  //s=sum((wei*(aspu[-1]<aspu[1])-paspu)^2)/(B-1)
  //se[length(Ga)+1]=sqrt(s/B)


  
  Rcpp::List res;  
  
  //res["wga"] = wga;
  
  //res["spu"] = spu;
  
  //res["pspu"] = pspu;

  res["Pspu"] = pspu.row(0);
  res["Paspu"] = Paspu;
  res["se"] = se;
  
  //res["Z"] = Z;
  
  //res["ww"] = ww;
  
  //res["ww2"] = ww2;
  
  //res["mius"] = mius;
  
  //res["aa"] = aa;
  
  //res["tt"] = tt;
  
  //res["pm"] = pm;

  return(res);
}









// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List aimpC0 (arma::vec stZ, arma::mat& D,int B, arma::vec Ga, arma::vec wga, int ddd, arma::mat& RR) {
  //SPU only
  
  const int n = D.n_rows;
  const int nga = Ga.size();
  
  arma::mat CovSsqrt = sqrtmat_sympd(D);
  
  //get miu's
  arma::vec mius(nga);
  
  for (int i = 0; i < nga; i++) {
    int ga = Ga[i];
    //    mius[i] = ga;
    if (ga < 100) {
      if(ga <= 2){
        mius[i] = pow(accu(pow(stZ,ga))/n,1.0/ga); 
      }else{
        mius[i] = pow(accu(pow(stZ,ga)),1.0/ga); 
      }
      //mius[i] = pow(abs(accu(pow(stZ,ga))),1.0/ga);    WRONG
    } else {
      arma::vec abstZ = abs(stZ);
      mius[i] = abstZ.max();
    }
  }
  
  
  //mius = mius * 0; //experiment 
  
  
  //conditional sqrt(cov)
  Rcpp::List conM;
  Rcpp::List conDsqrt;  
  
  for (int t = 0; t < n; t++) {
    arma::mat DD2 = sub2(D,t);
    
    conM[std::to_string(t)]=DD2.submat(0,t,n-2,t);
    
    arma::mat LL = sub1(D,t) - DD2.submat(0,t,n-2,t) * DD2.submat(0,t,n-2,t).t();
    conDsqrt[std::to_string(t)] = sqrtmat_sympd(LL);
  }  
  
  arma::vec wga0 = wga;
  
  //sample gamma   Ga   0:(nga-1)
  arma::vec aa = RcppArmadillo::sample(linspace(0,nga-1,nga), B, "T",wga0);  
  //arma::vec aa = RcppArmadillo::sample(Ga, B, "T");
  //sample index 0:(n-1)
  arma::vec tt = RcppArmadillo::sample(linspace(0,n-1,n), B, "T");  
  //sample sign 0:1
  arma::vec pm = RcppArmadillo::sample(linspace(0,1,2), B, "T");  
  pm = 2*pm-1;
  
  //simulate Z-scores
  arma::mat Z(B,n);
  
  for (int i = 0; i < B; i++) {
    int a = aa[i];
    int ga = Ga[a];
    
    if (ga <= 2) {
      Z.submat(i,0,i,n-1) = CovSsqrt * (RR.submat(i,0,i,n-1).t()) + mius[a] * pm[i];
      //Z.submat(i,0,i,n-1) = 0 * RR.submat(i,0,i,n-1);
    } else {
      Z.submat(i,0,i,n-1) = 0 * RR.submat(i,0,i,n-1);
      int t = tt[i];
      double yo1 = RR(i,t);

      Z(i,t) = yo1 + mius[a] * pm[i];
      
      arma::mat shit3 = conM[std::to_string(t)];
      
      arma::mat miunew = Z.submat(i,0,i,n-2) * 0  + shit3.t() * Z(i,t);  //row
      
      arma::mat shit4 = conDsqrt[std::to_string(t)];      
      
      arma::mat yo3 = shit4 * (sub3(RR.submat(i,0,i,n-1),t).t()) + miunew.t();   //col
      //arma::mat yo3 = shit4 * yo2.t() + miunew.t();   //col
      
      //insert
      if(t==0){
        Z.submat(i,1,i,n-1) = yo3;
      }else{
        if(t==n-1){
          Z.submat(i,0,i,n-2) = yo3;
        }else{
          Z.submat(i,0,i,t-1) = yo3.submat(0,0,t-1,0).t();
          Z.submat(i,t+1,i,n-1) = yo3.submat(t,0,n-2,0).t();
        }
      }
    }
  }  
  
  //weights
  arma::mat soD = inv(D);
  //arma::mat x = Z.t();
  arma::mat xx = pow(Z,2);
  
  arma::mat ww(B,nga);
  ww = ww * 0;
  
  for (int j = 0; j < nga; j++) {
    int ga = Ga[j];    
    if(ga<=2){
      arma::vec miu = rep(mius[j],n);
      
      arma::mat mRm = miu.t() * soD * miu;  
      arma::mat mRZ = miu.t() * soD * Z.t();  
      arma::mat w1 = -0.5*(mRm(0,0)-2*mRZ);
      w1 = 0.5 * exp(w1);
      arma::mat w2 = -0.5*(mRm(0,0)+2*mRZ);
      w2 = 0.5 * exp(w2);
      arma::mat w3 = w1 + w2;
      
      ww.col(j) = w3;

    }else{
      double br = mius[j];

      ww.col(j) = 0.5*(colSums(exp( (xx.t() - pow(Z-br,2).t()) /2))+colSums(exp( (xx.t() - pow(Z+br,2).t()) /2)))/n;
    }
    
  }
  
  //arma::mat ww2=ww;
  
  ww=1/(ww * wga);
  
  
  
  //p SPU
  arma::mat spu(B+1,nga);
  spu = spu*0;
  
  for (int j = 0; j < nga; j++) {
    int ga = Ga[j];    
    if(ga < 100){
      spu(0,j)=accu(pow(stZ,ga));
      spu.submat(1,j,B,j)=rowSums(pow(Z,ga));
    }else{
      spu(0,j)=max(abs(stZ));
      spu.submat(1,j,B,j)=rowMax(abs(Z));      
    }
  }
  
  arma::mat pspu(B+1,nga);  
  pspu = pspu*0;  
  
  //arma::mat oo(B+1,nga);
  for (int j = 0; j < nga; j++) {
    pspu(0,j)=accu(ww % islarger(abs(spu.submat(1,j,B,j)),abs(spu(0,j))))/B;
  }
  //s=sum((wei*(aspu[-1]<aspu[1])-paspu)^2)/(B-1)
  //se[length(Ga)+1]=sqrt(s/B)
  
  
  
  Rcpp::List res;  
  
  res["ww"] = ww;
  res["spu"] = spu;
  res["pspu"] = pspu;
  
 
  
  return(res);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List aimpW (arma::vec stZ, arma::mat& D,int B, arma::vec Ga, arma::vec wga, int ddd, arma::mat& RR, int div, double enlarge) {
  //update weights
  
  const int n = D.n_rows;
  const int nga = Ga.size();
  
  arma::mat CovSsqrt = sqrtmat_sympd(D);
  
  arma::vec wga1 = wga;
  
  arma::mat ww(B,1);
  ww = ww * 0;
  
  arma::mat spu(B+1,nga);
  spu = spu*0;
  
  int B1 = B/div;
  
  for (int kw = 0; kw < div; kw++) {
    
    arma::mat RR1 = RR.submat((B/div)*kw, 0, (B/div)*(kw+1)-1, n-1);
    
    Rcpp::List AS3 = aimpC0(stZ,D,B1,Ga,wga1,ddd,RR1);
    
    arma::mat wwB = AS3["ww"];
    arma::mat spuB = AS3["spu"];
    arma::mat pspuB = AS3["pspu"];
    
    ww.submat((B/div)*kw, 0, (B/div)*(kw+1)-1, 0) = wwB;
    
    spu.row(0) = spuB.row(0);
    
    spu.submat((B/div)*kw+1, 0, (B/div)*(kw+1), nga-1) = spuB.submat(1,0,B/div,nga-1); 
    
    arma::vec nono = pspuB.row(0);
    
    uvec hi = sort_index(nono);
    
    wga1[hi[0]] = wga1[hi[0]]*enlarge;
    
    wga1 = wga1/accu(wga1);
  } 
  
  arma::mat pspu(B+1,nga);  
  pspu = pspu*0;  
  
  //arma::mat oo(B+1,nga);
  for (int j = 0; j < nga; j++) {
    pspu(0,j)=accu(ww % islarger(abs(spu.submat(1,j,B,j)),abs(spu(0,j))))/B;
    
    arma::vec dieplz = -abs(spu.submat(1,j,B,j));
    
    //arma::vec hi(B);
    
    uvec hi = sort_index(dieplz);
    //add weights
    arma::vec dw(B);
    dw = dw * 0;
    pspu(hi[0]+1,j)=0;
    for (int k = 1; k < B; k++){
      dw[k] = dw[k-1] + ww[hi[k-1]];
      pspu(hi[k]+1,j)=dw[k]/(B-1);
    }
  }
  
  
  //p aSPU
  double aspu0 = min(pspu.row(0));
  arma::vec aspu(B);
  for (int i = 0; i < B; i++) {
    aspu[i] = min(pspu.row(i+1));
  }  
  
  double Paspu = accu(ww % islarger(-aspu,-aspu0))/B;
  
  
  //se
  arma::vec se(nga+1);
  for (int j = 0; j < nga; j++) {
    se[j]=accu( pow( ww % islarger(abs(spu.submat(1,j,B,j)),abs(spu(0,j))) - pspu(0,j), 2) )/B;
  }
  
  se[nga]=accu( pow( ww % islarger(-aspu,-aspu0) - Paspu, 2) )/B;
  
  se=pow(se/B,1.0/2);
  
  //s=sum((wei*(aspu[-1]<aspu[1])-paspu)^2)/(B-1)
  //se[length(Ga)+1]=sqrt(s/B)
  
  
  
  Rcpp::List res;  
  
  //res["wga1"] = wga1;
  
  //res["spu"] = spu;
  
  //res["pspu"] = pspu;
  
  res["Pspu"] = pspu.row(0);
  res["Paspu"] = Paspu;
  res["se"] = se;

  
  //res["ww"] = ww;
  
  return(res);
}
