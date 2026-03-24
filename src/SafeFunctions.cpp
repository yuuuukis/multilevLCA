#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// const double ZEROFLOAT = 1e-300;

const double dLowerProb = 1e-5;
const double dLowestProb = 1e-20;
const double dUpperProb = 1 - dLowerProb;

const double dZeroBound = 1e-10;
// const double dLowerBound = -1e5;
const double dLowerBound = -1e2;


// double zerocheck(double value){
//   if(abs(value) <= ZEROFLOAT){
//     value = sign(value)*ZEROFLOAT;
//   }
//   return value;
// }

double LogOfSum(arma::vec vLog) {
  
  double c    = max(vLog);
  double dLogSum = c + log(sum(exp(vLog - c)));
  
  return dLogSum;
  
}

double probcheck(double prob){
  if(prob < dLowerProb){
    prob = dLowerProb;
  }
  if(prob > dUpperProb){
    prob = dUpperProb;
  }
  return prob;
}

arma::mat GammaCheck(arma::mat mGamma, int iJ) {
  
  int i;
  int j;
  
  for (i = 0; i < iJ; i++) {
    for (j = 0; j < iJ; j++) {
      if (mGamma(i, j) < dLowerProb) {
        mGamma(i, j) = dLowerProb;
      }
      if (mGamma(i, j) > dUpperProb) {
        mGamma(i, j) = dUpperProb;
      }
    }
  }
  
  for (i = 0; i < iJ; i++) {
    mGamma.row(i) = mGamma.row(i) / accu(mGamma.row(i));
  }
  
  return mGamma;
  
}

arma::vec OmegaCheck(arma::vec vOmega, int iK) {
  
  int i;
  for (i = 0; i < iK; i++) {
    
    if (vOmega(i) < dLowerProb) {
      vOmega(i) = dLowerProb;
    }
    if (vOmega(i) > dUpperProb) {
      vOmega(i) = dUpperProb;
    }
    
  }
  
  vOmega = vOmega / accu(vOmega);
  
  return vOmega;
}


arma::vec PostCheck(arma::vec vPost, int iK) {
  
  int i;
  for (i = 0; i < iK; i++) {
    
    if (vPost(i) < dLowestProb) {
      vPost(i) = dLowerProb;
    }
    if (vPost(i) > (1.0 - dLowestProb)) {
      vPost(i) = 1.0 - dLowestProb;
    }
    
  }
  
  vPost = vPost / accu(vPost);
  
  return vPost;
}

double LogDensityCheck(double dLogDensity) {
  
  if (std::isnan(dLogDensity)) {
    dLogDensity = dLowerBound;
  }
  
  if (dLogDensity < dLowerBound) {
    
    dLogDensity = dLowerBound;
    
  }
  
  return dLogDensity;
  
}

double DensityCheck(double dDensity) {
  
  if (std::isnan(dDensity)) {
    dDensity = dZeroBound;
  }
  
  if (dDensity < dZeroBound) {
    dDensity = dZeroBound;
  }
  
  return dDensity;
  
}


//[[Rcpp::export]]
double MixtDensityScale(arma::vec vOmega, arma::vec vD_log, int iJ){
  
  arma::vec wp_log = log(vOmega) + vD_log;
  
  double dK = max(wp_log );
  
  arma::vec wp_log_scaled = wp_log - dK;
  
  double dLK = 0;
  int j;
  for (j = 0; j < iJ; j++) {
    dLK += exp(wp_log_scaled(j));
  }
  
  double dLLK = dK + log(dLK);
  
  dLLK = LogDensityCheck(dLLK);
  
  return dLLK;
}




