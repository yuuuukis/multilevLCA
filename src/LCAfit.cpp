#define DARMA_USE_CURRENT
#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "SafeFunctions.h"
#include "Utils.h"

using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
List LCA_fast_poly_includeall(arma::mat mY, arma::mat mDesign, arma::ivec ivFreq, int iK, arma::mat mU, arma::ivec ivItemcat, arma::uvec first_poly, arma::vec reord_user,  int maxIter = 1e3, double tol = 1e-8, int reord = 0){
  // mY is equal to the number of observed response patterns x iH
  // ivItempos tells the number of categories for each item 
  //
  int iNtot    = accu(ivFreq);
  int iN = mY.n_rows;
  int iH    = mY.n_cols;
  int iV    = ivItemcat.n_elem;
  int NT    = iN;
  int h,j,k,n,p,v;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iV);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPhi(iH,iK);
  arma::vec pg(iK);
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  int ifooDcat = 0;
  arma::vec LLKSeries(maxIter);
  // 
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    // step M
    for(k = 0; k < iK; k++){
      pg(k) = accu(mU.col(k)%ivFreq)/iNtot;
      for(h = 0; h < iH; h++){
        mPhi(h,k) = accu(mU.col(k)%ivFreq%mY.col(h)%mDesign.col(h))/accu(mU.col(k)%ivFreq%mDesign.col(h));
        mPhi(h,k) = probcheck(mPhi(h,k));
      }
    }
    
    pg = OmegaCheck(pg, iK);
    // step E
    for(n = 0; n < iN; n++){
      vHdY.zeros();
      for(k = 0; k < iK; k++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            if(mDesign(n,ifooDcat)==1){
              mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
            }
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mDesign(n,ifooDcat)==1){
                if(mY(n,ifooDcat) > 0.0){
                  mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
                }
              }
              ifooDcat += 1;
            }
          }
          vHdY(k) += mdY(n,k,v);
        }
      }
      vLLK(n) = MixtDensityScale(pg, vHdY, iK);
      for(k = 0;k < iK; k++){
        mU(n,k) = exp(log(pg(k)) + vHdY(k) - vLLK(n));
      }
    }
    LLKSeries(iter) = accu(vLLK%ivFreq);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::vec vPhisum = sum(mPhi.rows(first_poly)).t();
    arma::uvec order = sort_index(vPhisum,"descending");
    //
    int ireord_check = 0;
    for(k=0; k< iK; k++){
      if(reord_user(k) != k){
        ireord_check = 1;
      }
    }
    if(ireord_check == 1){
      arma::uvec order_foo = order;
      int ifoo_reord_user;
      for(k=0; k< iK; k++){
        ifoo_reord_user = reord_user(k);
        order_foo(k) = order(ifoo_reord_user);
      }
      order = order_foo;
    }
    //
    int ifoo = 0;
    arma::mat mPhi_sorted = mPhi;
    arma::vec pg_sorted = pg;
    arma::mat mU_sorted = mU;
    for(k=0; k< iK; k++){
      ifoo               = order(k);
      mPhi_sorted.col(k) = mPhi.col(ifoo);
      pg_sorted(k)       = pg(ifoo);
      mU_sorted.col(k)   = mU.col(ifoo);
    }
    mPhi = mPhi_sorted;
    pg = pg_sorted;
    mU = mU_sorted;
  }
  arma::vec alphafoo(iK);
  alphafoo = log(pg/pg(0));
  arma::vec alpha(iK-1);
  alpha = alphafoo.subvec(1,iK-1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);

  arma::mat gamma = zeros(nfreepar_res,iK);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(k = 0; k < iK; k++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          gamma(iRoll,k) = log(mPhi(iItemfoo,k)/(1-mPhi(iItemfoo,k)));
        } else{
          gamma.col(k).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(k).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,k));
        }
      }else{
        if(ivItemcat(v)==1){
          gamma(iRoll,k) = log(mPhi(0,k)/(1-mPhi(0,k)));
        } else{
          gamma.col(k).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(k).subvec(1,ivItemcat(v)-1)/mPhi(0,k));
        }
      }
      if(ivItemcat(v)==1){
        iRoll += 1;
      } else{
        iRoll += (ivItemcat(v)-1);
      }
    }
  }
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==1){
      ivItemcat(v) = 2;
    }
  }

  iItemfoo = 0;
  arma::vec parvec = join_cols(alpha,vectorise(gamma));

  double BIC, AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(iNtot)*(iK*nfreepar_res + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iK*nfreepar_res + (iK - 1));
  double Terr = accu(-pg%trunc_log(pg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = sum(sum(-mU%mlogU,1)%ivFreq)/iNtot;
  double R2entr = (Terr - Perr)/Terr;

  // classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }

  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j)%ivFreq);
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }

  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iNtot);


  // Computing the score
  arma::mat mPg_Score(iN,iK-1);
  for(k =1; k< iK; k++){
    mPg_Score.col(k-1) = (mU.col(k) - pg(k))%ivFreq;
  }
  arma::mat mGamma_Score=zeros(iN,nfreepar_res*iK);
  iRoll = 0;
  iItemfoo = 0;
  for(k = 0; k < iK; k++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mGamma_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2)%repmat(mU.col(k)%ivFreq,1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(k).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mGamma_Score.col(iRoll) = mDesign.col(iItemfoo)%mU.col(k)%ivFreq%(mY.col(iItemfoo) - mPhi(iItemfoo,k));
          iItemfoo += ivItemcat(v)-1;
        }
      }else{
        if(ivItemcat(v)>2){
          mGamma_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(1,ivItemcat(v)-1)%repmat(mU.col(k)%ivFreq,1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(k).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mGamma_Score.col(iRoll) = mDesign.col(0)%mU.col(k)%ivFreq%(mY.col(0) - mPhi(0,k));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);
    }
  }
  iItemfoo = 0;

  arma::mat mScore = join_rows(mPg_Score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs =  sqrt(Varmat.diag());
  
  List EMout;
  EMout["mU"]            = mU;
  EMout["mPhi"]          = mPhi;
  EMout["pg"]            = pg;
  EMout["alphas"]        = alpha;
  EMout["gamma"]         = gamma;
  EMout["mScore"] =mScore;
  EMout["Varmat"] =Varmat;
  EMout["SEs"] =SEs;
  EMout["LLKSeries"]     = LLKSeries;
  EMout["eps"]           = eps;
  EMout["iter"]          = iter;
  EMout["BIC"]           = BIC;
  EMout["AIC"]           = AIC;
  EMout["R2entr"]        = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["freq"]          = ivFreq;
  EMout["parvec"]    = parvec;
  
  return EMout;
}

// [[Rcpp::export]]
List LCAcov_poly_includeall(arma::mat mY, arma::mat mDesign, arma::mat mZ, int iK, arma::mat mPhi, arma::mat mBeta, arma::mat mStep1Var, arma::ivec ivItemcat, int fixed = 0, int maxIter = 1e3, double tol = 1e-8, double NRtol = 1e-6, int NRmaxit = 100){
  // mU must be n x iK
  // ivItemcat is the vector of number of categories for each item
  // fixed = 1 if measurement model parameters should not be updated
  // mStep1Var to be used only if measurement model parameters are kept fixed at the input value
  // sample = 0 if standard matrix formula for Var(Beta) is to be computed; sample = 1 for the Monte Carlo approach
  // mBeta is the iP+1 x iK-1 matrix of structural regression coefficients
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int iP    = mBeta.n_rows;
  int iV    = ivItemcat.n_elem;
  int h,j,k,n,p,v;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iV);
  arma::mat mLogDensY = zeros(NT,iK);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPg = ones(NT,iK);
  arma::mat mU = zeros(NT,iK);
  for(k = 0; k < (iK - 1); k++){
    mPg.col(k + 1) = exp(mZ*mBeta.col(k));
  }
  arma::vec sumPg = sum(mPg,1);
  for(n = 0; n < iN; n++){
    mPg.row(n) = mPg.row(n)/sumPg(n);
  }
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  int ifooDcat = 0;
  // Compute log densities
  if(fixed > 0){
    for(n = 0; n < iN; n++){
      for(k = 0; k < iK; k++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            if(mDesign(n,ifooDcat)==1){
              mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
            }
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mDesign(n,ifooDcat)==1){
                if(mY(n,ifooDcat) > 0.0){
                  mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
                }
              }
              ifooDcat += 1;
            }
          }
          mLogDensY(n,k) += mdY(n,k,v);
        }
      }
    }
  }
  
  // 
  arma::mat mbeta = mBeta.t();
  arma::mat mbeta_score(iN,(iK-1)*iP);
  arma::vec LLKSeries(maxIter);
  List NR_step;
  while(eps > tol && iter<maxIter){
    // 
    // compute log densities
    if(fixed == 0){
      mLogDensY.zeros();
      for(n = 0; n < iN; n++){
        for(k = 0; k < iK; k++){
          ifooDcat = 0;
          for(v = 0; v< iV;v++){
            if(ivItemcat(v)==2){
              if(mDesign(n,ifooDcat)==1){
                mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
              }
              ifooDcat += 1;
            }else{
              for(p = 0; p < ivItemcat(v); p++){
                if(mDesign(n,ifooDcat)==1){
                  if(mY(n,ifooDcat) > 0.0){
                    mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
                  }
                }
                ifooDcat += 1;
              }
            }
            mLogDensY(n,k) += mdY(n,k,v);
          }
        }
      }
    }
    // E step
    for(n = 0; n < iN; n++){
      vLLK(n) = MixtDensityScale(mPg.row(n).t(), mLogDensY.row(n).t(), iK);
      for(k = 0; k < iK; k++){
        mU(n,k) = exp(log(mPg(n,k)) + mLogDensY(n,k) - vLLK(n));
      } 
    }
    // step M
    if(fixed == 0){
      for(k = 0; k < iK; k++){
        for(h = 0; h < iH; h++){
          mPhi(h,k) = accu(mU.col(k)%mY.col(h)%mDesign.col(h))/accu(mU.col(k)%mDesign.col(h));
          mPhi(h,k) = probcheck(mPhi(h,k));
        }
      }
    }
    NR_step = NR_step_covIT(mZ, mbeta, mU, NRtol, NRmaxit);
    arma::mat w_i = NR_step["w_i"];
    mPg = w_i;
    arma::mat mbeta_next = NR_step["beta"];
    arma::mat mbeta_score_curr = NR_step["mSbeta"];
    mbeta_score = mbeta_score_curr;
    mbeta = mbeta_next;
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);
  
  arma::mat gamma = zeros(nfreepar_res,iK);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(k = 0; k < iK; k++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          gamma(iRoll,k) = log(mPhi(iItemfoo,k)/(1-mPhi(iItemfoo,k)));
        } else{
          gamma.col(k).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(k).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,k));
        }
      }else{
        if(ivItemcat(v)==1){
          gamma(iRoll,k) = log(mPhi(0,k)/(1-mPhi(0,k)));
        } else{
          gamma.col(k).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(k).subvec(1,ivItemcat(v)-1)/mPhi(0,k));
        }
      }
      if(ivItemcat(v)==1){
        iRoll += 1;
      } else{
        iRoll += (ivItemcat(v)-1);
      }
    }
  }
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==1){
      ivItemcat(v) = 2;
    }
  }
  iItemfoo = 0;
  
  // 
  // 
  // mbeta is iK-1 x iP+1
  // 
  arma::vec parvec = join_cols(vectorise(mbeta.t()),vectorise(gamma));
  
  double BIC,AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(NT)*1.0*(iK*nfreepar_res + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2.0*(iK*nfreepar_res + (iK - 1));
  // double Terr = accu(-mPg%log(mPg));
  arma::mat mlogU = trunc_log(mU);
  arma::vec vPg = mean(mU).t();
  double Terr = iN*accu(-vPg%log(vPg));
  mlogU.elem(find_nonfinite(mlogU) ).zeros();
  mlogU.elem(find_nan(mlogU) ).zeros();
  double Perr = accu(-mU%mlogU);
  double R2entr = (Terr - Perr)/Terr;
  
  /// classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j));
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iN);
  
  // 
  // Computing SEs
  // simultaneous estimation
  
  // Computing the score
  
  
  arma::mat mGamma_Score=zeros(iN,nfreepar_res*iK);
  
  iRoll = 0;
  iItemfoo = 0;
  for(k = 0; k < iK; k++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mGamma_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2)%repmat(mU.col(k),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(k).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mGamma_Score.col(iRoll) = mDesign.col(iItemfoo)%mU.col(k)%(mY.col(iItemfoo) - mPhi(iItemfoo,k));
          iItemfoo += ivItemcat(v)-1;
        }
      }else{
        if(ivItemcat(v)>2){
          mGamma_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(1,ivItemcat(v)-1)%repmat(mU.col(k),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(k).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mGamma_Score.col(iRoll) = mDesign.col(0)%mU.col(k)%(mY.col(0) - mPhi(0,k));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);
    }
  }
  iItemfoo = 0;
  
  
  
  // Expected information matrix 
  
  arma::mat mScore = join_rows(mbeta_score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  
  // Matrix formula
  int uncondLatpars   = iK - 1;
  int parsfree        = (iK - 1)*iP;
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + (iK*nfreepar_res)-1,uncondLatpars + (iK*nfreepar_res)-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16); 
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + (iK*nfreepar_res)-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixed == 1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  
  // if(sample==1){
  //   // this subroutine samples from the approximate sampling distribution of the first step estimator
  //   
  // }
  // 
  List EMout;
  EMout["mU"]        = mU;
  EMout["mPhi"]      = mPhi;
  EMout["mPg"]        = mPg;
  EMout["beta"]    = mbeta.t();
  EMout["gamma"]     = gamma;
  EMout["SEs_unc"]     = SEs_unc;
  EMout["SEs_cor"]     = SEs_cor;
  EMout["mQ"]          = mQ;
  EMout["mV2"]         = mV2;
  EMout["Varmat_unc"]     = Varmat;
  EMout["Varmat_cor"]     = mVar_corr;
  EMout["LLKSeries"] = LLKSeries;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIC"]       = BIC;
  EMout["AIC"]       = AIC;
  EMout["R2entr"]    = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["parvec"]    = parvec;
  
  return EMout;
}  

// [[Rcpp::export]]
List LCAcov_poly(arma::mat mY, arma::mat mZ, int iK, arma::mat mPhi, arma::mat mBeta, arma::mat mStep1Var, arma::ivec ivItemcat, int fixed = 0, int maxIter = 1e3, double tol = 1e-8, double NRtol = 1e-6, int NRmaxit = 100){
  // mU must be n x iK
  // ivItemcat is the vector of number of categories for each item
  // fixed = 1 if measurement model parameters should not be updated
  // mStep1Var to be used only if measurement model parameters are kept fixed at the input value
  // sample = 0 if standard matrix formula for Var(Beta) is to be computed; sample = 1 for the Monte Carlo approach
  // mBeta is the iP+1 x iK-1 matrix of structural regression coefficients
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int iP    = mBeta.n_rows;
  int iV    = ivItemcat.n_elem;
  int h,j,k,n,p,v;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iV);
  arma::mat mLogDensY = zeros(NT,iK);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPg = ones(NT,iK);
  arma::mat mU = zeros(NT,iK);
  for(k = 0; k < (iK - 1); k++){
    mPg.col(k + 1) = exp(mZ*mBeta.col(k));
  }
  arma::vec sumPg = sum(mPg,1);
  for(n = 0; n < iN; n++){
    mPg.row(n) = mPg.row(n)/sumPg(n);
  }
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  int ifooDcat = 0;
  // Compute log densities
  if(fixed > 0){
    for(n = 0; n < iN; n++){
      for(k = 0; k < iK; k++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mY(n,ifooDcat) > 0.0){
                mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
              }
              ifooDcat += 1;
            }
          }
          mLogDensY(n,k) += mdY(n,k,v);
        }
      }
    }
  }
  
  // 
  arma::mat mbeta = mBeta.t();
  arma::mat mbeta_score(iN,(iK-1)*iP);
  arma::vec LLKSeries(maxIter);
  List NR_step;
  while(eps > tol && iter<maxIter){
    // 
    // compute log densities
    if(fixed == 0){
      mLogDensY.zeros();
      for(n = 0; n < iN; n++){
        for(k = 0; k < iK; k++){
          ifooDcat = 0;
          for(v = 0; v< iV;v++){
            if(ivItemcat(v)==2){
              mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
              ifooDcat += 1;
            }else{
              for(p = 0; p < ivItemcat(v); p++){
                if(mY(n,ifooDcat) > 0.0){
                  mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
                }
                ifooDcat += 1;
              }
            }
            mLogDensY(n,k) += mdY(n,k,v);
          }
        }
      }
    }
    // E step
    for(n = 0; n < iN; n++){
      vLLK(n) = MixtDensityScale(mPg.row(n).t(), mLogDensY.row(n).t(), iK);
      for(k = 0; k < iK; k++){
        mU(n,k) = exp(log(mPg(n,k)) + mLogDensY(n,k) - vLLK(n));
      } 
    }
    // step M
    if(fixed == 0){
      for(k = 0; k < iK; k++){
        for(h = 0; h < iH; h++){
          mPhi(h,k) = accu(mU.col(k)%mY.col(h))/accu(mU.col(k));
          mPhi(h,k) = probcheck(mPhi(h,k));
        }
      }
    }
    NR_step = NR_step_covIT(mZ, mbeta, mU, NRtol, NRmaxit);
    arma::mat w_i = NR_step["w_i"];
    mPg = w_i;
    arma::mat mbeta_next = NR_step["beta"];
    arma::mat mbeta_score_curr = NR_step["mSbeta"];
    mbeta_score = mbeta_score_curr;
    mbeta = mbeta_next;
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);
  arma::mat gamma = zeros(nfreepar_res,iK);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(k = 0; k < iK; k++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          gamma(iRoll,k) = log(mPhi(iItemfoo,k)/(1-mPhi(iItemfoo,k)));
        } else{
          gamma.col(k).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(k).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,k));
        }
      }else{
        if(ivItemcat(v)==1){
          gamma(iRoll,k) = log(mPhi(0,k)/(1-mPhi(0,k)));
        } else{
          gamma.col(k).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(k).subvec(1,ivItemcat(v)-1)/mPhi(0,k));
        }
      }
      if(ivItemcat(v)==1){
        iRoll += 1;
      } else{
        iRoll += (ivItemcat(v)-1);
      }
    }
  }
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==1){
      ivItemcat(v) = 2;
    }
  }
  iItemfoo = 0;
  
  // 
  // 
  // mbeta is iK-1 x iP+1
  // 
  arma::vec parvec = join_cols(vectorise(mbeta.t()),vectorise(gamma));
  
  double BIC,AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(NT)*1.0*(iK*nfreepar_res + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2.0*(iK*nfreepar_res + (iK - 1));
  // double Terr = accu(-mPg%log(mPg));
  arma::mat mlogU = trunc_log(mU);
  arma::vec vPg = mean(mU).t();
  double Terr = iN*accu(-vPg%log(vPg));
  mlogU.elem(find_nonfinite(mlogU) ).zeros();
  mlogU.elem(find_nan(mlogU) ).zeros();
  double Perr = accu(-mU%mlogU);
  double R2entr = (Terr - Perr)/Terr;
  
  /// classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j));
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iN);
  
  // 
  // Computing SEs
  // simultaneous estimation
  
  // Computing the score
  
  
  arma::mat mGamma_Score=zeros(iN,nfreepar_res*iK);
  
  iRoll = 0;
  iItemfoo = 0;
  for(k = 0; k < iK; k++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mGamma_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) = 
            repmat(mU.col(k),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(k).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mGamma_Score.col(iRoll) = mU.col(k)%(mY.col(iItemfoo) - mPhi(iItemfoo,k));
          iItemfoo += ivItemcat(v)-1;
        }
        
      }else{
        if(ivItemcat(v)>2){
          mGamma_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) = 
            repmat(mU.col(k),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(k).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mGamma_Score.col(iRoll) = mU.col(k)%(mY.col(0) - mPhi(0,k));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);  
    }
  } 
  iItemfoo = 0;
  
  
  
  // Expected information matrix 
  
  arma::mat mScore = join_rows(mbeta_score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  
  // Matrix formula
  int uncondLatpars   = iK - 1;
  int parsfree        = (iK - 1)*iP;
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + (iK*nfreepar_res)-1,uncondLatpars + (iK*nfreepar_res)-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16); 
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + (iK*nfreepar_res)-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixed == 1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  
  // if(sample==1){
  //   // this subroutine samples from the approximate sampling distribution of the first step estimator
  //   
  // }
  // 
  List EMout;
  EMout["mU"]        = mU;
  EMout["mPhi"]      = mPhi;
  EMout["mPg"]        = mPg;
  EMout["beta"]    = mbeta.t();
  EMout["gamma"]     = gamma;
  EMout["SEs_unc"]     = SEs_unc;
  EMout["SEs_cor"]     = SEs_cor;
  EMout["mQ"]          = mQ;
  EMout["mV2"]         = mV2;
  EMout["Varmat_unc"]     = Varmat;
  EMout["Varmat_cor"]     = mVar_corr;
  EMout["LLKSeries"] = LLKSeries;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIC"]       = BIC;
  EMout["AIC"]       = AIC;
  EMout["R2entr"]    = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["parvec"]    = parvec;
  
  return EMout;
}  

//[[Rcpp::export]]
List LCA_fast_poly(arma::mat mY, arma::ivec ivFreq, int iK, arma::mat mU, arma::ivec ivItemcat, arma::uvec first_poly, arma::vec reord_user,  int maxIter = 1e3, double tol = 1e-8, int reord = 0){
  // mY is equal to the number of observed response patterns x iH
  // ivItempos tells the number of categories for each item 
  //
  int iNtot    = accu(ivFreq);
  int iN = mY.n_rows;
  int iH    = mY.n_cols;
  int iV    = ivItemcat.n_elem;
  int NT    = iN;
  int h,j,k,n,p,v;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iV);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPhi(iH,iK);
  arma::vec pg(iK);
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  int ifooDcat = 0;
  arma::vec LLKSeries(maxIter);
  // 
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    // step M
    for(k = 0; k < iK; k++){
      pg(k) = accu(mU.col(k)%ivFreq)/iNtot;
      for(h = 0; h < iH; h++){
        mPhi(h,k) = accu(mU.col(k)%ivFreq%mY.col(h))/accu(mU.col(k)%ivFreq);
        mPhi(h,k) = probcheck(mPhi(h,k));
      }
    }
    
    pg = OmegaCheck(pg, iK);
    // step E
    for(n = 0; n < iN; n++){
      vHdY.zeros();
      for(k = 0; k < iK; k++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mY(n,ifooDcat) > 0.0){
                mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
              }
              ifooDcat += 1;
            }
          }
          vHdY(k) += mdY(n,k,v);
        }
      }
      vLLK(n) = MixtDensityScale(pg, vHdY, iK);
      for(k = 0;k < iK; k++){
        mU(n,k) = exp(log(pg(k)) + vHdY(k) - vLLK(n));
      }
    }
    LLKSeries(iter) = accu(vLLK%ivFreq);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::vec vPhisum = sum(mPhi.rows(first_poly)).t();
    arma::uvec order = sort_index(vPhisum,"descending");
    //
    int ireord_check = 0;
    for(k=0; k< iK; k++){
      if(reord_user(k) != k){
        ireord_check = 1;
      }
    }
    if(ireord_check == 1){
      arma::uvec order_foo = order;
      int ifoo_reord_user;
      for(k=0; k< iK; k++){
        ifoo_reord_user = reord_user(k);
        order_foo(k) = order(ifoo_reord_user);
      }
      order = order_foo;
    }
    //
    int ifoo = 0;
    arma::mat mPhi_sorted = mPhi;
    arma::vec pg_sorted = pg;
    arma::mat mU_sorted = mU;
    for(k=0; k< iK; k++){
      ifoo               = order(k);
      mPhi_sorted.col(k) = mPhi.col(ifoo);
      pg_sorted(k)       = pg(ifoo);
      mU_sorted.col(k)   = mU.col(ifoo); 
    }
    mPhi = mPhi_sorted;
    pg = pg_sorted;
    mU = mU_sorted;
  }
  arma::vec alphafoo(iK);
  alphafoo = log(pg/pg(0));
  arma::vec alpha(iK-1);
  alpha = alphafoo.subvec(1,iK-1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);

  arma::mat gamma = zeros(nfreepar_res,iK);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(k = 0; k < iK; k++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          gamma(iRoll,k) = log(mPhi(iItemfoo,k)/(1-mPhi(iItemfoo,k)));
        } else{
          gamma.col(k).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(k).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,k));
        }
      }else{
        if(ivItemcat(v)==1){
          gamma(iRoll,k) = log(mPhi(0,k)/(1-mPhi(0,k)));
        } else{
          gamma.col(k).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(k).subvec(1,ivItemcat(v)-1)/mPhi(0,k));
        }
      }
      if(ivItemcat(v)==1){
        iRoll += 1;
      } else{
        iRoll += (ivItemcat(v)-1);
      }
    }
  }
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==1){
      ivItemcat(v) = 2;
    }
  }
  
  iItemfoo = 0;
  arma::vec parvec = join_cols(alpha,vectorise(gamma));

  double BIC, AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(iNtot)*(iK*nfreepar_res + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iK*nfreepar_res + (iK - 1));
  double Terr = accu(-pg%trunc_log(pg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = sum(sum(-mU%mlogU,1)%ivFreq)/iNtot;
  double R2entr = (Terr - Perr)/Terr;

  // classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }

  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j)%ivFreq);
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }

  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iNtot);


  // Computing the score
  arma::mat mPg_Score(iN,iK-1);
  for(k =1; k< iK; k++){
    mPg_Score.col(k-1) = (mU.col(k) - pg(k))%ivFreq;
  }
  arma::mat mGamma_Score=zeros(iN,nfreepar_res*iK);
  iRoll = 0;
  iItemfoo = 0;
  for(k = 0; k < iK; k++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
        mGamma_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
          repmat(mU.col(k)%ivFreq,1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(k).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mGamma_Score.col(iRoll) = mU.col(k)%ivFreq%(mY.col(iItemfoo) - mPhi(iItemfoo,k));
          iItemfoo += ivItemcat(v)-1;
        }

      }else{
        if(ivItemcat(v)>2){
        mGamma_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
        repmat(mU.col(k)%ivFreq,1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(k).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mGamma_Score.col(iRoll) = mU.col(k)%ivFreq%(mY.col(0) - mPhi(0,k));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);
    }
  }
  iItemfoo = 0;

  arma::mat mScore = join_rows(mPg_Score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs =  sqrt(Varmat.diag());
  
  List EMout;
  EMout["mU"]            = mU;
  EMout["mPhi"]          = mPhi;
  EMout["pg"]            = pg;
  EMout["alphas"]        = alpha;
  EMout["gamma"]         = gamma;
  EMout["mScore"] =mScore;
  EMout["Varmat"] =Varmat;
  EMout["SEs"] =SEs;
  EMout["LLKSeries"]     = LLKSeries;
  EMout["eps"]           = eps;
  EMout["iter"]          = iter;
  EMout["BIC"]           = BIC;
  EMout["AIC"]           = AIC;
  EMout["R2entr"]        = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["freq"]          = ivFreq;
  EMout["parvec"]    = parvec;
  
  return EMout;
}

//[[Rcpp::export]]
double LCA_LLK(arma::vec parvec, arma::mat mY, int iK){
  // mU must be n*Ti x iK
  // works only with dichotomous items
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int h;
  int k;
  int n;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iH);
  arma::vec vLLK = zeros(NT,1);
  arma::vec pg=ones(iK,1);
  pg.subvec(1,iK-1) = exp(parvec.subvec(0,iK-2));
  pg = pg/accu(pg);
  
  int foopar = iK-1;
  arma::mat mPhi = reshape(conv_to<mat>::from(parvec(span(foopar,foopar + iH*iK -1))),iH,iK);
  mPhi = exp(mPhi)/(1.0 + exp(mPhi));
  
  arma::vec vHdY(iK);
  
  for(n = 0; n < iN; n++){
    vHdY.zeros();
    for(k = 0; k < iK; k++){
      for(h=0; h< iH;h++){
        mdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
        vHdY(k) += mdY(n,k,h);
      }
    }
    vLLK(n) = MixtDensityScale(pg, vHdY, iK);
  }
  
  double log_like = accu(vLLK);
  
  return log_like;
}

//[[Rcpp::export]]
arma::vec LCA_LLK_j(arma::vec parvec, arma::mat mY, int iK){
  // mU must be n*Ti x iK
  // works only with dichotomous items
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int h;
  int k;
  int n;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iH);
  arma::vec vLLK = zeros(NT,1);
  arma::vec pg=ones(iK,1);
  pg.subvec(1,iK-1) = exp(parvec.subvec(0,iK-2));
  pg = pg/accu(pg);
  
  int foopar = iK-1;
  arma::mat mPhi = reshape(conv_to<mat>::from(parvec(span(foopar,foopar + iH*iK -1))),iH,iK);
  mPhi = exp(mPhi)/(1.0 + exp(mPhi));
  
  arma::vec vHdY(iK);
  
  for(n = 0; n < iN; n++){
    vHdY.zeros();
    for(k = 0; k < iK; k++){
      for(h=0; h< iH;h++){
        mdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
        vHdY(k) += mdY(n,k,h);
      }
    }
    vLLK(n) = MixtDensityScale(pg, vHdY, iK);
  }
  
  return vLLK;
}

// [[Rcpp::export]]
List LCAcov(arma::mat mY, arma::mat mZ, int iK, arma::mat mPhi, arma::mat mBeta, arma::mat mStep1Var, int fixed = 0, int maxIter = 1e3, double tol = 1e-8, double NRtol = 1e-6, int NRmaxit = 100){
  // mU must be n x iK
  //
  // fixed = 1 if measurement model parameters should not be updated
  // mStep1Var to be used only if measurement model parameters are kept fixed at the input value
  // sample = 0 if standard matrix formula for Var(Beta) is to be computed; sample = 1 for the Monte Carlo approach
  // mBeta is the iP+1 x iK-1 matrix of structural regression coefficients
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int iP    = mBeta.n_rows;
  int h,j,k,n;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iH);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPg = ones(NT,iK);
  arma::mat mU = zeros(NT,iK);
  for(k = 0; k < (iK - 1); k++){
    mPg.col(k + 1) = exp(mZ*mBeta.col(k));
  }
  arma::vec sumPg = sum(mPg,1);
  for(n = 0; n < iN; n++){
    mPg.row(n) = mPg.row(n)/sumPg(n);
  }
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  // arma::mat mbeta = zeros(iK-1,iP);
  arma::mat mbeta = mBeta.t();
  arma::mat mbeta_score(iN,(iK-1)*iP);
  arma::vec LLKSeries(maxIter);
  List NR_step;
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    
    // step E
    
    for(n = 0; n < iN; n++){
      vHdY.zeros();
      for(k = 0; k < iK; k++){
        for(h=0; h< iH;h++){
          mdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
          vHdY(k) += mdY(n,k,h);
        }
      }
      vLLK(n) = MixtDensityScale(mPg.row(n).t(), vHdY, iK);
      for(k = 0;k < iK; k++){
        mU(n,k) = exp(log(mPg(n,k)) + vHdY(k) - vLLK(n));
      }
    }
    
    // step M
    if(fixed==0){
      for(k = 0; k < iK; k++){
        // pg(k) = mean(mU.col(k));
        for(h = 0; h < iH; h++){
          mPhi(h,k) = accu(mU.col(k)%mY.col(h))/accu(mU.col(k));
          mPhi(h,k) = probcheck(mPhi(h,k));
        }
      }
    }
    if(fixed==1){
      for(k = 0; k < iK; k++){
        // pg(k) = mean(mU.col(k));
        for(h = 0; h < iH; h++){
          mPhi(h,k) = probcheck(mPhi(h,k));
        }
      }
    }
    NR_step = NR_step_covIT(mZ, mbeta, mU, NRtol, NRmaxit);
    arma::mat w_i = NR_step["w_i"];
    mPg = w_i;
    arma::mat mbeta_next = NR_step["beta"];
    arma::mat mbeta_score_curr = NR_step["mSbeta"];
    mbeta_score = mbeta_score_curr;
    mbeta = mbeta_next;
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::mat gamma(iH,iK);
  gamma = log(mPhi/(1.0 - mPhi));
  // 
  // mbeta is iK-1 x iP+1
  // 
  arma::vec parvec = join_cols(vectorise(mbeta.t()),vectorise(gamma));
  
  double BIC,AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(NT)*(iH*iK + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iH*iK + (iK - 1));
  // double Terr = accu(-mPg%log(mPg));
  arma::mat mlogU = trunc_log(mU);
  arma::vec vPg = mean(mU).t();
  double Terr = iN*accu(-vPg%log(vPg));
  mlogU.elem(find_nonfinite(mlogU) ).zeros();
  mlogU.elem(find_nan(mlogU) ).zeros();
  double Perr = accu(-mU%mlogU);
  double R2entr = (Terr - Perr)/Terr;
  
  /// classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j));
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iN);
  
  // 
  // Computing SEs
  // simultaneous estimation
  
  // Computing the score
  
  
  arma::mat mGamma_Score=zeros(iN,iH*iK);
  int iroll = 0;
  for(k = 0; k < iK; k++){
    for(h = 0; h < iH; h++){
      mGamma_Score.col(iroll) = (mU.col(k)%(mY.col(h) - mPhi(h,k)));
      iroll += 1;
    }
  }
  
  iroll=0;
  
  // Expected information matrix 
  
  arma::mat mScore = join_rows(mbeta_score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  
  // Matrix formula
  int uncondLatpars   = iK - 1;
  int parsfree        = (iK - 1)*iP;
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + (iK*iH)-1,uncondLatpars + (iK*iH)-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,3.308722e-24);
  // arma::mat mJmatInv  = pinv(mJmat); 
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + (iK*iH)-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  mVar_corr           = symmatu(mVar_corr);
  arma::vec SEs_cor =  SEs_unc;
  if(fixed == 1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  
  // if(sample==1){
  //   // this subroutine samples from the approximate sampling distribution of the first step estimator
  //   
  // }
  // 
  List EMout;
  EMout["mU"]        = mU;
  EMout["mPhi"]      = mPhi;
  EMout["mPg"]        = mPg;
  EMout["beta"]    = mbeta.t();
  EMout["gamma"]     = gamma;
  EMout["SEs_unc"]     = SEs_unc;
  EMout["SEs_cor"]     = SEs_cor;
  EMout["mQ"]          = mQ;
  EMout["mV2"]         = mV2;
  EMout["Varmat_unc"]     = Varmat;
  EMout["Varmat_cor"]     = mVar_corr;
  EMout["LLKSeries"] = LLKSeries;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIC"]       = BIC;
  EMout["AIC"]       = AIC;
  EMout["R2entr"]    = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["parvec"]    = parvec;
  
  return EMout;
}

//[[Rcpp::export]]
List LCA(arma::mat mY, int iK, arma::mat mU, arma::vec reord_user, int maxIter = 1e3, double tol = 1e-8, int reord = 0){
  // mU must be n*Ti x iK
  //
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int h,k,j,n;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iH);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPhi(iH,iK);
  arma::vec pg(iK);
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  arma::vec LLKSeries(maxIter);
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    // step M
    for(k = 0; k < iK; k++){
      pg(k) = mean(mU.col(k));
      for(h = 0; h < iH; h++){
        mPhi(h,k) = accu(mU.col(k)%mY.col(h))/accu(mU.col(k));
        mPhi(h,k) = probcheck(mPhi(h,k));
      }
    }
    
    pg = OmegaCheck(pg, iK);
    // step E
    
    for(n = 0; n < iN; n++){
      vHdY.zeros();
      for(k = 0; k < iK; k++){
        for(h=0; h< iH;h++){
          mdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
          vHdY(k) += mdY(n,k,h);
        }
      }
      vLLK(n) = MixtDensityScale(pg, vHdY, iK);
      for(k = 0;k < iK; k++){
        mU(n,k) = exp(log(pg(k)) + vHdY(k) - vLLK(n));
      }
    }
    
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::vec vPhisum = sum(mPhi).t();
    arma::uvec order = sort_index(vPhisum,"descending");
    //
    int ireord_check = 0;
    for(k=0; k< iK; k++){
      if(reord_user(k) != k){
        ireord_check = 1;
      }
    }
    if(ireord_check == 1){
      arma::uvec order_foo = order;
      int ifoo_reord_user;
      for(k=0; k< iK; k++){
        ifoo_reord_user = reord_user(k);
        order_foo(k) = order(ifoo_reord_user);
      }
      order = order_foo;
    }
    //
    int ifoo = 0;
    arma::mat mPhi_sorted = mPhi;
    arma::vec pg_sorted = pg;
    arma::mat mU_sorted = mU;
    for(k=0; k< iK; k++){
      ifoo               = order(k);
      mPhi_sorted.col(k) = mPhi.col(ifoo);
      pg_sorted(k)       = pg(ifoo);
      mU_sorted.col(k)   = mU.col(ifoo); 
    }
    mPhi = mPhi_sorted;
    pg = pg_sorted;
    mU = mU_sorted;
  }
  arma::vec alphafoo(iK);
  alphafoo = log(pg/pg(0));
  arma::vec alpha(iK-1);
  alpha = alphafoo.subvec(1,iK-1);
  arma::mat gamma(iH,iK);
  gamma = log(mPhi/(1.0 - mPhi));
  
  arma::vec parvec = join_cols(alpha,vectorise(gamma));
  
  double BIC, AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(NT)*(iH*iK + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iH*iK + (iK - 1));
  double Terr = accu(-pg%trunc_log(pg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = mean(sum(-mU%mlogU,1));
  double R2entr = (Terr - Perr)/Terr;
  
  // classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j));
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iN);
  
  
  List EMout;
  EMout["mU"]            = mU;
  EMout["mPhi"]          = mPhi;
  EMout["pg"]            = pg;
  EMout["alphas"]        = alpha;
  EMout["gamma"]         = gamma;
  EMout["LLKSeries"]     = LLKSeries;
  EMout["eps"]           = eps;
  EMout["iter"]          = iter;
  EMout["BIC"]           = BIC;
  EMout["AIC"]           = AIC;
  EMout["R2entr"]        = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["parvec"]    = parvec;
  
  return EMout;
}

//[[Rcpp::export]]
List LCA_fast(arma::mat mY, arma::ivec ivFreq, int iK, arma::mat mU, arma::vec reord_user, int maxIter = 1e3, double tol = 1e-8, int reord = 0){
  // mY is equal to the number of observed response patterns x iH
  //
  int iNtot    = accu(ivFreq);
  int iN = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int h,k,j,n;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iH);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPhi(iH,iK);
  arma::vec pg(iK);
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  arma::vec LLKSeries(maxIter);
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    // step M
    for(k = 0; k < iK; k++){
      pg(k) = accu(mU.col(k)%ivFreq)/iNtot;
      for(h = 0; h < iH; h++){
        mPhi(h,k) = accu(mU.col(k)%ivFreq%mY.col(h))/accu(mU.col(k)%ivFreq);
        mPhi(h,k) = probcheck(mPhi(h,k));
      }
    }
    
    pg = OmegaCheck(pg, iK);
    // step E
    
    for(n = 0; n < iN; n++){
      vHdY.zeros();
      for(k = 0; k < iK; k++){
        for(h=0; h< iH;h++){
          mdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
          vHdY(k) += mdY(n,k,h);
        }
      }
      vLLK(n) = MixtDensityScale(pg, vHdY, iK);
      for(k = 0;k < iK; k++){
        mU(n,k) = exp(log(pg(k)) + vHdY(k) - vLLK(n));
      }
    }
    
    LLKSeries(iter) = accu(vLLK%ivFreq);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::vec vPhisum = sum(mPhi).t();
    arma::uvec order = sort_index(vPhisum,"descending");
    //
    int ireord_check = 0;
    for(k=0; k< iK; k++){
      if(reord_user(k) != k){
        ireord_check = 1;
      }
    }
    if(ireord_check == 1){
      arma::uvec order_foo = order;
      int ifoo_reord_user;
      for(k=0; k< iK; k++){
        ifoo_reord_user = reord_user(k);
        order_foo(k) = order(ifoo_reord_user);
      }
      order = order_foo;
    }
    //
    int ifoo = 0;
    arma::mat mPhi_sorted = mPhi;
    arma::vec pg_sorted = pg;
    arma::mat mU_sorted = mU;
    for(k=0; k< iK; k++){
      ifoo               = order(k);
      mPhi_sorted.col(k) = mPhi.col(ifoo);
      pg_sorted(k)       = pg(ifoo);
      mU_sorted.col(k)   = mU.col(ifoo); 
    }
    mPhi = mPhi_sorted;
    pg = pg_sorted;
    mU = mU_sorted;
  }
  arma::vec alphafoo(iK);
  alphafoo = log(pg/pg(0));
  arma::vec alpha(iK-1);
  alpha = alphafoo.subvec(1,iK-1);
  arma::mat gamma(iH,iK);
  gamma = log(mPhi/(1.0 - mPhi));
  
  arma::vec parvec = join_cols(alpha,vectorise(gamma));
  
  double BIC, AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(iNtot)*(iH*iK + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iH*iK + (iK - 1));
  double Terr = accu(-pg%trunc_log(pg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = sum(sum(-mU%mlogU,1)%ivFreq)/iNtot;
  double R2entr = (Terr - Perr)/Terr;
  
  // classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j)%ivFreq);
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iNtot);
  
  
  // Computing the score
  
  arma::mat mPg_Score(iN,iK-1);
  for(k =1; k< iK; k++){
    mPg_Score.col(k-1) = (mU.col(k) - pg(k))%ivFreq;
  }
  
  arma::mat mGamma_Score=zeros(iN,iH*iK);
  int iroll = 0;
  for(k = 0; k < iK; k++){
    for(h = 0; h < iH; h++){
      mGamma_Score.col(iroll) = (mU.col(k)%(mY.col(h) - mPhi(h,k)))%ivFreq;
      iroll += 1;
    }
  }
  
  iroll=0;
  
  arma::mat mScore = join_rows(mPg_Score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs =  sqrt(Varmat.diag());
  
  List EMout;
  EMout["mU"]            = mU;
  EMout["mPhi"]          = mPhi;
  EMout["pg"]            = pg;
  EMout["alphas"]        = alpha;
  EMout["gamma"]         = gamma;
  EMout["mScore"] =mScore;
  EMout["Varmat"] =Varmat;
  EMout["SEs"] =SEs;
  EMout["LLKSeries"]     = LLKSeries;
  EMout["eps"]           = eps;
  EMout["iter"]          = iter;
  EMout["BIC"]           = BIC;
  EMout["AIC"]           = AIC;
  EMout["R2entr"]        = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["freq"]          = ivFreq;
  EMout["parvec"]    = parvec;
  
  return EMout;
}

