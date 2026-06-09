#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "SafeFunctions.h"
#include "Utils.h"

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
List MLTLCA_covWfixedlowhigh_poly_includeall(arma::mat mY, arma::mat mDesign, arma::mat mZ, arma::mat mZh, arma::vec vNj, arma::mat mDelta_start, arma::cube cGamma_start, arma::mat mPhi_start, arma::mat mStep1Var, arma::ivec ivItemcat, int maxIter = 1e3, double tol = 1e-8, int fixedpars = 0, int nsteps = 1, double NRtol = 1e-6, int NRmaxit = 100){
  // mY is iJ*sum(nj) x K
  // mZ is iJ*sum(nj) x iP, where iP is the number of iP-1 covariates (can even be 1) + an intercept term
  // mZh is iJ x iPh, where iPh is the number of iP-1 covariates + an intercept term
  // a column vector of ones to include the intercept
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // cPi is iN x iT x iM
  // cGamma is iT-1 x iP x iM
  // mDelta is (iM -1) x iPh
  // first Class is taken as reference
  // ivItemcat is the vector of number of categories for each item
  // fixedpars = 0 for one-step estimator
  // fixedpars = 1 for "pure" two-step (mPhi fixed, estimated using LCAfit or MLTLCA)
  // fixedpars = 2 for two-stage (fixing vOmega and mPhi, estimated using MLTLCA)
  
  int iV    = ivItemcat.n_elem;
  int n,j,k,m,t,p,v;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int iPh = mZh.n_cols;
  int iT = cGamma_start.n_rows + 1 ;
  int iM = cGamma_start.n_slices;
  arma::vec vOnes = ones(iJ,1);
  arma::cube cPi     = ones(iN,iT,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma_start.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  // 
  int iNparfoo = 1+(iM-1)+(iP-1);
  arma::mat mGamma_fixslope = zeros(iNparfoo,iT-1);
  mGamma_fixslope.row(0).fill(0.0);
  for(t = 0; t < (iT-1); t++){
    for(m = 1; m < iM; m++){
      mGamma_fixslope(m,t) = cGamma_start(t,0,m);
    }
    for(p = 0; p < (iP-1); p++){
      mGamma_fixslope.col(t).subvec(iM, iM + iP-2) = cGamma_start.slice(0).row(t).subvec(1,iP-1).t();
    }
  }
  // 
  arma::mat mOmega   = ones(iJ,iM);
  arma::mat mDelta_st_foo = mDelta_start.t();
  for(m = 1; m < iM; m++){
    mOmega.col(m) = exp(mZh*mDelta_st_foo.col(m-1));
  }
  for(j = 0; j < iJ; j++){
    mOmega.row(j) = mOmega.row(j)/accu(mOmega.row(j));
  }
  
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPXag    = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPXag   = mPXag;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::mat mPMX = zeros(iN,iT*iM);
  // 
  arma::mat mZrep = zeros(iM*iN,iNparfoo);
  arma::mat matPX = zeros(iM*iN,iT);
  arma::vec vPW_N = vectorise(mPW_N);
  // 
  int ifoopar = iN;
  mZrep.col(0).fill(1.0);
  for(m = 1; m < iM; m++){
    mZrep.col(m).subvec(ifoopar, ifoopar + iN - 1).fill(1.0);
    ifoopar += iN;
  }
  mZrep.cols(iM,iM +iP-2) = repmat(mZ.cols(1,iP-1),iM,1);
  // 
  arma::cube cLogdY = zeros(iN,iT,iV);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::mat mDelta    = mDelta_start;
  arma::mat mLogOmega = log(mOmega);
  arma::mat mDelta_Next = mDelta;
  
  arma::cube cGamma = cGamma_start;
  arma::mat mPhi = mPhi_start;
  arma::cube cGamma_Next = cGamma;
  arma::mat mGamma = zeros((iT-1)*iM,iP);
  arma::cube cGamma_Info(iP,iP,iM*(iT-1));
  arma::cube cDelta_Info(iPh,iPh,iM-1);
  arma::mat mDelta_Score(iN,iPh*(iM-1));
  arma::mat mGamma_Score(iN,(iT-1)*iP*iM);
  arma::mat mOmega_Next = mOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  int ifooDcat = 0;
  arma::vec foovec = zeros(iN,1);
  // 
  if(fixedpars>0){
    // compute log-densities
    for(n = 0; n < iN; n++){
      // 
      for(t = 0; t < iT; t++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            if(mDesign(n,ifooDcat)==1){
              cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
            }
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mDesign(n,ifooDcat)==1){
                if(mY(n,ifooDcat) > 0.0){
                  cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                }
              }
              ifooDcat += 1;
            }
          }
          mLogdKY(n,t) +=cLogdY(n,t,v);
        }
      }
    }
    // 
  }
  
  List NR_step;
  List NR_step_Delta;
  while(eps > tol && iter<maxIter){
    if(fixedpars == 0){
      // compute log-densities
      mLogdKY.zeros();
      for(n = 0; n < iN; n++){
        // 
        for(t = 0; t < iT; t++){
          ifooDcat = 0;
          for(v = 0; v< iV;v++){
            if(ivItemcat(v)==2){
              if(mDesign(n,ifooDcat)==1){
                cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              }
              ifooDcat += 1;
            }else{
              for(p = 0; p < ivItemcat(v); p++){
                if(mDesign(n,ifooDcat)==1){
                  if(mY(n,ifooDcat) > 0.0){
                    cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                  }
                }
                ifooDcat += 1;
              }
            }
            mLogdKY(n,t) +=cLogdY(n,t,v);
          }
        }
      }
    }
    //
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPXag(j,m) = accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
        mlogPW(j,m) = mLogOmega(j,m) + mlogPXag(j,m);
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1);
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX  = exp(clogPMX);
    mPXag = exp(mlogPXag);
    // M step 
    //
    NR_step = NR_step_covIT_wei(mZrep, mGamma_fixslope.t(), matPX, vPW_N, NRtol, NRmaxit);
    arma::mat mGamma_foo = NR_step["beta"];
    arma::mat mGammaScore_foo = NR_step["mSbeta"];
    arma::cube cGammaInfo_foo = NR_step["ibeta"];
    mGamma_Score = mGammaScore_foo;
    // cGamma is iT-1 x iP x iM 
    mGamma_fixslope = mGamma_foo.t();
    cGamma_Next.fill(0.0);
    cGamma_Next.slice(0).col(0) = mGamma_fixslope.row(0).t();
    for(m = 1; m < iM; m++){
      for(t = 0; t < (iT-1); t++){
        cGamma_Next(t,0,m) += mGamma_fixslope(m,t);
      }
    }
    for(m = 0; m < iM; m++){
      for(t = 0; t < (iT-1); t++){
        cGamma_Next.slice(m).row(t).subvec(1,iP-1) = mGamma_fixslope.col(t).subvec(iM,iNparfoo-1).t();
      }
    }
    
    for(m = 0; m < iM; m++){
      for(n = 0; n< iN; n++){
        for(t = 1; t < iT; t++){
          cPi_foo(n,t,m) = exp(accu(mZ.row(n)%cGamma_Next.slice(m).row(t-1)));
        }
        cPi_foo.slice(m).row(n) = cPi_foo.slice(m).row(n)/accu(cPi_foo.slice(m).row(n));
      }
    }
    // 
    if(fixedpars != 2){
      NR_step_Delta = NR_step_covIT(mZh, mDelta, mPW, NRtol, NRmaxit);
      // NR_step_Delta = NR_step_covIT_wei(mZh, mDelta,mPW, vOnes, NRtol, NRmaxit);
      arma::mat mDelta_foo = NR_step_Delta["beta"];
      arma::cube cDeltaInfo_foo = NR_step_Delta["ibeta"];
      arma::mat mDeltaScore_foo = NR_step_Delta["mSbeta"];
      arma::mat mOmega_foo = NR_step_Delta["w_i"];
      mOmega_Next = mOmega_foo;
      mDelta_Next = mDelta_foo;
      cDelta_Info = cDeltaInfo_foo;
      mDelta_Score = mDeltaScore_foo;
    }
    if(fixedpars == 0){
      mPhi_Next.zeros();
      mPMsumX = sum(cPMX,2);
      for(t = 0; t < iT; t++){
        for(k = 0; k < iK; k++){
          mPhi_Next(k,t) = accu(mPMsumX.col(t)%mY.col(k)%mDesign.col(k))/accu(mPMsumX.col(t)%mDesign.col(k));
          mPhi_Next(k,t) = probcheck(mPhi_Next(k,t));
        }
      }
    }
    if(fixedpars > 0){
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
          }
          mPMsumX.col(t) = foovec;
        }
      }
    }
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    
    cGamma    = cGamma_Next;
    mDelta    = mDelta_Next;
    cPi       = cPi_foo;
    cLogPi    = log(cPi);
    mOmega    = mOmega_Next;
    mLogOmega = log(mOmega);
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);
  
  double BIClow;
  double BIChigh;
  BIClow  = -2.0*LLKSeries(iter-1) + log(iN)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  
  double AIC;
  AIC = -2.0*LLKSeries(iter-1) + 2.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  
  // computing log-linear parameters
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  
  
  
  arma::mat mBeta = zeros(nfreepar_res,iT);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(t = 0; t < iT; t++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(iItemfoo,t)/(1-mPhi(iItemfoo,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,t));
        }
      }else{
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(0,t)/(1-mPhi(0,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(1,ivItemcat(v)-1)/mPhi(0,t));
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
  
  arma::vec parvec = join_cols(join_cols(vectorise(mDelta),vectorise(mGamma_fixslope)),vectorise(mBeta));
  
  arma::mat mBeta_Score=zeros(iN,nfreepar_res*iT);
  iRoll = 0;
  iItemfoo = 0;
  for(t = 0; t < iT; t++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2)%repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mDesign.col(iItemfoo)%mPMsumX.col(t)%(mY.col(iItemfoo) - mPhi(iItemfoo,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }else{
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(1,ivItemcat(v)-1)%repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(t).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mDesign.col(0)%mPMsumX.col(t)%(mY.col(0) - mPhi(0,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);
    }
  }
  iItemfoo = 0;
  
  arma::mat mDelta_Score_out(iN,(iM-1)*iPh);
  foo= 0;
  for(j = 0; j < iJ; j++){
    mDelta_Score_out.rows(foo,foo + vNj(j)-1)= repmat(mDelta_Score.row(j),vNj(j),1); 
    foo = foo + vNj(j);
  }
  arma::mat mScore = join_rows(repmat(mDelta_Score_out,iM,1),mGamma_Score);
  if(nsteps != 3){
    mScore = join_rows(mScore,repmat(mBeta_Score,iM,1));
  }
  if(nsteps == 3){
    nfreepar_res = 0;
  }
  arma::mat Infomat = mScore.t()*mScore/(iN*iM);
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/(iN*iM);
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  // asymptotic SEs correction
  int uncondLatpars   = (iM-1) + (iT-1)*iM;
  // int parsfree        = (iT - 1)*iP*iM + (iM - 1)*iPh;
  int parsfree        = 1 + (iM - 1) + (iT - 1)*iP + (iM - 1)*iPh;
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + nfreepar_res-1,uncondLatpars + nfreepar_res-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16); 
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + nfreepar_res-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixedpars==1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  
  if(nsteps==3){
    mV2.eye();
    mQ.eye();
    mVar_corr.eye();
    SEs_cor =  SEs_unc;
  }
  
  List EMout;
  EMout["mPhi"]        = mPhi;
  EMout["cGamma"]      = cGamma;
  EMout["mGamma_fixslope"] = mGamma_fixslope;
  EMout["mGamma_Score"] = mGamma_Score;
  EMout["mDelta"]      = mDelta;
  EMout["cDelta_Info"] = cDelta_Info;
  EMout["mDelta_Score"] = mDelta_Score;
  EMout["cPi"]       = cPi;
  EMout["mOmega"]    = mOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIChigh"]       = BIChigh;
  EMout["BIClow"]       = BIClow;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["AIC"]       = AIC;  
  EMout["mScore"]    = mScore;
  EMout["Infomat"]   = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs_unc"]   = SEs_unc;
  EMout["SEs_cor"]   = SEs_cor;
  EMout["mV2"]       = mV2;
  EMout["mQ"]        = mQ;
  EMout["mVar_corr"] = mVar_corr;
  
  return EMout;
  
}

// [[Rcpp::export]]
List MLTLCA_covWfixedlowhigh_poly(arma::mat mY, arma::mat mZ, arma::mat mZh, arma::vec vNj, arma::mat mDelta_start, arma::cube cGamma_start, arma::mat mPhi_start, arma::mat mStep1Var, arma::ivec ivItemcat, int maxIter = 1e3, double tol = 1e-8, int fixedpars = 0, int nsteps = 1, double NRtol = 1e-6, int NRmaxit = 100){
  // mY is iJ*sum(nj) x K
  // mZ is iJ*sum(nj) x iP, where iP is the number of iP-1 covariates (can even be 1) + an intercept term
  // mZh is iJ x iPh, where iPh is the number of iP-1 covariates + an intercept term
  // a column vector of ones to include the intercept
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // cPi is iN x iT x iM
  // cGamma is iT-1 x iP x iM
  // mDelta is (iM -1) x iPh
  // first Class is taken as reference
  // ivItemcat is the vector of number of categories for each item
  // fixedpars = 0 for one-step estimator
  // fixedpars = 1 for "pure" two-step (mPhi fixed, estimated using LCAfit or MLTLCA)
  // fixedpars = 2 for two-stage (fixing vOmega and mPhi, estimated using MLTLCA)
  
  int iV    = ivItemcat.n_elem;
  int n,j,k,m,t,p,v;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int iPh = mZh.n_cols;
  int iT = cGamma_start.n_rows + 1 ;
  int iM = cGamma_start.n_slices;
  arma::vec vOnes = ones(iJ,1);
  arma::cube cPi     = ones(iN,iT,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma_start.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  // 
  int iNparfoo = 1+(iM-1)+(iP-1);
  arma::mat mGamma_fixslope = zeros(iNparfoo,iT-1);
  mGamma_fixslope.row(0).fill(0.0);
  for(t = 0; t < (iT-1); t++){
    for(m = 1; m < iM; m++){
      mGamma_fixslope(m,t) = cGamma_start(t,0,m);
    }
    for(p = 0; p < (iP-1); p++){
      mGamma_fixslope.col(t).subvec(iM, iM + iP-2) = cGamma_start.slice(0).row(t).subvec(1,iP-1).t();
    }
  }
  // 
  arma::mat mOmega   = ones(iJ,iM);
  arma::mat mDelta_st_foo = mDelta_start.t();
  for(m = 1; m < iM; m++){
    mOmega.col(m) = exp(mZh*mDelta_st_foo.col(m-1));
  }
  for(j = 0; j < iJ; j++){
    mOmega.row(j) = mOmega.row(j)/accu(mOmega.row(j));
  }
  
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPXag    = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPXag   = mPXag;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::mat mPMX = zeros(iN,iT*iM);
  // 
  arma::mat mZrep = zeros(iM*iN,iNparfoo);
  arma::mat matPX = zeros(iM*iN,iT);
  arma::vec vPW_N = vectorise(mPW_N);
  // 
  int ifoopar = iN;
  mZrep.col(0).fill(1.0);
  for(m = 1; m < iM; m++){
    mZrep.col(m).subvec(ifoopar, ifoopar + iN - 1).fill(1.0);
    ifoopar += iN;
  }
  mZrep.cols(iM,iM +iP-2) = repmat(mZ.cols(1,iP-1),iM,1);
  //
  arma::cube cLogdY = zeros(iN,iT,iV);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::mat mDelta    = mDelta_start;
  arma::mat mLogOmega = log(mOmega);
  arma::mat mDelta_Next = mDelta;
  
  arma::cube cGamma = cGamma_start;
  arma::mat mPhi = mPhi_start;
  arma::cube cGamma_Next = cGamma;
  arma::mat mGamma = zeros((iT-1)*iM,iP);
  arma::cube cGamma_Info(iP,iP,iM*(iT-1));
  arma::cube cDelta_Info(iPh,iPh,iM-1);
  arma::mat mDelta_Score(iN,iPh*(iM-1));
  arma::mat mGamma_Score(iN,(iT-1)*iP*iM);
  arma::mat mOmega_Next = mOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  int ifooDcat = 0;
  arma::vec foovec = zeros(iN,1);
  // 
  if(fixedpars>0){
    // compute log-densities
    for(n = 0; n < iN; n++){
      // 
      for(t = 0; t < iT; t++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mY(n,ifooDcat) > 0.0){
                cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              }
              ifooDcat += 1;
            }
          }
          mLogdKY(n,t) +=cLogdY(n,t,v);
        }
      }
    }
    // 
  }
  
  List NR_step;
  List NR_step_Delta;
  while(eps > tol && iter<maxIter){
    if(fixedpars == 0){
      // compute log-densities
      mLogdKY.zeros();
      for(n = 0; n < iN; n++){
        // 
        for(t = 0; t < iT; t++){
          ifooDcat = 0;
          for(v = 0; v< iV;v++){
            if(ivItemcat(v)==2){
              cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              ifooDcat += 1;
            }else{
              for(p = 0; p < ivItemcat(v); p++){
                if(mY(n,ifooDcat) > 0.0){
                  cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                }
                ifooDcat += 1;
              }
            }
            mLogdKY(n,t) +=cLogdY(n,t,v);
          }
        }
      }
    }
    //
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPXag(j,m) = accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
        mlogPW(j,m) = mLogOmega(j,m) + mlogPXag(j,m);
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1);
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX  = exp(clogPMX);
    mPXag = exp(mlogPXag);
    // M step 
    //
    NR_step = NR_step_covIT_wei(mZrep, mGamma_fixslope.t(), matPX, vPW_N, NRtol, NRmaxit);
    arma::mat mGamma_foo = NR_step["beta"];
    arma::mat mGammaScore_foo = NR_step["mSbeta"];
    arma::cube cGammaInfo_foo = NR_step["ibeta"];
    
    mGamma_Score = mGammaScore_foo;
    // cGamma is iT-1 x iP x iM 
    mGamma_fixslope = mGamma_foo.t();
    cGamma_Next.fill(0.0);
    cGamma_Next.slice(0).col(0) = mGamma_fixslope.row(0).t();
    for(m = 1; m < iM; m++){
      for(t = 0; t < (iT-1); t++){
        cGamma_Next(t,0,m) += mGamma_fixslope(m,t);
      }
    }
    for(m = 0; m < iM; m++){
      for(t = 0; t < (iT-1); t++){
        cGamma_Next.slice(m).row(t).subvec(1,iP-1) = mGamma_fixslope.col(t).subvec(iM,iNparfoo-1).t();
      }
    }
    
    for(m = 0; m < iM; m++){
      for(n = 0; n< iN; n++){
        for(t = 1; t < iT; t++){
          cPi_foo(n,t,m) = exp(accu(mZ.row(n)%cGamma_Next.slice(m).row(t-1)));
        }
        cPi_foo.slice(m).row(n) = cPi_foo.slice(m).row(n)/accu(cPi_foo.slice(m).row(n));
      }
    }
    // 
    if(fixedpars != 2){
      NR_step_Delta = NR_step_covIT(mZh, mDelta, mPW, NRtol, NRmaxit);
      // NR_step_Delta = NR_step_covIT_wei(mZh, mDelta,mPW, vOnes, NRtol, NRmaxit);
      arma::mat mDelta_foo = NR_step_Delta["beta"];
      arma::cube cDeltaInfo_foo = NR_step_Delta["ibeta"];
      arma::mat mDeltaScore_foo = NR_step_Delta["mSbeta"];
      arma::mat mOmega_foo = NR_step_Delta["w_i"];
      mOmega_Next = mOmega_foo;
      mDelta_Next = mDelta_foo;
      cDelta_Info = cDeltaInfo_foo;
      mDelta_Score = mDeltaScore_foo;
    }
    if(fixedpars == 0){
      mPhi_Next.zeros();
      mPMsumX = sum(cPMX,2);
      for(t = 0; t < iT; t++){
        for(k = 0; k < iK; k++){
          mPhi_Next(k,t) = accu(mPMsumX.col(t)%mY.col(k))/accu(mPMsumX.col(t));
          mPhi_Next(k,t) = probcheck(mPhi_Next(k,t));
        }
      }
    }
    if(fixedpars > 0){
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
          }
          mPMsumX.col(t) = foovec;
        }
      }
    }
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    
    cGamma    = cGamma_Next;
    mDelta    = mDelta_Next;
    cPi       = cPi_foo;
    cLogPi    = log(cPi);
    mOmega    = mOmega_Next;
    mLogOmega = log(mOmega);
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);
  
  double BIClow;
  double BIChigh;
  BIClow  = -2.0*LLKSeries(iter-1) + log(iN)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  
  double AIC;
  AIC = -2.0*LLKSeries(iter-1) + 2.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  
  // computing log-linear parameters
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  
  
  
  arma::mat mBeta = zeros(nfreepar_res,iT);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(t = 0; t < iT; t++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(iItemfoo,t)/(1-mPhi(iItemfoo,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,t));
        }
      }else{
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(0,t)/(1-mPhi(0,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(1,ivItemcat(v)-1)/mPhi(0,t));
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
  
  arma::vec parvec = join_cols(join_cols(vectorise(mDelta),vectorise(mGamma_fixslope)),vectorise(mBeta));
  
  arma::mat mBeta_Score=zeros(iN,nfreepar_res*iT);
  iRoll = 0;
  iItemfoo = 0;
  for(t = 0; t < iT; t++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) = 
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(iItemfoo) - mPhi(iItemfoo,t));
          iItemfoo += ivItemcat(v)-1;
        }
        
      }else{
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) = 
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(t).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(0) - mPhi(0,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);  
    }
  } 
  iItemfoo = 0;
  
  arma::mat mDelta_Score_out(iN,(iM-1)*iPh);
  foo= 0;
  for(j = 0; j < iJ; j++){
    mDelta_Score_out.rows(foo,foo + vNj(j)-1)= repmat(mDelta_Score.row(j),vNj(j),1); 
    foo = foo + vNj(j);
  }
  arma::mat mScore = join_rows(repmat(mDelta_Score_out,iM,1),mGamma_Score);
  if(nsteps != 3){
    mScore = join_rows(mScore,repmat(mBeta_Score,iM,1));
  }
  if(nsteps == 3){
    nfreepar_res = 0;
  }
  arma::mat Infomat = mScore.t()*mScore/(iN*iM);
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/(iN*iM);
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  // asymptotic SEs correction
  int uncondLatpars   = (iM-1) + (iT-1)*iM;
  // int parsfree        = (iT - 1)*iP*iM + (iM - 1)*iPh;
  int parsfree        = 1 + (iM - 1) + (iT - 1)*iP + (iM - 1)*iPh;
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + nfreepar_res-1,uncondLatpars + nfreepar_res-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16); 
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + nfreepar_res-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixedpars==1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  
  if(nsteps==3){
    mV2.eye();
    mQ.eye();
    mVar_corr.eye();
    SEs_cor =  SEs_unc;
  }
  
  List EMout;
  EMout["mPhi"]        = mPhi;
  EMout["cGamma"]      = cGamma;
  EMout["mGamma_fixslope"] = mGamma_fixslope;
  EMout["mGamma_Score"] = mGamma_Score;
  EMout["mDelta"]      = mDelta;
  EMout["cDelta_Info"] = cDelta_Info;
  EMout["mDelta_Score"] = mDelta_Score;
  EMout["cPi"]       = cPi;
  EMout["mOmega"]    = mOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIChigh"]       = BIChigh;
  EMout["BIClow"]       = BIClow;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["AIC"]       = AIC;  
  EMout["mScore"]    = mScore;
  EMout["Infomat"]   = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs_unc"]   = SEs_unc;
  EMout["SEs_cor"]   = SEs_cor;
  EMout["mV2"]       = mV2;
  EMout["mQ"]        = mQ;
  EMout["mVar_corr"] = mVar_corr;
  
  return EMout;
  
}

// [[Rcpp::export]]
List MLTLCA_covWfixed_poly_includeall(arma::mat mY, arma::mat mDesign, arma::mat mZ, arma::vec vNj, arma::vec vOmega_start, arma::cube cGamma_start, arma::mat mPhi_start, arma::mat mStep1Var, arma::ivec ivItemcat, int maxIter = 1e3, double tol = 1e-8, int fixedpars = 0, int nsteps = 1, double NRtol = 1e-6, int NRmaxit = 100){
  // mY is iJ*sum(nj) x K
  // mZ is iJ*sum(nj) x iP, where iP is the number of iP-1 covariates (can even be 1) +
  // a column vector of ones to include the intercept
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // cPi is iN x iT x iM
  // cGamma is iT-1 x iP x iM but here we assume random intercepts and fixed slopes
  // first Class is taken as reference
  // ivItemcat is the vector of number of categories for each item
  // fixedpars = 0 for one-step estimator
  // fixedpars = 1 for "pure" two-step (mPhi fixed, estimated using LCAfit or MLTLCA)
  // fixedpars = 2 for two-stage (fixing vOmega and mPhi, estimated using MLTLCA)
  
  int iV    = ivItemcat.n_elem;
  int n,j,k,m,t,p,v;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int iT = cGamma_start.n_rows + 1 ;
  int iM = cGamma_start.n_slices;
  arma::cube cPi     = ones(iN,iT,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma_start.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  // 
  int iNparfoo = 1+(iM-1)+(iP-1);
  arma::mat mGamma_fixslope = zeros(iNparfoo,iT-1);
  mGamma_fixslope.row(0).fill(0.0);
  for(t = 0; t < (iT-1); t++){
    for(m = 1; m < iM; m++){
      mGamma_fixslope(m,t) = cGamma_start(t,0,m);
    }
    for(p = 0; p < (iP-1); p++){
      mGamma_fixslope.col(t).subvec(iM, iM + iP-2) = cGamma_start.slice(0).row(t).subvec(1,iP-1).t();
    }
  }
  //
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPXag    = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPXag   = mPXag;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::mat mPMX = zeros(iN,iT*iM);
  // 
  arma::mat mZrep = zeros(iM*iN,iNparfoo);
  arma::mat matPX = zeros(iM*iN,iT);
  arma::vec vPW_N = vectorise(mPW_N);
  // 
  int ifoopar = iN;
  mZrep.col(0).fill(1.0);
  for(m = 1; m < iM; m++){
    mZrep.col(m).subvec(ifoopar, ifoopar + iN - 1).fill(1.0);
    ifoopar += iN;
  }
  mZrep.cols(iM,iM +iP-2) = repmat(mZ.cols(1,iP-1),iM,1);
  // 
  arma::cube cLogdY = zeros(iN,iT,iV);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::vec vOmega = vOmega_start;
  arma::cube cGamma = cGamma_start;
  arma::mat mPhi = mPhi_start;
  arma::vec vLogOmega = log(vOmega);
  arma::cube cGamma_Next = cGamma;
  arma::mat mGamma = zeros((iT-1)*iM,iP);
  arma::mat mGamma_Score(iN*iM,1 + (iM-1) + (iT-1)*iP);
  arma::vec vOmega_Next = vOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  int ifooDcat = 0;
  arma::vec foovec = zeros(iN,1);
  //
  if(fixedpars>0){
    // compute log-densities
    for(n = 0; n < iN; n++){
      // 
      for(t = 0; t < iT; t++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(mDesign(n,v)==1){
            if(ivItemcat(v)==2){
              cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              ifooDcat += 1;
            }else{
              for(p = 0; p < ivItemcat(v); p++){
                if(mY(n,ifooDcat) > 0.0){
                  cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                }
                ifooDcat += 1;
              }
            }
            mLogdKY(n,t) +=cLogdY(n,t,v);
          }
        }
      }
    }
    // 
  }
  List NR_step;
  while(eps > tol && iter<maxIter){
    if(fixedpars == 0){
      // compute log-densities
      mLogdKY.zeros();
      for(n = 0; n < iN; n++){
        // 
        for(t = 0; t < iT; t++){
          ifooDcat = 0;
          for(v = 0; v< iV;v++){
            if(mDesign(n,v)==1){
              if(ivItemcat(v)==2){
                cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                ifooDcat += 1;
              }else{
                for(p = 0; p < ivItemcat(v); p++){
                  if(mY(n,ifooDcat) > 0.0){
                    cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                  }
                  ifooDcat += 1;
                }
              }
              mLogdKY(n,t) +=cLogdY(n,t,v);
            }
          }
        }
      }
    }
    //
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    ifoopar = 0;
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
        matPX.rows(ifoopar,ifoopar + iN-1) = cPX.slice(m);
        ifoopar = iN;
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPXag(j,m) = accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
        mlogPW(j,m) = vLogOmega(m) + mlogPXag(j,m);
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1);
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX  = exp(clogPMX);
    mPXag = exp(mlogPXag);
    // 
    vPW_N = vectorise(mPW_N);
    // 
    // M step 
    //
    NR_step = NR_step_covIT_wei(mZrep, mGamma_fixslope.t(), matPX, vPW_N, NRtol, NRmaxit);
    arma::mat mGamma_foo = NR_step["beta"];
    arma::mat mGammaScore_foo = NR_step["mSbeta"];
    arma::cube cGammaInfo_foo = NR_step["ibeta"];
    mGamma_Score = mGammaScore_foo;
    // cGamma is iT-1 x iP x iM 
    mGamma_fixslope = mGamma_foo.t();
    cGamma_Next.fill(0.0);
    cGamma_Next.slice(0).col(0) = mGamma_fixslope.row(0).t();
    for(m = 1; m < iM; m++){
      for(t = 0; t < (iT-1); t++){
        cGamma_Next(t,0,m) += mGamma_fixslope(m,t);
      }
    }
    for(m = 0; m < iM; m++){
      for(t = 0; t < (iT-1); t++){
        cGamma_Next.slice(m).row(t).subvec(1,iP-1) = mGamma_fixslope.col(t).subvec(iM,iNparfoo-1).t();
      }
    }
    
    for(m = 0; m < iM; m++){
      for(n = 0; n< iN; n++){
        for(t = 1; t < iT; t++){
          cPi_foo(n,t,m) = exp(accu(mZ.row(n)%cGamma_Next.slice(m).row(t-1)));
        }
        cPi_foo.slice(m).row(n) = cPi_foo.slice(m).row(n)/accu(cPi_foo.slice(m).row(n));
      }
    }
    
    if(fixedpars != 2){
      vOmega_Next = sum(mPW,0).t();
      vOmega_Next = vOmega_Next/accu(vOmega_Next);
    }
    if(fixedpars == 0){
      mPhi_Next.zeros();
      mPMsumX = sum(cPMX,2);
      for(t = 0; t < iT; t++){
        for(k = 0; k < iK; k++){
          mPhi_Next(k,t) = accu(mPMsumX.col(t)%mY.col(k)%mDesign.col(k))/accu(mPMsumX.col(t)%mDesign.col(k));
          mPhi_Next(k,t) = probcheck(mPhi_Next(k,t));
        }
      }
    }
    if(fixedpars > 0){
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
          }
          mPMsumX.col(t) = foovec;
        }
      }
    }
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    
    cGamma    = cGamma_Next;
    cPi       = cPi_foo;
    cLogPi    = log(cPi);
    vOmega    = vOmega_Next;
    vLogOmega = log(vOmega);
    
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  
  // computing log-linear parameters
  arma::vec vDeltafoo(iM);
  vDeltafoo = log(vOmega/vOmega(0));
  arma::vec vDelta(iM-1);
  vDelta = vDeltafoo.subvec(1,iM-1);
  int nfreepar_res = sum(ivItemcat_red);
  arma::mat mBeta = zeros(nfreepar_res,iT);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(t = 0; t < iT; t++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(iItemfoo,t)/(1-mPhi(iItemfoo,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,t));
        }
      }else{
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(0,t)/(1-mPhi(0,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(1,ivItemcat(v)-1)/mPhi(0,t));
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
  
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  
  arma::vec parvec = join_cols(join_cols(vDelta,vectorise(mGamma_fixslope)),vectorise(mBeta));
  
  double npars = parvec.n_elem;
  double BIClow;
  double BIChigh;
  double AIC;
  
  BIClow  = -2.0*LLKSeries(iter-1) + log(iN)*1.0*npars;
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*1.0*npars;
  AIC     = -2.0*LLKSeries(iter-1) + 2.0*npars;
  
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  
  arma::mat mOmega_Score(iN,iM-1);
  for(m =1; m< iM; m++){
    mOmega_Score.col(m-1) = mPW_N.col(m)*(1.0 - vOmega(m));
  }
  arma::mat mBeta_Score=zeros(iN,nfreepar_res*iT);
  iRoll = 0;
  iItemfoo = 0;
  for(t = 0; t < iT; t++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(iItemfoo) - mPhi(iItemfoo,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }else{
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(t).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(0) - mPhi(0,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);
    }
  }
  iItemfoo = 0;
  //
  
  // computing the score for cGamma
  arma::mat mScore = join_rows(repmat(mOmega_Score,iM,1),mGamma_Score);
  if(nsteps != 3){
    mScore = join_rows(mScore,repmat(mBeta_Score,iM,1));
  }
  
  arma::mat Infomat = mScore.t()*mScore/(iN*iM);
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/(iN*iM);
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  // asymptotic SEs correction
  int uncondLatpars   = (iM-1) + (iT-1)*iM;
  // int parsfree        = (iT - 1)*iP*iM + (iM - 1);
  int parsfree        = 1 + (iM - 1) + (iT - 1)*iP + (iM - 1);
  if(nsteps == 3){
    nfreepar_res = 0;
  }
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + nfreepar_res-1,uncondLatpars + nfreepar_res-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16);
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + nfreepar_res-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixedpars==1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  if(nsteps==3){
    mV2.eye();
    mQ.eye();
    mVar_corr.eye();
    SEs_cor =  SEs_unc;
  }
  
  List EMout;
  EMout["mPhi"]            = mPhi;
  EMout["cGamma"]          = cGamma;
  EMout["mGamma_fixslope"] = mGamma_fixslope;
  EMout["mGamma_Score"] = mGamma_Score;
  EMout["mOmega_Score"] = mOmega_Score;
  EMout["mBeta_Score"] = mBeta_Score;
  EMout["cPi"]       = cPi;
  EMout["vOmega"]    = vOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["vDelta"] = vDelta;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIChigh"]       = BIChigh;
  EMout["AIC"]       = AIC;
  EMout["BIClow"]       = BIClow;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  EMout["mScore"]    = mScore;
  EMout["Infomat"]   = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs_unc"]   = SEs_unc;
  EMout["SEs_cor"]   = SEs_cor;
  EMout["mV2"]       = mV2;
  EMout["mQ"]        = mQ;
  EMout["mVar_corr"] = mVar_corr;
  
  return EMout;
  
}

// [[Rcpp::export]]
List MLTLCA_covWfixed_poly(arma::mat mY, arma::mat mZ, arma::vec vNj, arma::vec vOmega_start, arma::cube cGamma_start, arma::mat mPhi_start, arma::mat mStep1Var, arma::ivec ivItemcat, int maxIter = 1e3, double tol = 1e-8, int fixedpars = 0, int nsteps = 1, double NRtol = 1e-6, int NRmaxit = 100){
  // mY is iJ*sum(nj) x K
  // mZ is iJ*sum(nj) x iP, where iP is the number of iP-1 covariates (can even be 1) +
  // a column vector of ones to include the intercept
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // cPi is iN x iT x iM
  // cGamma is iT-1 x iP x iM but here we assume random intercepts and fixed slopes
  // first Class is taken as reference
  // ivItemcat is the vector of number of categories for each item
  // fixedpars = 0 for one-step estimator
  // fixedpars = 1 for "pure" two-step (mPhi fixed, estimated using LCAfit or MLTLCA)
  // fixedpars = 2 for two-stage (fixing vOmega and mPhi, estimated using MLTLCA)
  
  int iV    = ivItemcat.n_elem;
  int n,j,k,m,t,p,v;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int iT = cGamma_start.n_rows + 1 ;
  int iM = cGamma_start.n_slices;
  arma::cube cPi     = ones(iN,iT,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma_start.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  // 
  int iNparfoo = 1+(iM-1)+(iP-1);
  arma::mat mGamma_fixslope = zeros(iNparfoo,iT-1);
  mGamma_fixslope.row(0).fill(0.0);
  for(t = 0; t < (iT-1); t++){
    for(m = 1; m < iM; m++){
      mGamma_fixslope(m,t) = cGamma_start(t,0,m);
    }
    for(p = 0; p < (iP-1); p++){
      mGamma_fixslope.col(t).subvec(iM, iM + iP-2) = cGamma_start.slice(0).row(t).subvec(1,iP-1).t();
    }
  }
  //
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPXag    = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPXag   = mPXag;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::mat mPMX = zeros(iN,iT*iM);
  // 
  arma::mat mZrep = zeros(iM*iN,iNparfoo);
  arma::mat matPX = zeros(iM*iN,iT);
  arma::vec vPW_N = vectorise(mPW_N);
  // 
  int ifoopar = iN;
  mZrep.col(0).fill(1.0);
  for(m = 1; m < iM; m++){
    mZrep.col(m).subvec(ifoopar, ifoopar + iN - 1).fill(1.0);
    ifoopar += iN;
  }
  mZrep.cols(iM,iM +iP-2) = repmat(mZ.cols(1,iP-1),iM,1);
  // 
  arma::cube cLogdY = zeros(iN,iT,iV);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::vec vOmega = vOmega_start;
  arma::cube cGamma = cGamma_start;
  arma::mat mPhi = mPhi_start;
  arma::vec vLogOmega = log(vOmega);
  arma::cube cGamma_Next = cGamma;
  arma::mat mGamma = zeros((iT-1)*iM,iP);
  arma::mat mGamma_Score(iN*iM,1 + (iM-1) + (iT-1)*iP);
  arma::vec vOmega_Next = vOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  int ifooDcat = 0;
  arma::vec foovec = zeros(iN,1);
  // 
  if(fixedpars>0){
    // compute log-densities
    for(n = 0; n < iN; n++){
      // 
      for(t = 0; t < iT; t++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mY(n,ifooDcat) > 0.0){
                cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              }
              ifooDcat += 1;
            }
          }
          mLogdKY(n,t) +=cLogdY(n,t,v);
        }
      }
    }
  }
  List NR_step;
  while(eps > tol && iter<maxIter){
    if(fixedpars==0){
      // compute log-densities
      mLogdKY.zeros();
      for(n = 0; n < iN; n++){
        // 
        for(t = 0; t < iT; t++){
          ifooDcat = 0;
          for(v = 0; v< iV;v++){
            if(ivItemcat(v)==2){
              cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              ifooDcat += 1;
            }else{
              for(p = 0; p < ivItemcat(v); p++){
                if(mY(n,ifooDcat) > 0.0){
                  cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                }
                ifooDcat += 1;
              }
            }
            mLogdKY(n,t) +=cLogdY(n,t,v);
          }
        }
      }
    }
    //
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    ifoopar = 0;
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
        matPX.rows(ifoopar,ifoopar + iN-1) = cPX.slice(m);
        ifoopar = iN;
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPXag(j,m) = accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
        mlogPW(j,m) = vLogOmega(m) + mlogPXag(j,m);
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1);
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX  = exp(clogPMX);
    mPXag = exp(mlogPXag);
    // 
    vPW_N = vectorise(mPW_N);
    // 
    // M step 
    //
    NR_step = NR_step_covIT_wei(mZrep, mGamma_fixslope.t(), matPX, vPW_N, NRtol, NRmaxit);
    arma::mat mGamma_foo = NR_step["beta"];
    arma::mat mGammaScore_foo = NR_step["mSbeta"];
    arma::cube cGammaInfo_foo = NR_step["ibeta"];
    
    mGamma_Score = mGammaScore_foo;
    // cGamma is iT-1 x iP x iM 
    mGamma_fixslope = mGamma_foo.t();
    cGamma_Next.fill(0.0);
    cGamma_Next.slice(0).col(0) = mGamma_fixslope.row(0).t();
    for(m = 1; m < iM; m++){
      for(t = 0; t < (iT-1); t++){
        cGamma_Next(t,0,m) += mGamma_fixslope(m,t);
      }
    }
    for(m = 0; m < iM; m++){
      for(t = 0; t < (iT-1); t++){
        cGamma_Next.slice(m).row(t).subvec(1,iP-1) = mGamma_fixslope.col(t).subvec(iM,iNparfoo-1).t();
      }
    }
    
    for(m = 0; m < iM; m++){
      for(n = 0; n< iN; n++){
        for(t = 1; t < iT; t++){
          cPi_foo(n,t,m) = exp(accu(mZ.row(n)%cGamma_Next.slice(m).row(t-1)));
        }
        cPi_foo.slice(m).row(n) = cPi_foo.slice(m).row(n)/accu(cPi_foo.slice(m).row(n));
      }
    }
    
    if(fixedpars != 2){
      vOmega_Next = sum(mPW,0).t();
      vOmega_Next = vOmega_Next/accu(vOmega_Next);
    }
    if(fixedpars == 0){
      mPhi_Next.zeros();
      mPMsumX = sum(cPMX,2);
      for(t = 0; t < iT; t++){
        for(k = 0; k < iK; k++){
          mPhi_Next(k,t) = accu(mPMsumX.col(t)%mY.col(k))/accu(mPMsumX.col(t));
          mPhi_Next(k,t) = probcheck(mPhi_Next(k,t));
        }
      }
    }
    if(fixedpars > 0){
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
          }
          mPMsumX.col(t) = foovec;
        }
      }
    }
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    
    cGamma    = cGamma_Next;
    cPi       = cPi_foo;
    cLogPi    = log(cPi);
    vOmega    = vOmega_Next;
    vLogOmega = log(vOmega);
    
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  
  // computing log-linear parameters
  arma::vec vDeltafoo(iM);
  vDeltafoo = log(vOmega/vOmega(0));
  arma::vec vDelta(iM-1);
  vDelta = vDeltafoo.subvec(1,iM-1);
  int nfreepar_res = sum(ivItemcat_red);
  arma::mat mBeta = zeros(nfreepar_res,iT);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(t = 0; t < iT; t++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(iItemfoo,t)/(1-mPhi(iItemfoo,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,t));
        }
      }else{
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(0,t)/(1-mPhi(0,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(1,ivItemcat(v)-1)/mPhi(0,t));
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
  
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  
  arma::vec parvec = join_cols(join_cols(vDelta,vectorise(mGamma_fixslope)),vectorise(mBeta));
  
  double npars = parvec.n_elem;
  double BIClow;
  double BIChigh;
  double AIC;
  
  BIClow  = -2.0*LLKSeries(iter-1) + log(iN)*1.0*npars;
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*1.0*npars;
  AIC     = -2.0*LLKSeries(iter-1) + 2.0*npars;
  
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  
  arma::mat mOmega_Score(iN,iM-1);
  for(m =1; m< iM; m++){
    mOmega_Score.col(m-1) = mPW_N.col(m)*(1.0 - vOmega(m));
  }
  arma::mat mBeta_Score=zeros(iN,nfreepar_res*iT);
  iRoll = 0;
  iItemfoo = 0;
  for(t = 0; t < iT; t++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(iItemfoo) - mPhi(iItemfoo,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }else{
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(t).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(0) - mPhi(0,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);
    }
  }
  iItemfoo = 0;
  //
  
  // computing the score for cGamma
  arma::mat mScore = join_rows(repmat(mOmega_Score,iM,1),mGamma_Score);
  if(nsteps != 3){
    mScore = join_rows(mScore,repmat(mBeta_Score,iM,1));
  }
  
  arma::mat Infomat = mScore.t()*mScore/(iN*iM);
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/(iN*iM);
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  // asymptotic SEs correction
  int uncondLatpars   = (iM-1) + (iT-1)*iM;
  // int parsfree        = (iT - 1)*iP*iM + (iM - 1);
  int parsfree        = 1 + (iM - 1) + (iT - 1)*iP + (iM - 1);
  if(nsteps == 3){
    nfreepar_res = 0;
  }
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + nfreepar_res-1,uncondLatpars + nfreepar_res-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16);
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + nfreepar_res-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixedpars==1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  if(nsteps==3){
    mV2.eye();
    mQ.eye();
    mVar_corr.eye();
    SEs_cor =  SEs_unc;
  }
  
  List EMout;
  EMout["mPhi"]            = mPhi;
  EMout["cGamma"]          = cGamma;
  EMout["mGamma_fixslope"] = mGamma_fixslope;
  EMout["mGamma_Score"] = mGamma_Score;
  EMout["mOmega_Score"] = mOmega_Score;
  EMout["mBeta_Score"] = mBeta_Score;
  EMout["cPi"]       = cPi;
  EMout["vOmega"]    = vOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["vDelta"] = vDelta;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIChigh"]       = BIChigh;
  EMout["AIC"]       = AIC;
  EMout["BIClow"]       = BIClow;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  EMout["mScore"]    = mScore;
  EMout["Infomat"]   = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs_unc"]   = SEs_unc;
  EMout["SEs_cor"]   = SEs_cor;
  EMout["mV2"]       = mV2;
  EMout["mQ"]        = mQ;
  EMout["mVar_corr"] = mVar_corr;
  
  return EMout;
  
}

// [[Rcpp::export]]
List MLTLCA_covlowhigh_poly_includeall(arma::mat mY, arma::mat mDesign, arma::mat mZ, arma::mat mZh, arma::vec vNj, arma::mat mDelta_start, arma::cube cGamma_start, arma::mat mPhi_start, arma::mat mStep1Var, arma::ivec ivItemcat, int maxIter = 1e3, double tol = 1e-8, int fixedpars = 0, int nsteps = 1, double NRtol = 1e-6, int NRmaxit = 100){
  // mY is iJ*sum(nj) x K
  // mZ is iJ*sum(nj) x iP, where iP is the number of iP-1 covariates (can even be 1) + an intercept term
  // mZh is iJ x iPh, where iPh is the number of iP-1 covariates + an intercept term
  // a column vector of ones to include the intercept
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // cPi is iN x iT x iM
  // cGamma is iT-1 x iP x iM
  // mDelta is (iM -1) x iPh
  // first Class is taken as reference
  // ivItemcat is the vector of number of categories for each item
  // fixedpars = 0 for one-step estimator
  // fixedpars = 1 for "pure" two-step (mPhi fixed, estimated using LCAfit or MLTLCA)
  // fixedpars = 2 for two-stage (fixing vOmega and mPhi, estimated using MLTLCA)
  
  int iV    = ivItemcat.n_elem;
  int n,j,k,m,t,p,v;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int iPh = mZh.n_cols;
  int iT = cGamma_start.n_rows + 1 ;
  int iM = cGamma_start.n_slices;
  arma::vec vOnes = ones(iJ,1);
  arma::cube cPi     = ones(iN,iT,iM);
  arma::mat mOmega   = ones(iJ,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma_start.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  arma::mat mDelta_st_foo = mDelta_start.t();
  for(m = 1; m < iM; m++){
    mOmega.col(m) = exp(mZh*mDelta_st_foo.col(m-1));
  }
  for(j = 0; j < iJ; j++){
    mOmega.row(j) = mOmega.row(j)/accu(mOmega.row(j));
  }
  
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPXag    = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPXag   = mPXag;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::mat mPMX = zeros(iN,iT*iM);
  arma::cube cLogdY = zeros(iN,iT,iV);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::mat mDelta    = mDelta_start;
  arma::mat mLogOmega = log(mOmega);
  arma::mat mDelta_Next = mDelta;
  
  arma::cube cGamma = cGamma_start;
  arma::mat mPhi = mPhi_start;
  arma::cube cGamma_Next = cGamma;
  arma::mat mGamma = zeros((iT-1)*iM,iP);
  arma::cube cGamma_Info(iP,iP,iM*(iT-1));
  arma::cube cDelta_Info(iPh,iPh,iM-1);
  arma::mat mDelta_Score(iN,iPh*(iM-1));
  arma::mat mGamma_Score(iN,(iT-1)*iP*iM);
  arma::mat mOmega_Next = mOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  int ifooDcat = 0;
  arma::vec foovec = zeros(iN,1);
  // 
  if(fixedpars>0){
    // compute log-densities
    for(n = 0; n < iN; n++){
      // 
      for(t = 0; t < iT; t++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            if(mDesign(n,ifooDcat)==1){
              cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
            }
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mDesign(n,ifooDcat)==1){
                if(mY(n,ifooDcat) > 0.0){
                  cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                }
              }
              ifooDcat += 1;
            }
          }
          mLogdKY(n,t) +=cLogdY(n,t,v);
        }
      }
    }
    // 
  }
  
  List NR_step;
  List NR_step_Delta;
  while(eps > tol && iter<maxIter){
    if(fixedpars == 0){
      // compute log-densities
      mLogdKY.zeros();
      for(n = 0; n < iN; n++){
        // 
        for(t = 0; t < iT; t++){
          ifooDcat = 0;
          for(v = 0; v< iV;v++){
            if(ivItemcat(v)==2){
              if(mDesign(n,ifooDcat)==1){
                cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              }
              ifooDcat += 1;
            }else{
              for(p = 0; p < ivItemcat(v); p++){
                if(mDesign(n,ifooDcat)==1){
                  if(mY(n,ifooDcat) > 0.0){
                    cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                  }
                }
                ifooDcat += 1;
              }
            }
            mLogdKY(n,t) +=cLogdY(n,t,v);
          }
        }
      }
    }
    //
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPXag(j,m) = accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
        mlogPW(j,m) = mLogOmega(j,m) + mlogPXag(j,m);
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1);
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX  = exp(clogPMX);
    mPXag = exp(mlogPXag);
    // M step 
    //
    int iFoo1 = 0;
    int iFoo2 = 0;
    for(m = 0; m < iM; m++){
      NR_step = NR_step_covIT_wei(mZ, cGamma.slice(m), cPX.slice(m), mPW_N.col(m), NRtol, NRmaxit);
      arma::mat mGamma_foo = NR_step["beta"];
      arma::cube cGammaInfo_foo = NR_step["ibeta"];
      arma::mat mGammaScore_foo = NR_step["mSbeta"];
      arma::mat mPi_foo = NR_step["w_i"];
      cPi_foo.slice(m) = mPi_foo;
      cGamma_Next.slice(m) = mGamma_foo;
      cGamma_Info.slices(iFoo1,iFoo1 + iT-2) = cGammaInfo_foo;
      mGamma_Score.cols(iFoo2, iFoo2 + ((iT-1)*iP) -1) = mGammaScore_foo;
      iFoo1 = iFoo1 + iT-1;
      iFoo2 = iFoo2 + ((iT-1)*iP);
    }
    // 
    if(fixedpars != 2){
      NR_step_Delta = NR_step_covIT(mZh, mDelta, mPW, NRtol, NRmaxit);
      // NR_step_Delta = NR_step_covIT_wei(mZh, mDelta,mPW, vOnes, NRtol, NRmaxit);
      arma::mat mDelta_foo = NR_step_Delta["beta"];
      arma::cube cDeltaInfo_foo = NR_step_Delta["ibeta"];
      arma::mat mDeltaScore_foo = NR_step_Delta["mSbeta"];
      arma::mat mOmega_foo = NR_step_Delta["w_i"];
      mOmega_Next = mOmega_foo;
      mDelta_Next = mDelta_foo;
      cDelta_Info = cDeltaInfo_foo;
      mDelta_Score = mDeltaScore_foo;
    }
    if(fixedpars == 0){
      mPhi_Next.zeros();
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t)%mDesign.col(k));
            foovec = foovec + (cPMX.slice(m).col(t)%mDesign.col(k));
            mPhi_Next(k,t) += accu(cPMX.slice(m).col(t)%mDesign.col(k)%mY.col(k));
          }
          mPhi_Next(k,t) = probcheck(mPhi_Next(k,t)/foo);
          mPMsumX.col(t) = foovec;
        }
      }
    }
    if(fixedpars > 0){
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
          }
          mPMsumX.col(t) = foovec;
        }
      }
    }
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    
    cGamma    = cGamma_Next;
    mDelta    = mDelta_Next;
    cPi = cPi_foo;
    cLogPi    = log(cPi);
    mOmega    = mOmega_Next;
    mLogOmega = log(mOmega);
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);
  
  double BIClow;
  double BIChigh;
  BIClow  = -2.0*LLKSeries(iter-1) + log(iN)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  
  double AIC;
  AIC = -2.0*LLKSeries(iter-1) + 2.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  
  // computing log-linear parameters
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  
  
  
  arma::mat mBeta = zeros(nfreepar_res,iT);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(t = 0; t < iT; t++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(iItemfoo,t)/(1-mPhi(iItemfoo,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,t));
        }
      }else{
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(0,t)/(1-mPhi(0,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(1,ivItemcat(v)-1)/mPhi(0,t));
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
  
  arma::vec parvec = join_cols(join_cols(vectorise(mDelta),vectorise(cGamma)),vectorise(mBeta));
  
  arma::mat mBeta_Score=zeros(iN,nfreepar_res*iT);
  iRoll = 0;
  iItemfoo = 0;
  for(t = 0; t < iT; t++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2)%repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mDesign.col(iItemfoo)%mPMsumX.col(t)%(mY.col(iItemfoo) - mPhi(iItemfoo,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }else{
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(1,ivItemcat(v)-1)%repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(t).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mDesign.col(0)%mPMsumX.col(t)%(mY.col(0) - mPhi(0,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);
    }
  }
  iItemfoo = 0;
  
  arma::mat mDelta_Score_out(iN,(iM-1)*iPh);
  foo= 0;
  for(j = 0; j < iJ; j++){
    mDelta_Score_out.rows(foo,foo + vNj(j)-1)= repmat(mDelta_Score.row(j),vNj(j),1); 
    foo = foo + vNj(j);
  }
  arma::mat mScore = join_rows(mDelta_Score_out,mGamma_Score);
  if(nsteps != 3){
    mScore = join_rows(mScore,mBeta_Score);
  }
  if(nsteps == 3){
    nfreepar_res = 0;
  }
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  // asymptotic SEs correction
  int uncondLatpars   = (iM-1) + (iT-1)*iM;
  int parsfree        = (iT - 1)*iP*iM + (iM - 1)*iPh;
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + nfreepar_res-1,uncondLatpars + nfreepar_res-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16); 
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + nfreepar_res-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixedpars==1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  
  if(nsteps==3){
    mV2.eye();
    mQ.eye();
    mVar_corr.eye();
    SEs_cor =  SEs_unc;
  }
  
  List EMout;
  EMout["mPhi"]        = mPhi;
  EMout["cGamma"]      = cGamma;
  EMout["cGamma_Info"] = cGamma_Info;
  EMout["mGamma_Score"] = mGamma_Score;
  EMout["mDelta"]      = mDelta;
  EMout["cDelta_Info"] = cDelta_Info;
  EMout["mDelta_Score"] = mDelta_Score;
  EMout["cPi"]       = cPi;
  EMout["mOmega"]    = mOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIChigh"]       = BIChigh;
  EMout["BIClow"]       = BIClow;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["AIC"]       = AIC;  
  EMout["mScore"]    = mScore;
  EMout["Infomat"]   = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs_unc"]   = SEs_unc;
  EMout["SEs_cor"]   = SEs_cor;
  EMout["mV2"]       = mV2;
  EMout["mQ"]        = mQ;
  EMout["mVar_corr"] = mVar_corr;
  
  return EMout;
  
}

// [[Rcpp::export]]
List MLTLCA_cov_poly_includeall(arma::mat mY, arma::mat mDesign, arma::mat mZ, arma::vec vNj, arma::vec vOmega_start, arma::cube cGamma_start, arma::mat mPhi_start, arma::mat mStep1Var, arma::ivec ivItemcat, int maxIter = 1e3, double tol = 1e-8, int fixedpars = 0, int nsteps = 1, double NRtol = 1e-6, int NRmaxit = 100){
  // mY is iJ*sum(nj) x K
  // mZ is iJ*sum(nj) x iP, where iP is the number of iP-1 covariates (can even be 1) +
  // a column vector of ones to include the intercept
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // cPi is iN x iT x iM
  // cGamma is iT-1 x iP x iM
  // first Class is taken as reference
  // ivItemcat is the vector of number of categories for each item
  // fixedpars = 0 for one-step estimator
  // fixedpars = 1 for "pure" two-step (mPhi fixed, estimated using LCAfit or MLTLCA)
  // fixedpars = 2 for two-stage (fixing vOmega and mPhi, estimated using MLTLCA)
  
  int iV    = ivItemcat.n_elem;
  int n,j,k,m,t,p,v;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int iT = cGamma_start.n_rows + 1 ;
  int iM = cGamma_start.n_slices;
  arma::cube cPi     = ones(iN,iT,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma_start.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPXag    = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPXag   = mPXag;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::mat mPMX = zeros(iN,iT*iM);
  arma::cube cLogdY = zeros(iN,iT,iV);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::vec vOmega = vOmega_start;
  arma::cube cGamma = cGamma_start;
  arma::mat mPhi = mPhi_start;
  arma::vec vLogOmega = log(vOmega);
  arma::cube cGamma_Next = cGamma;
  arma::mat mGamma = zeros((iT-1)*iM,iP);
  arma::cube cGamma_Info(iP,iP,iM*(iT-1));
  arma::mat mGamma_Score(iN,(iT-1)*iP*iM);
  arma::vec vOmega_Next = vOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  int ifooDcat = 0;
  arma::vec foovec = zeros(iN,1);
  // 
  if(fixedpars>0){
    // compute log-densities
    for(n = 0; n < iN; n++){
      // 
      for(t = 0; t < iT; t++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          // if(mDesign(n,v)==1){
          //   if(ivItemcat(v)==2){
          //     cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
          //     ifooDcat += 1;
          //   }else{
          //     for(p = 0; p < ivItemcat(v); p++){
          //       if(mY(n,ifooDcat) > 0.0){
          //         cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
          //       }
          //       ifooDcat += 1;
          //     }
          //   }
          //   mLogdKY(n,t) +=cLogdY(n,t,v);
          // }
          if(ivItemcat(v)==2){
            if(mDesign(n,ifooDcat)==1){
              cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
            }
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mDesign(n,ifooDcat)==1){
                if(mY(n,ifooDcat) > 0.0){
                  cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                }
              }
              ifooDcat += 1;
            }
          }
          mLogdKY(n,t) +=cLogdY(n,t,v);
        }
      }
    }
  }
  
  List NR_step;
  while(eps > tol && iter<maxIter){
    if(fixedpars==0){
      // compute log-densities
      mLogdKY.zeros();
      for(n = 0; n < iN; n++){
        // 
        for(t = 0; t < iT; t++){
          ifooDcat = 0;
          for(v = 0; v< iV;v++){
            if(ivItemcat(v)==2){
              if(mDesign(n,ifooDcat)==1){
                cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              }
              ifooDcat += 1;
            }else{
              for(p = 0; p < ivItemcat(v); p++){
                if(mDesign(n,ifooDcat)==1){
                  if(mY(n,ifooDcat) > 0.0){
                    cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                  }
                }
                ifooDcat += 1;
              }
            }
            mLogdKY(n,t) +=cLogdY(n,t,v);
          }
        }
      }
    }
    //
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPXag(j,m) = accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
        mlogPW(j,m) = vLogOmega(m) + mlogPXag(j,m);
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1);
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX  = exp(clogPMX);
    mPXag = exp(mlogPXag);
    // M step 
    //
    int iFoo1 = 0;
    int iFoo2 = 0;
    for(m = 0; m < iM; m++){
      NR_step = NR_step_covIT_wei(mZ, cGamma.slice(m), cPX.slice(m), mPW_N.col(m), NRtol, NRmaxit);
      arma::mat mGamma_foo = NR_step["beta"];
      arma::cube cGammaInfo_foo = NR_step["ibeta"];
      arma::mat mGammaScore_foo = NR_step["mSbeta"];
      arma::mat mPi_foo = NR_step["w_i"];
      cPi_foo.slice(m) = mPi_foo;
      cGamma_Next.slice(m) = mGamma_foo;
      cGamma_Info.slices(iFoo1,iFoo1 + iT-2) = cGammaInfo_foo;
      mGamma_Score.cols(iFoo2, iFoo2 + ((iT-1)*iP) -1) = mGammaScore_foo;
      iFoo1 = iFoo1 + iT-1;
      iFoo2 = iFoo2 + ((iT-1)*iP);
    }
    
    if(fixedpars != 2){
      vOmega_Next = sum(mPW,0).t();
      vOmega_Next = vOmega_Next/accu(vOmega_Next);
    }
    if(fixedpars == 0){
      mPhi_Next.zeros();
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t)%mDesign.col(k));
            foovec = foovec + (cPMX.slice(m).col(t)%mDesign.col(k));
            mPhi_Next(k,t) += accu(cPMX.slice(m).col(t)%mDesign.col(k)%mY.col(k));
          }
          mPhi_Next(k,t) = probcheck(mPhi_Next(k,t)/foo);
          mPMsumX.col(t) = foovec;
        }
      }
    }
    if(fixedpars > 0){
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
          }
          mPMsumX.col(t) = foovec;
        }
      }
    }
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    
    cGamma    = cGamma_Next;
    cPi = cPi_foo;
    cLogPi    = log(cPi);
    vOmega    = vOmega_Next;
    vLogOmega = log(vOmega);
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);
  double BIClow;
  double BIChigh;
  double AIC;
  BIClow  = -2.0*LLKSeries(iter-1) + log(iN)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0));
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0));
  AIC     = -2.0*LLKSeries(iter-1) + 2.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0));
  
  // computing log-linear parameters
  arma::vec vDeltafoo(iM);
  vDeltafoo = log(vOmega/vOmega(0));
  arma::vec vDelta(iM-1);
  vDelta = vDeltafoo.subvec(1,iM-1);
  
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  
  
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  arma::mat mBeta = zeros(nfreepar_res,iT);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(t = 0; t < iT; t++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(iItemfoo,t)/(1-mPhi(iItemfoo,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,t));
        }
      }else{
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(0,t)/(1-mPhi(0,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(1,ivItemcat(v)-1)/mPhi(0,t));
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
  
  arma::vec parvec = join_cols(join_cols(vDelta,vectorise(cGamma)),vectorise(mBeta));
  
  arma::mat mOmega_Score(iN,iM-1);
  for(m =1; m< iM; m++){
    mOmega_Score.col(m-1) = mPW_N.col(m)*(1.0 - vOmega(m));
  }
  arma::mat mBeta_Score=zeros(iN,nfreepar_res*iT);
  iRoll = 0;
  iItemfoo = 0;
  for(t = 0; t < iT; t++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2)%repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mDesign.col(iItemfoo)%mPMsumX.col(t)%(mY.col(iItemfoo) - mPhi(iItemfoo,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }else{
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(1,ivItemcat(v)-1)%repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(t).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mDesign.col(0)%mPMsumX.col(t)%(mY.col(0) - mPhi(0,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);
    }
  }
  iItemfoo = 0;
  //
  
  arma::mat mScore = join_rows(mOmega_Score,mGamma_Score);
  if(nsteps != 3){
    mScore = join_rows(mScore,mBeta_Score);
  }
  
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  // asymptotic SEs correction
  int uncondLatpars   = (iM-1) + (iT-1)*iM;
  int parsfree        = (iT - 1)*iP*iM + (iM - 1);
  if(nsteps == 3){
    nfreepar_res = 0;
  }
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + nfreepar_res-1,uncondLatpars + nfreepar_res-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16);
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + nfreepar_res-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixedpars==1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  if(nsteps==3){
    mV2.eye();
    mQ.eye();
    mVar_corr.eye();
    SEs_cor =  SEs_unc;
  }
  
  List EMout;
  EMout["mPhi"]        = mPhi;
  EMout["cGamma"]      = cGamma;
  EMout["cGamma_Info"] = cGamma_Info;
  EMout["mGamma_Score"] = mGamma_Score;
  EMout["mOmega_Score"] = mOmega_Score;
  EMout["mBeta_Score"] = mBeta_Score;
  EMout["cPi"]       = cPi;
  EMout["vOmega"]    = vOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["vDelta"] = vDelta;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIChigh"]       = BIChigh;
  EMout["AIC"]       = AIC;
  EMout["BIClow"]       = BIClow;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  EMout["mScore"]    = mScore;
  EMout["Infomat"]   = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs_unc"]   = SEs_unc;
  EMout["SEs_cor"]   = SEs_cor;
  EMout["mV2"]       = mV2;
  EMout["mQ"]        = mQ;
  EMout["mVar_corr"] = mVar_corr;
  
  return EMout;
  
}

// [[Rcpp::export]]
List MLTLCA_poly_includeall(arma::mat mY, arma::mat mDesign, arma::vec vNj, arma::vec vOmega, arma::mat mPi, arma::mat mPhi, arma::ivec ivItemcat, arma::uvec first_poly, arma::vec reord_user, arma::vec reord_user_high, int maxIter = 1e3, double tol = 1e-8, int reord = 1){
  // mY is iJ*sum(nj) x K
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // mPi is iT x iM
  // ivItemcat is the vector of number of categories for each item
  // mDesign should be iN x iV [number of responses non dichotomized]
  
  int iV    = ivItemcat.n_elem;
  int n,j,k,m,t,p,v;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = 1;
  
  int iT = mPi.n_rows;
  int iM = mPi.n_cols;
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::cube cLogdY = zeros(iN,iT,iV);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::mat mLogPi = log(mPi);
  arma::vec vLogOmega = log(vOmega);
  arma::mat mPi_Next = mPi;
  arma::vec vOmega_Next = vOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  int ifooDcat = 0;
  arma::vec foovec = zeros(iN,1);
  
  while(eps > tol && iter<maxIter){
    // compute log-densities
    mLogdKY.zeros();
    for(n = 0; n < iN; n++){
      // 
      for(t = 0; t < iT; t++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            if(mDesign(n,ifooDcat)==1){
              cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
            }
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mDesign(n,ifooDcat)==1){
                if(mY(n,ifooDcat) > 0.0){
                  cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                }
              }
              ifooDcat += 1;
            }
          }
          mLogdKY(n,t) +=cLogdY(n,t,v);
        }
      }
    }
    
    // 
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = mLogPi(t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(mPi.col(m),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPW(j,m) = vLogOmega(m) + accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1); 
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX = exp(clogPMX);
    // 
    // M step
    //
    mPhi_Next.zeros();
    mPi_Next = sum(cPMX,0);
    vOmega_Next = sum(mPW,0).t();
    vOmega_Next = vOmega_Next/accu(vOmega_Next);
    for(m = 0; m < iM; m++){
      mPi_Next.col(m) = mPi_Next.col(m)/accu(mPi_Next.col(m)); 
    }
    
    mPhi_Next.zeros();
    for(k = 0; k < iK; k++){
      for(t = 0; t < iT; t++){
        foo = 0.0;
        foovec.zeros();
        for(m = 0; m < iM; m++){
          foo += accu(cPMX.slice(m).col(t)%mDesign.col(k));
          foovec = foovec + (cPMX.slice(m).col(t)%mDesign.col(k));
          mPhi_Next(k,t) += accu(cPMX.slice(m).col(t)%mDesign.col(k)%mY.col(k));
        }
        mPhi_Next(k,t) = probcheck(mPhi_Next(k,t)/foo);
        mPMsumX.col(t) = foovec;
      }
    }
    for(m = 0; m < iM; m++){
      mPi_Next.col(m) = OmegaCheck(mPi_Next.col(m), iT);
    }
    vOmega_Next = OmegaCheck(vOmega_Next, iM);
    
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    mPi       = mPi_Next;
    vOmega    = vOmega_Next;
    mLogPi    = log(mPi);
    vLogOmega = log(vOmega);
    
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::cube cPhi_foo = zeros(iK,iT,iM);
    // int fp  = first_poly.n_elem;
    int fp  = iT;
    arma::mat mPhisum = zeros(fp,iM);
    for(k = 0; k < iK; k++){
      for(m = 0; m < iM; m++){ 
        for(t = 0; t < iT; t++){
          cPhi_foo(k,t,m) = accu(cPMX.slice(m).col(t)%mDesign.col(k)%mY.col(k))/accu(cPMX.slice(m).col(t)%mDesign.col(k));
        }
      }
    }
    
    
    arma::umat mLow_order(iT,iM);
    for(m = 0; m < iM; m++){
      mPhisum.col(m) = sum(cPhi_foo.slice(m).rows(first_poly)).t();
      mLow_order.col(m) = sort_index(mPhisum.col(m),"descending"); 
    }
    arma::vec vPhisum = sum(mPhi.rows(first_poly)).t();
    arma::uvec low_order = sort_index(vPhisum,"descending");
    arma::uvec high_order = sort_index(vOmega,"descending");
    int ireord_check = 0;
    for(t = 0; t < iT; t++){
      if(reord_user(t) != t){
        ireord_check = 1;
      }
    }
    if(ireord_check == 1){
      arma::umat mLow_order_foo = mLow_order;
      arma::uvec low_order_foo = low_order;
      int ifoo_reord_user;
      for(t = 0; t < iT; t++){
        ifoo_reord_user = reord_user(t);
        low_order_foo(t) = low_order(ifoo_reord_user);
        for(m = 0; m < iM; m++){
          mLow_order_foo(t,m) = mLow_order(ifoo_reord_user,m);
        }
      }
      low_order = low_order_foo;
      mLow_order = mLow_order_foo;
    }
    //
    int ireord_check_high = 0;
    for(m = 0; m < iM; m++){
      if(reord_user_high(m) != m){
        ireord_check_high = 1;
      }
    }
    if(ireord_check_high == 1){
      arma::umat mLow_order_foo = mLow_order;
      arma::uvec high_order_foo = high_order;
      int ifoo_reord_user_high;
      for(m = 0; m < iM; m++){
        ifoo_reord_user_high = reord_user_high(m);
        high_order_foo(m) = high_order(ifoo_reord_user_high);
        for(t = 0; t < iT; t++){
          mLow_order_foo(t,m) = mLow_order(t,ifoo_reord_user_high);
        }
      }
      high_order = high_order_foo;
      mLow_order = mLow_order_foo;
    }
    //
    int ifoo = 0;
    int ifoo_high = 0;
    arma::vec vOmega_sorted = vOmega;
    arma::mat mPhi_sorted = mPhi;
    arma::mat mPi_sorted = mPi;
    arma::mat mPW_sorted = mPW;
    arma::mat mPW_N_sorted = mPW_N;
    arma::cube cPX_sorted = cPX;
    arma::cube cPMX_sorted = cPMX;
    arma::cube clogPX_sorted = cPX;
    arma::cube clogPMX_sorted = cPMX;
    arma::mat mSumPX_sorted = mSumPX;
    arma::mat mPMsumX_sorted = mPMsumX;
    arma::cube cLogdY_sorted = cLogdY;
    arma::mat mLogdKY_sorted = mLogdKY;
    for(t=0; t< iT; t++){
      ifoo               = low_order(t);
      mPhi_sorted.col(t) = mPhi.col(ifoo);
      mPMsumX_sorted.col(t)  = mPMsumX.col(ifoo);
      mLogdKY_sorted.col(t)  = mLogdKY.col(ifoo);
      for(m = 0; m < iM; m++){
        ifoo =  mLow_order(t,m);
        mPi_sorted(t,m)    = mPi(ifoo,m);
        cPX_sorted.slice(m).col(t)  = cPX.slice(m).col(ifoo);
        clogPX_sorted.slice(m).col(t)  = clogPX.slice(m).col(ifoo);
        cPMX_sorted.slice(m).col(t)  = cPMX.slice(m).col(ifoo);
        clogPMX_sorted.slice(m).col(t)  = clogPMX.slice(m).col(ifoo);
      }
    }
    for(m = 0; m < iM; m++){
      ifoo_high = high_order(m);
      vOmega_sorted(m) = vOmega(ifoo_high);
      mPi_sorted.col(m)  = mPi.col(ifoo_high);
      cPX_sorted.slice(m)  = cPX.slice(ifoo_high);
      clogPX_sorted.slice(m)  = clogPX.slice(ifoo_high);
      cPMX_sorted.slice(m)  = cPMX.slice(ifoo_high);
      clogPMX_sorted.slice(m)  = clogPMX.slice(ifoo_high);
      mSumPX_sorted.col(m)  = mSumPX.col(ifoo_high);
    }
    vOmega = vOmega_sorted;
    mPhi = mPhi_sorted;
    mPi = mPi_sorted;
    cPX = cPX_sorted;
    clogPX = clogPX_sorted;
    cPMX = cPMX_sorted;
    clogPMX = clogPMX_sorted;
    mPMsumX = mPMsumX_sorted;
    mLogdKY = mLogdKY_sorted;
    mSumPX = mSumPX_sorted;
  }
  
  
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);
  
  double BIClow;
  double BIChigh;
  double AIC;
  BIClow = -2.0*LLKSeries(iter-1) + log(iN)*1.0*(iT*nfreepar_res + (iT - 1.0)*iM + (iM - 1.0));
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*1.0*(iT*nfreepar_res + (iT - 1.0)*iM + (iM - 1.0));
  AIC = -2.0*LLKSeries(iter-1) + 2.0*(iT*nfreepar_res + (iT - 1.0)*iM + (iM - 1.0));
  
  // computing log-linear parameters
  arma::vec vDeltafoo(iM);
  vDeltafoo = log(vOmega/vOmega(0));
  arma::vec vDelta = vDeltafoo.subvec(1,iM-1);
  arma::mat mGammafoo(iT,iM);
  arma::mat mGamma(iT-1,iM);
  for(m=0; m < iM; m++){
    mGammafoo.col(m) = log(mPi.col(m)/mPi(0,m));
    mGamma.col(m) = mGammafoo.col(m).subvec(1,iT-1);
  }
  //
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  
  
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  
  arma::mat mBeta = zeros(nfreepar_res,iT);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(t = 0; t < iT; t++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(iItemfoo,t)/(1-mPhi(iItemfoo,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,t));
        }
      }else{
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(0,t)/(1-mPhi(0,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(1,ivItemcat(v)-1)/mPhi(0,t));
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
  arma::vec parvec = join_cols(join_cols(vDelta,vectorise(mGamma)),vectorise(mBeta));
  
  // Computing the score
  
  arma::mat mOmega_Score(iN,iM-1);
  for(m =1; m< iM; m++){
    mOmega_Score.col(m-1) = mPW_N.col(m)*(1.0 - vOmega(m));
  }
  arma::mat mGamma_Score(iN,(iT-1)*iP*iM);
  int iFoo2 = 0;
  arma::mat mGammaScore_foo = zeros(iN,(iT-1)*iP);
  for(m = 0; m < iM; m++){
    mGammaScore_foo.fill(0.0);
    for(t = 1; t < iT; t++){
      mGammaScore_foo.col(t-1) = (cPX.slice(m).col(t) - mPi(t,m))%(mPW_N.col(m));
    }
    
    mGamma_Score.cols(iFoo2, iFoo2 + ((iT-1)*iP) -1) = mGammaScore_foo;
    iFoo2 = iFoo2 + ((iT-1)*iP);
  }
  iFoo2 = 0;
  
  arma::mat mBeta_Score=zeros(iN,nfreepar_res*iT);
  iRoll = 0;
  iItemfoo = 0;
  for(t = 0; t < iT; t++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2)%repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mDesign.col(iItemfoo)%mPMsumX.col(t)%(mY.col(iItemfoo) - mPhi(iItemfoo,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }else{
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            mDesign.cols(1,ivItemcat(v)-1)%repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(t).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mDesign.col(0)%mPMsumX.col(t)%(mY.col(0) - mPhi(0,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);
    }
  }
  iItemfoo = 0;
  
  
  arma::mat mScore = join_rows(join_rows(mOmega_Score,mGamma_Score),mBeta_Score);
  arma::mat Infomat = mScore.t()*mScore/(1.0*iN);
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/(1.0*iN);
  // arma::mat Varmat = inv(Infomat)/(1.0*iN);
  arma::vec SEs =  sqrt(Varmat.diag());
  
  List EMout;
  EMout["mScore"]    = mScore;
  EMout["Infomat"]    = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs"]    = SEs;
  EMout["mPhi"]      = mPhi;
  EMout["mPi"]       = mPi;
  EMout["vOmega"]    = vOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["mPXmarg"]   = mPXmarg;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["vDelta"] = vDelta;
  EMout["mGamma"] = mGamma;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIClow"]       = BIClow;
  EMout["BIChigh"]       = BIChigh;
  EMout["AIC"]       = AIC;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  
  return EMout;
}

// [[Rcpp::export]]
List MLTLCA_covlowhigh_poly(arma::mat mY, arma::mat mZ, arma::mat mZh, arma::vec vNj, arma::mat mDelta_start, arma::cube cGamma_start, arma::mat mPhi_start, arma::mat mStep1Var, arma::ivec ivItemcat, int maxIter = 1e3, double tol = 1e-8, int fixedpars = 0, int nsteps = 1, double NRtol = 1e-6, int NRmaxit = 100){
  // mY is iJ*sum(nj) x K
  // mZ is iJ*sum(nj) x iP, where iP is the number of iP-1 covariates (can even be 1) + an intercept term
  // mZh is iJ x iPh, where iPh is the number of iP-1 covariates + an intercept term
  // a column vector of ones to include the intercept
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // cPi is iN x iT x iM
  // cGamma is iT-1 x iP x iM
  // mDelta is (iM -1) x iPh
  // first Class is taken as reference
  // ivItemcat is the vector of number of categories for each item
  // fixedpars = 0 for one-step estimator
  // fixedpars = 1 for "pure" two-step (mPhi fixed, estimated using LCAfit or MLTLCA)
  // fixedpars = 2 for two-stage (fixing vOmega and mPhi, estimated using MLTLCA)
  
  int iV    = ivItemcat.n_elem;
  int n,j,k,m,t,p,v;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int iPh = mZh.n_cols;
  int iT = cGamma_start.n_rows + 1 ;
  int iM = cGamma_start.n_slices;
  arma::vec vOnes = ones(iJ,1);
  arma::cube cPi     = ones(iN,iT,iM);
  arma::mat mOmega   = ones(iJ,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma_start.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  arma::mat mDelta_st_foo = mDelta_start.t();
  for(m = 1; m < iM; m++){
    mOmega.col(m) = exp(mZh*mDelta_st_foo.col(m-1));
  }
  for(j = 0; j < iJ; j++){
    mOmega.row(j) = mOmega.row(j)/accu(mOmega.row(j));
  }
  
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPXag    = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPXag   = mPXag;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::mat mPMX = zeros(iN,iT*iM);
  arma::cube cLogdY = zeros(iN,iT,iV);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::mat mDelta    = mDelta_start;
  arma::mat mLogOmega = log(mOmega);
  arma::mat mDelta_Next = mDelta;
  
  arma::cube cGamma = cGamma_start;
  arma::mat mPhi = mPhi_start;
  arma::cube cGamma_Next = cGamma;
  arma::mat mGamma = zeros((iT-1)*iM,iP);
  arma::cube cGamma_Info(iP,iP,iM*(iT-1));
  arma::cube cDelta_Info(iPh,iPh,iM-1);
  arma::mat mDelta_Score(iN,iPh*(iM-1));
  arma::mat mGamma_Score(iN,(iT-1)*iP*iM);
  arma::mat mOmega_Next = mOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  int ifooDcat = 0;
  arma::vec foovec = zeros(iN,1);
  // 
  if(fixedpars>0){
    // compute log-densities
    for(n = 0; n < iN; n++){
      // 
      for(t = 0; t < iT; t++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mY(n,ifooDcat) > 0.0){
                cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              }
              ifooDcat += 1;
            }
          }
          mLogdKY(n,t) +=cLogdY(n,t,v);
        }
      }
    }
    // 
  }
  
  List NR_step;
  List NR_step_Delta;
  while(eps > tol && iter<maxIter){
    if(fixedpars == 0){
      // compute log-densities
      mLogdKY.zeros();
      for(n = 0; n < iN; n++){
        // 
        for(t = 0; t < iT; t++){
          ifooDcat = 0;
          for(v = 0; v< iV;v++){
            if(ivItemcat(v)==2){
              cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              ifooDcat += 1;
            }else{
              for(p = 0; p < ivItemcat(v); p++){
                if(mY(n,ifooDcat) > 0.0){
                  cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                }
                ifooDcat += 1;
              }
            }
            mLogdKY(n,t) +=cLogdY(n,t,v);
          }
        }
      }
    }
    //
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPXag(j,m) = accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
        mlogPW(j,m) = mLogOmega(j,m) + mlogPXag(j,m);
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1);
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX  = exp(clogPMX);
    mPXag = exp(mlogPXag);
    // M step 
    //
    int iFoo1 = 0;
    int iFoo2 = 0;
    for(m = 0; m < iM; m++){
      NR_step = NR_step_covIT_wei(mZ, cGamma.slice(m), cPX.slice(m), mPW_N.col(m), NRtol, NRmaxit);
      arma::mat mGamma_foo = NR_step["beta"];
      arma::cube cGammaInfo_foo = NR_step["ibeta"];
      arma::mat mGammaScore_foo = NR_step["mSbeta"];
      arma::mat mPi_foo = NR_step["w_i"];
      cPi_foo.slice(m) = mPi_foo;
      cGamma_Next.slice(m) = mGamma_foo;
      cGamma_Info.slices(iFoo1,iFoo1 + iT-2) = cGammaInfo_foo;
      mGamma_Score.cols(iFoo2, iFoo2 + ((iT-1)*iP) -1) = mGammaScore_foo;
      iFoo1 = iFoo1 + iT-1;
      iFoo2 = iFoo2 + ((iT-1)*iP);
    }
    // 
    if(fixedpars != 2){
      NR_step_Delta = NR_step_covIT(mZh, mDelta, mPW, NRtol, NRmaxit);
      // NR_step_Delta = NR_step_covIT_wei(mZh, mDelta,mPW, vOnes, NRtol, NRmaxit);
      arma::mat mDelta_foo = NR_step_Delta["beta"];
      arma::cube cDeltaInfo_foo = NR_step_Delta["ibeta"];
      arma::mat mDeltaScore_foo = NR_step_Delta["mSbeta"];
      arma::mat mOmega_foo = NR_step_Delta["w_i"];
      mOmega_Next = mOmega_foo;
      mDelta_Next = mDelta_foo;
      cDelta_Info = cDeltaInfo_foo;
      mDelta_Score = mDeltaScore_foo;
    }
    if(fixedpars == 0){
      mPhi_Next.zeros();
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
            mPhi_Next(k,t) += accu(cPMX.slice(m).col(t)%mY.col(k));
          }
          mPhi_Next(k,t) = probcheck(mPhi_Next(k,t)/foo);
          mPMsumX.col(t) = foovec;
        }
      }
    }
    if(fixedpars > 0){
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
          }
          mPMsumX.col(t) = foovec;
        }
      }
    }
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    
    cGamma    = cGamma_Next;
    mDelta    = mDelta_Next;
    cPi = cPi_foo;
    cLogPi    = log(cPi);
    mOmega    = mOmega_Next;
    mLogOmega = log(mOmega);
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);
  
  double BIClow;
  double BIChigh;
  BIClow  = -2.0*LLKSeries(iter-1) + log(iN)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  
  double AIC;
  AIC = -2.0*LLKSeries(iter-1) + 2.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0)*iPh);
  
  // computing log-linear parameters
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  
  
  
  arma::mat mBeta = zeros(nfreepar_res,iT);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(t = 0; t < iT; t++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(iItemfoo,t)/(1-mPhi(iItemfoo,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,t));
        }
      }else{
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(0,t)/(1-mPhi(0,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(1,ivItemcat(v)-1)/mPhi(0,t));
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
  
  arma::vec parvec = join_cols(join_cols(vectorise(mDelta),vectorise(cGamma)),vectorise(mBeta));
  
  arma::mat mBeta_Score=zeros(iN,nfreepar_res*iT);
  iRoll = 0;
  iItemfoo = 0;
  for(t = 0; t < iT; t++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) = 
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(iItemfoo) - mPhi(iItemfoo,t));
          iItemfoo += ivItemcat(v)-1;
        }
        
      }else{
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) = 
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(t).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(0) - mPhi(0,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);  
    }
  } 
  iItemfoo = 0;
  
  arma::mat mDelta_Score_out(iN,(iM-1)*iPh);
  foo= 0;
  for(j = 0; j < iJ; j++){
    mDelta_Score_out.rows(foo,foo + vNj(j)-1)= repmat(mDelta_Score.row(j),vNj(j),1); 
    foo = foo + vNj(j);
  }
  arma::mat mScore = join_rows(mDelta_Score_out,mGamma_Score);
  if(nsteps != 3){
    mScore = join_rows(mScore,mBeta_Score);
  }
  if(nsteps == 3){
    nfreepar_res = 0;
  }
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  // asymptotic SEs correction
  int uncondLatpars   = (iM-1) + (iT-1)*iM;
  int parsfree        = (iT - 1)*iP*iM + (iM - 1)*iPh;
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + nfreepar_res-1,uncondLatpars + nfreepar_res-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16); 
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + nfreepar_res-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixedpars==1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  
  if(nsteps==3){
    mV2.eye();
    mQ.eye();
    mVar_corr.eye();
    SEs_cor =  SEs_unc;
  }
  
  List EMout;
  EMout["mPhi"]        = mPhi;
  EMout["cGamma"]      = cGamma;
  EMout["cGamma_Info"] = cGamma_Info;
  EMout["mGamma_Score"] = mGamma_Score;
  EMout["mDelta"]      = mDelta;
  EMout["cDelta_Info"] = cDelta_Info;
  EMout["mDelta_Score"] = mDelta_Score;
  EMout["cPi"]       = cPi;
  EMout["mOmega"]    = mOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIChigh"]       = BIChigh;
  EMout["BIClow"]       = BIClow;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["AIC"]       = AIC;  
  EMout["mScore"]    = mScore;
  EMout["Infomat"]   = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs_unc"]   = SEs_unc;
  EMout["SEs_cor"]   = SEs_cor;
  EMout["mV2"]       = mV2;
  EMout["mQ"]        = mQ;
  EMout["mVar_corr"] = mVar_corr;
  
  return EMout;
  
}

// [[Rcpp::export]]
List MLTLCA_cov_poly(arma::mat mY, arma::mat mZ, arma::vec vNj, arma::vec vOmega_start, arma::cube cGamma_start, arma::mat mPhi_start, arma::mat mStep1Var, arma::ivec ivItemcat, int maxIter = 1e3, double tol = 1e-8, int fixedpars = 0, int nsteps = 1, double NRtol = 1e-6, int NRmaxit = 100){
  // mY is iJ*sum(nj) x K
  // mZ is iJ*sum(nj) x iP, where iP is the number of iP-1 covariates (can even be 1) +
  // a column vector of ones to include the intercept
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // cPi is iN x iT x iM
  // cGamma is iT-1 x iP x iM
  // first Class is taken as reference
  // ivItemcat is the vector of number of categories for each item
  // fixedpars = 0 for one-step estimator
  // fixedpars = 1 for "pure" two-step (mPhi fixed, estimated using LCAfit or MLTLCA)
  // fixedpars = 2 for two-stage (fixing vOmega and mPhi, estimated using MLTLCA)
  
  int iV    = ivItemcat.n_elem;
  int n,j,k,m,t,p,v;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int iT = cGamma_start.n_rows + 1 ;
  int iM = cGamma_start.n_slices;
  arma::cube cPi     = ones(iN,iT,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma_start.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPXag    = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPXag   = mPXag;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::mat mPMX = zeros(iN,iT*iM);
  arma::cube cLogdY = zeros(iN,iT,iV);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::vec vOmega = vOmega_start;
  arma::cube cGamma = cGamma_start;
  arma::mat mPhi = mPhi_start;
  arma::vec vLogOmega = log(vOmega);
  arma::cube cGamma_Next = cGamma;
  arma::mat mGamma = zeros((iT-1)*iM,iP);
  arma::cube cGamma_Info(iP,iP,iM*(iT-1));
  arma::mat mGamma_Score(iN,(iT-1)*iP*iM);
  arma::vec vOmega_Next = vOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  int ifooDcat = 0;
  arma::vec foovec = zeros(iN,1);
  // 
  if(fixedpars>0){
    // compute log-densities
    for(n = 0; n < iN; n++){
      // 
      for(t = 0; t < iT; t++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mY(n,ifooDcat) > 0.0){
                cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              }
              ifooDcat += 1;
            }
          }
          mLogdKY(n,t) +=cLogdY(n,t,v);
        }
      }
    }
  }
  
  List NR_step;
  while(eps > tol && iter<maxIter){
    if(fixedpars==0){
      // compute log-densities
      mLogdKY.zeros();
      for(n = 0; n < iN; n++){
        // 
        for(t = 0; t < iT; t++){
          ifooDcat = 0;
          for(v = 0; v< iV;v++){
            if(ivItemcat(v)==2){
              cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              ifooDcat += 1;
            }else{
              for(p = 0; p < ivItemcat(v); p++){
                if(mY(n,ifooDcat) > 0.0){
                  cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
                }
                ifooDcat += 1;
              }
            }
            mLogdKY(n,t) +=cLogdY(n,t,v);
          }
        }
      }
    }
    //
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPXag(j,m) = accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
        mlogPW(j,m) = vLogOmega(m) + mlogPXag(j,m);
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1);
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX  = exp(clogPMX);
    mPXag = exp(mlogPXag);
    // M step 
    //
    int iFoo1 = 0;
    int iFoo2 = 0;
    for(m = 0; m < iM; m++){
      NR_step = NR_step_covIT_wei(mZ, cGamma.slice(m), cPX.slice(m), mPW_N.col(m), NRtol, NRmaxit);
      arma::mat mGamma_foo = NR_step["beta"];
      arma::cube cGammaInfo_foo = NR_step["ibeta"];
      arma::mat mGammaScore_foo = NR_step["mSbeta"];
      arma::mat mPi_foo = NR_step["w_i"];
      cPi_foo.slice(m) = mPi_foo;
      cGamma_Next.slice(m) = mGamma_foo;
      cGamma_Info.slices(iFoo1,iFoo1 + iT-2) = cGammaInfo_foo;
      mGamma_Score.cols(iFoo2, iFoo2 + ((iT-1)*iP) -1) = mGammaScore_foo;
      iFoo1 = iFoo1 + iT-1;
      iFoo2 = iFoo2 + ((iT-1)*iP);
    }
    
    if(fixedpars != 2){
      vOmega_Next = sum(mPW,0).t();
      vOmega_Next = vOmega_Next/accu(vOmega_Next);
    }
    if(fixedpars == 0){
      mPhi_Next.zeros();
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
            mPhi_Next(k,t) += accu(cPMX.slice(m).col(t)%mY.col(k));
          }
          mPhi_Next(k,t) = probcheck(mPhi_Next(k,t)/foo);
          mPMsumX.col(t) = foovec;
        }
      }
    }
    if(fixedpars > 0){
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
          }
          mPMsumX.col(t) = foovec;
        }
      }
    }
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    
    cGamma    = cGamma_Next;
    cPi = cPi_foo;
    cLogPi    = log(cPi);
    vOmega    = vOmega_Next;
    vLogOmega = log(vOmega);
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);
  double BIClow;
  double BIChigh;
  double AIC;
  BIClow  = -2.0*LLKSeries(iter-1) + log(iN)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0));
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*1.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0));
  AIC     = -2.0*LLKSeries(iter-1) + 2.0*(iT*nfreepar_res + (iT - 1.0)*iP*iM + (iM - 1.0));
  
  // computing log-linear parameters
  arma::vec vDeltafoo(iM);
  vDeltafoo = log(vOmega/vOmega(0));
  arma::vec vDelta(iM-1);
  vDelta = vDeltafoo.subvec(1,iM-1);
  
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  
  
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  arma::mat mBeta = zeros(nfreepar_res,iT);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(t = 0; t < iT; t++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(iItemfoo,t)/(1-mPhi(iItemfoo,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,t));
        }
      }else{
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(0,t)/(1-mPhi(0,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(1,ivItemcat(v)-1)/mPhi(0,t));
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
  
  arma::vec parvec = join_cols(join_cols(vDelta,vectorise(cGamma)),vectorise(mBeta));
  
  arma::mat mOmega_Score(iN,iM-1);
  for(m =1; m< iM; m++){
    mOmega_Score.col(m-1) = mPW_N.col(m)*(1.0 - vOmega(m));
  }
  arma::mat mBeta_Score=zeros(iN,nfreepar_res*iT);
  iRoll = 0;
  iItemfoo = 0;
  for(t = 0; t < iT; t++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(iItemfoo) - mPhi(iItemfoo,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }else{
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(t).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(0) - mPhi(0,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);
    }
  }
  iItemfoo = 0;
  //
  
  arma::mat mScore = join_rows(mOmega_Score,mGamma_Score);
  if(nsteps != 3){
    mScore = join_rows(mScore,mBeta_Score);
  }
  
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  // asymptotic SEs correction
  int uncondLatpars   = (iM-1) + (iT-1)*iM;
  int parsfree        = (iT - 1)*iP*iM + (iM - 1);
  if(nsteps == 3){
    nfreepar_res = 0;
  }
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + nfreepar_res-1,uncondLatpars + nfreepar_res-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16);
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + nfreepar_res-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixedpars==1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  if(nsteps==3){
    mV2.eye();
    mQ.eye();
    mVar_corr.eye();
    SEs_cor =  SEs_unc;
  }
  
  List EMout;
  EMout["mPhi"]        = mPhi;
  EMout["cGamma"]      = cGamma;
  EMout["cGamma_Info"] = cGamma_Info;
  EMout["mGamma_Score"] = mGamma_Score;
  EMout["mOmega_Score"] = mOmega_Score;
  EMout["mBeta_Score"] = mBeta_Score;
  EMout["cPi"]       = cPi;
  EMout["vOmega"]    = vOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["vDelta"] = vDelta;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIChigh"]       = BIChigh;
  EMout["AIC"]       = AIC;
  EMout["BIClow"]       = BIClow;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  EMout["mScore"]    = mScore;
  EMout["Infomat"]   = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs_unc"]   = SEs_unc;
  EMout["SEs_cor"]   = SEs_cor;
  EMout["mV2"]       = mV2;
  EMout["mQ"]        = mQ;
  EMout["mVar_corr"] = mVar_corr;
  
  return EMout;
  
}

// [[Rcpp::export]]
List MLTLCA_poly(arma::mat mY, arma::vec vNj, arma::vec vOmega, arma::mat mPi, arma::mat mPhi, arma::ivec ivItemcat, arma::uvec first_poly, arma::vec reord_user, arma::vec reord_user_high, int maxIter = 1e3, double tol = 1e-8, int reord = 1){
  // mY is iJ*sum(nj) x K
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // mPi is iT x iM
  // ivItemcat is the vector of number of categories for each item
  
  int iV    = ivItemcat.n_elem;
  int n,j,k,m,t,p,v;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = 1;
  
  int iT = mPi.n_rows;
  int iM = mPi.n_cols;
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::cube cLogdY = zeros(iN,iT,iV);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::mat mLogPi = log(mPi);
  arma::vec vLogOmega = log(vOmega);
  arma::mat mPi_Next = mPi;
  arma::vec vOmega_Next = vOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  int ifooDcat = 0;
  arma::vec foovec = zeros(iN,1);
  
  while(eps > tol && iter<maxIter){
    // compute log-densities
    mLogdKY.zeros();
    for(n = 0; n < iN; n++){
      // 
      for(t = 0; t < iT; t++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mY(n,ifooDcat) > 0.0){
                cLogdY(n,t,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,t), isize);
              }
              ifooDcat += 1;
            }
          }
          mLogdKY(n,t) +=cLogdY(n,t,v);
        }
      }
    }
    // 
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = mLogPi(t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(mPi.col(m),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPW(j,m) = vLogOmega(m) + accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1); 
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX = exp(clogPMX);
    // 
    // M step
    //
    mPhi_Next.zeros();
    mPi_Next = sum(cPMX,0);
    vOmega_Next = sum(mPW,0).t();
    vOmega_Next = vOmega_Next/accu(vOmega_Next);
    for(m = 0; m < iM; m++){
      mPi_Next.col(m) = mPi_Next.col(m)/accu(mPi_Next.col(m)); 
    }
    
    mPhi_Next.zeros();
    for(k = 0; k < iK; k++){
      for(t = 0; t < iT; t++){
        foo = 0.0;
        foovec.zeros();
        for(m = 0; m < iM; m++){
          foo += accu(cPMX.slice(m).col(t));
          foovec = foovec + cPMX.slice(m).col(t);
          mPhi_Next(k,t) += accu(cPMX.slice(m).col(t)%mY.col(k));
        }
        mPhi_Next(k,t) = probcheck(mPhi_Next(k,t)/foo);
        mPMsumX.col(t) = foovec;
      }
    }
    for(m = 0; m < iM; m++){
      mPi_Next.col(m) = OmegaCheck(mPi_Next.col(m), iT);
    }
    vOmega_Next = OmegaCheck(vOmega_Next, iM);
    
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    mPi       = mPi_Next;
    vOmega    = vOmega_Next;
    mLogPi    = log(mPi);
    vLogOmega = log(vOmega);
    
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::cube cPhi_foo = zeros(iK,iT,iM);
    // int fp  = first_poly.n_elem;
    int fp  = iT;
    arma::mat mPhisum = zeros(fp,iM);
    for(k = 0; k < iK; k++){
      for(m = 0; m < iM; m++){ 
        for(t = 0; t < iT; t++){
          cPhi_foo(k,t,m) = accu(cPMX.slice(m).col(t)%mY.col(k))/accu(cPMX.slice(m).col(t));
        }
      }
    }
    
    
    arma::umat mLow_order(iT,iM);
    for(m = 0; m < iM; m++){
      mPhisum.col(m) = sum(cPhi_foo.slice(m).rows(first_poly)).t();
      mLow_order.col(m) = sort_index(mPhisum.col(m),"descending"); 
    }
    arma::vec vPhisum = sum(mPhi.rows(first_poly)).t();
    arma::uvec low_order = sort_index(vPhisum,"descending");
    arma::uvec high_order = sort_index(vOmega,"descending");
    int ireord_check = 0;
    for(t = 0; t < iT; t++){
      if(reord_user(t) != t){
        ireord_check = 1;
      }
    }
    if(ireord_check == 1){
      arma::umat mLow_order_foo = mLow_order;
      arma::uvec low_order_foo = low_order;
      int ifoo_reord_user;
      for(t = 0; t < iT; t++){
        ifoo_reord_user = reord_user(t);
        low_order_foo(t) = low_order(ifoo_reord_user);
        for(m = 0; m < iM; m++){
          mLow_order_foo(t,m) = mLow_order(ifoo_reord_user,m);
        }
      }
      low_order = low_order_foo;
      mLow_order = mLow_order_foo;
    }
    //
    int ireord_check_high = 0;
    for(m = 0; m < iM; m++){
      if(reord_user_high(m) != m){
        ireord_check_high = 1;
      }
    }
    if(ireord_check_high == 1){
      arma::umat mLow_order_foo = mLow_order;
      arma::uvec high_order_foo = high_order;
      int ifoo_reord_user_high;
      for(m = 0; m < iM; m++){
        ifoo_reord_user_high = reord_user_high(m);
        high_order_foo(m) = high_order(ifoo_reord_user_high);
        for(t = 0; t < iT; t++){
          mLow_order_foo(t,m) = mLow_order(t,ifoo_reord_user_high);
        }
      }
      high_order = high_order_foo;
      mLow_order = mLow_order_foo;
    }
    //
    int ifoo = 0;
    int ifoo_high = 0;
    arma::vec vOmega_sorted = vOmega;
    arma::mat mPhi_sorted = mPhi;
    arma::mat mPi_sorted = mPi;
    arma::mat mPW_sorted = mPW;
    arma::mat mPW_N_sorted = mPW_N;
    arma::cube cPX_sorted = cPX;
    arma::cube cPMX_sorted = cPMX;
    arma::cube clogPX_sorted = cPX;
    arma::cube clogPMX_sorted = cPMX;
    arma::mat mSumPX_sorted = mSumPX;
    arma::mat mPMsumX_sorted = mPMsumX;
    arma::cube cLogdY_sorted = cLogdY;
    arma::mat mLogdKY_sorted = mLogdKY;
    for(t=0; t< iT; t++){
      ifoo               = low_order(t);
      mPhi_sorted.col(t) = mPhi.col(ifoo);
      mPMsumX_sorted.col(t)  = mPMsumX.col(ifoo);
      mLogdKY_sorted.col(t)  = mLogdKY.col(ifoo);
      for(m = 0; m < iM; m++){
        ifoo =  mLow_order(t,m);
        mPi_sorted(t,m)    = mPi(ifoo,m);
        cPX_sorted.slice(m).col(t)  = cPX.slice(m).col(ifoo);
        clogPX_sorted.slice(m).col(t)  = clogPX.slice(m).col(ifoo);
        cPMX_sorted.slice(m).col(t)  = cPMX.slice(m).col(ifoo);
        clogPMX_sorted.slice(m).col(t)  = clogPMX.slice(m).col(ifoo);
      }
    }
    for(m = 0; m < iM; m++){
      ifoo_high = high_order(m);
      vOmega_sorted(m) = vOmega(ifoo_high);
      mPi_sorted.col(m)  = mPi.col(ifoo_high);
      cPX_sorted.slice(m)  = cPX.slice(ifoo_high);
      clogPX_sorted.slice(m)  = clogPX.slice(ifoo_high);
      cPMX_sorted.slice(m)  = cPMX.slice(ifoo_high);
      clogPMX_sorted.slice(m)  = clogPMX.slice(ifoo_high);
      mSumPX_sorted.col(m)  = mSumPX.col(ifoo_high);
    }
    vOmega = vOmega_sorted;
    mPhi = mPhi_sorted;
    mPi = mPi_sorted;
    cPX = cPX_sorted;
    clogPX = clogPX_sorted;
    cPMX = cPMX_sorted;
    clogPMX = clogPMX_sorted;
    mPMsumX = mPMsumX_sorted;
    mLogdKY = mLogdKY_sorted;
    mSumPX = mSumPX_sorted;
  }
  
  
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);
  
  double BIClow;
  double BIChigh;
  double AIC;
  BIClow = -2.0*LLKSeries(iter-1) + log(iN)*1.0*(iT*nfreepar_res + (iT - 1.0)*iM + (iM - 1.0));
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*1.0*(iT*nfreepar_res + (iT - 1.0)*iM + (iM - 1.0));
  AIC = -2.0*LLKSeries(iter-1) + 2.0*(iT*nfreepar_res + (iT - 1.0)*iM + (iM - 1.0));
  
  // computing log-linear parameters
  arma::vec vDeltafoo(iM);
  vDeltafoo = log(vOmega/vOmega(0));
  arma::vec vDelta = vDeltafoo.subvec(1,iM-1);
  arma::mat mGammafoo(iT,iM);
  arma::mat mGamma(iT-1,iM);
  for(m=0; m < iM; m++){
    mGammafoo.col(m) = log(mPi.col(m)/mPi(0,m));
    mGamma.col(m) = mGammafoo.col(m).subvec(1,iT-1);
  }
  //
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  
  
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  
  arma::mat mBeta = zeros(nfreepar_res,iT);
  int iRoll = 0;
  int iItemfoo = 0;
  for(v  = 0; v < iV; v++){
    if(ivItemcat(v)==2){
      ivItemcat(v) = 1;
    }
  }
  for(t = 0; t < iT; t++){
    iRoll = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        iItemfoo = sum(ivItemcat.subvec(0,v-1));
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(iItemfoo,t)/(1-mPhi(iItemfoo,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + ivItemcat(v)-1)/mPhi(iItemfoo,t));
        }
      }else{
        if(ivItemcat(v)==1){
          mBeta(iRoll,t) = log(mPhi(0,t)/(1-mPhi(0,t)));
        } else{
          mBeta.col(t).subvec(iRoll,iRoll + (ivItemcat(v)-2)) = log(mPhi.col(t).subvec(1,ivItemcat(v)-1)/mPhi(0,t));
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
  arma::vec parvec = join_cols(join_cols(vDelta,vectorise(mGamma)),vectorise(mBeta));
  
  // Computing the score
  
  arma::mat mOmega_Score(iN,iM-1);
  for(m =1; m< iM; m++){
    mOmega_Score.col(m-1) = mPW_N.col(m)*(1.0 - vOmega(m));
  }
  arma::mat mGamma_Score(iN,(iT-1)*iP*iM);
  int iFoo2 = 0;
  arma::mat mGammaScore_foo = zeros(iN,(iT-1)*iP);
  for(m = 0; m < iM; m++){
    mGammaScore_foo.fill(0.0);
    for(t = 1; t < iT; t++){
      mGammaScore_foo.col(t-1) = (cPX.slice(m).col(t) - mPi(t,m))%(mPW_N.col(m));
    }
    
    mGamma_Score.cols(iFoo2, iFoo2 + ((iT-1)*iP) -1) = mGammaScore_foo;
    iFoo2 = iFoo2 + ((iT-1)*iP);
  }
  iFoo2 = 0;
  
  arma::mat mBeta_Score=zeros(iN,nfreepar_res*iT);
  iRoll = 0;
  iItemfoo = 0;
  for(t = 0; t < iT; t++){
    iItemfoo = 0;
    for(v  = 0; v < iV; v++){
      if(v>0){
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2) - repmat(mPhi.col(t).subvec(iItemfoo + 1,iItemfoo + 1+ ivItemcat(v)-2).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(iItemfoo) - mPhi(iItemfoo,t));
          iItemfoo += ivItemcat(v)-1;
        }
        
      }else{
        if(ivItemcat(v)>2){
          mBeta_Score.cols(iRoll,iRoll + (ivItemcat(v)-2)) =
            repmat(mPMsumX.col(t),1,ivItemcat(v)-1)%(mY.cols(1,ivItemcat(v)-1) - repmat(mPhi.col(t).subvec(1,ivItemcat(v)-1).t(),iN,1));
          iItemfoo += ivItemcat(v);
        }else{
          mBeta_Score.col(iRoll) = mPMsumX.col(t)%(mY.col(0) - mPhi(0,t));
          iItemfoo += ivItemcat(v)-1;
        }
      }
      iRoll += (ivItemcat(v)-1);
    }
  }
  iItemfoo = 0;
  
  
  arma::mat mScore = join_rows(join_rows(mOmega_Score,mGamma_Score),mBeta_Score);
  arma::mat Infomat = mScore.t()*mScore/(1.0*iN);
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/(1.0*iN);
  // arma::mat Varmat = inv(Infomat)/(1.0*iN);
  arma::vec SEs =  sqrt(Varmat.diag());
  
  List EMout;
  EMout["mScore"]    = mScore;
  EMout["Infomat"]    = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs"]    = SEs;
  EMout["mPhi"]      = mPhi;
  EMout["mPi"]       = mPi;
  EMout["vOmega"]    = vOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["mPXmarg"]   = mPXmarg;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["vDelta"] = vDelta;
  EMout["mGamma"] = mGamma;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIClow"]       = BIClow;
  EMout["BIChigh"]       = BIChigh;
  EMout["AIC"]       = AIC;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  
  return EMout;
}

// [[Rcpp::export]]
List MLTLCA_covlowhigh(arma::mat mY, arma::mat mZ, arma::mat mZh, arma::vec vNj, arma::mat mDelta_start, arma::cube cGamma_start, arma::mat mPhi_start, arma::mat mStep1Var, int maxIter = 1e3, double tol = 1e-8, int fixedpars = 0, double NRtol = 1e-6, int NRmaxit = 100){
  // mY is iJ*sum(nj) x K
  // mZ is iJ*sum(nj) x iP, where iP is the number of iP-1 covariates (can even be 1) + an intercept term
  // mZh is iJ x iPh, where iPh is the number of iP-1 covariates + an intercept term
  // a column vector of ones to include the intercept
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // cPi is iN x iT x iM
  // cGamma is iT-1 x iP x iM
  // mDelta is (iM -1) x iPh
  // first Class is taken as reference
  // 
  // fixedpars = 0 for one-step estimator
  // fixedpars = 1 for "pure" two-step (mPhi fixed, estimated using LCAfit or MLTLCA)
  // fixedpars = 2 for two-stage (fixing vOmega and mPhi, estimated using MLTLCA)
  
  int n,j,k,m,t;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int iPh = mZh.n_cols;
  int iT = cGamma_start.n_rows + 1 ;
  int iM = cGamma_start.n_slices;
  arma::vec vOnes = ones(iJ,1);
  arma::cube cPi     = ones(iN,iT,iM);
  arma::mat mOmega   = ones(iJ,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma_start.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  //
  // for(k = 0; k < (iK - 1); k++){
  //   mPg.col(k + 1) = exp(mZ*mBeta.col(k));
  // }
  // arma::vec sumPg = sum(mPg,1);
  // for(n = 0; n < iN; n++){
  //   mPg.row(n) = mPg.row(n)/sumPg(n);
  // }
  //
  
  
  arma::mat mDelta_st_foo = mDelta_start.t();
  for(m = 1; m < iM; m++){
    mOmega.col(m) = exp(mZh*mDelta_st_foo.col(m-1));
  }
  for(j = 0; j < iJ; j++){
    mOmega.row(j) = mOmega.row(j)/accu(mOmega.row(j));
  }
  
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPXag    = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPXag   = mPXag;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::mat mPMX = zeros(iN,iT*iM);
  arma::cube cLogdY = zeros(iN,iT,iK);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::mat mDelta    = mDelta_start;
  arma::mat mLogOmega = log(mOmega);
  arma::mat mDelta_Next = mDelta;
  
  arma::cube cGamma = cGamma_start;
  arma::mat mPhi = mPhi_start;
  arma::cube cGamma_Next = cGamma;
  arma::mat mGamma = zeros((iT-1)*iM,iP);
  arma::cube cGamma_Info(iP,iP,iM*(iT-1));
  arma::cube cDelta_Info(iPh,iPh,iM-1);
  arma::mat mDelta_Score(iN,iPh*(iM-1));
  arma::mat mGamma_Score(iN,(iT-1)*iP*iM);
  arma::mat mOmega_Next = mOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  arma::vec foovec = zeros(iN,1);
  // 
  if(fixedpars>0){
    // compute log-densities
    for(n = 0; n < iN; n++){
      for(t = 0; t < iT; t++){
        for(k = 0; k < iK; k++){
          cLogdY(n,t,k) = Rf_dbinom(mY(n,k), size, mPhi(k,t), isize);
          mLogdKY(n,t) +=cLogdY(n,t,k);
        }
      }
    }
  }
  
  List NR_step;
  List NR_step_Delta;
  while(eps > tol && iter<maxIter){
    if(fixedpars == 0){
      // compute log-densities
      mLogdKY.zeros();
      for(n = 0; n < iN; n++){
        for(t = 0; t < iT; t++){
          for(k = 0; k < iK; k++){
            cLogdY(n,t,k) = Rf_dbinom(mY(n,k), size, mPhi(k,t), isize);
            mLogdKY(n,t) +=cLogdY(n,t,k);
          }
        }
      }
    }
    //
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPXag(j,m) = accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
        mlogPW(j,m) = mLogOmega(j,m) + mlogPXag(j,m);
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1);
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX  = exp(clogPMX);
    mPXag = exp(mlogPXag);
    // M step 
    //
    int iFoo1 = 0;
    int iFoo2 = 0;
    for(m = 0; m < iM; m++){
      NR_step = NR_step_covIT_wei(mZ, cGamma.slice(m), cPX.slice(m), mPW_N.col(m), NRtol, NRmaxit);
      arma::mat mGamma_foo = NR_step["beta"];
      arma::cube cGammaInfo_foo = NR_step["ibeta"];
      arma::mat mGammaScore_foo = NR_step["mSbeta"];
      arma::mat mPi_foo = NR_step["w_i"];
      cPi_foo.slice(m) = mPi_foo;
      cGamma_Next.slice(m) = mGamma_foo;
      cGamma_Info.slices(iFoo1,iFoo1 + iT-2) = cGammaInfo_foo;
      mGamma_Score.cols(iFoo2, iFoo2 + ((iT-1)*iP) -1) = mGammaScore_foo;
      iFoo1 = iFoo1 + iT-1;
      iFoo2 = iFoo2 + ((iT-1)*iP);
    }
    // 
    if(fixedpars != 2){
      NR_step_Delta = NR_step_covIT(mZh, mDelta, mPW, NRtol, NRmaxit);
      // NR_step_Delta = NR_step_covIT_wei(mZh, mDelta,mPW, vOnes, NRtol, NRmaxit);
      arma::mat mDelta_foo = NR_step_Delta["beta"];
      arma::cube cDeltaInfo_foo = NR_step_Delta["ibeta"];
      arma::mat mDeltaScore_foo = NR_step_Delta["mSbeta"];
      arma::mat mOmega_foo = NR_step_Delta["w_i"];
      mOmega_Next = mOmega_foo;
      mDelta_Next = mDelta_foo;
      cDelta_Info = cDeltaInfo_foo;
      mDelta_Score = mDeltaScore_foo;
    }
    if(fixedpars == 0){
      mPhi_Next.zeros();
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
            mPhi_Next(k,t) += accu(cPMX.slice(m).col(t)%mY.col(k));
          }
          mPhi_Next(k,t) = probcheck(mPhi_Next(k,t)/foo);
          mPMsumX.col(t) = foovec;
        }
      }
    }
    if(fixedpars > 0){
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
          }
          mPMsumX.col(t) = foovec;
        }
      }
    }
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    
    cGamma    = cGamma_Next;
    mDelta    = mDelta_Next;
    cPi = cPi_foo;
    cLogPi    = log(cPi);
    mOmega    = mOmega_Next;
    mLogOmega = log(mOmega);
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  double BIClow;
  double BIChigh;
  BIClow  = -2.0*LLKSeries(iter-1) + log(iN)*(iK*iT + (iT - 1.0)*iP*iM + (iM - 1.0));
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*(iK*iT + (iT - 1.0)*iP*iM + (iM - 1.0));
  
  double AIC;
  AIC = -2.0*LLKSeries(iter-1) + 2.0*(iK*iT + (iT - 1.0)*iP*iM + (iM - 1.0));
  
  // computing log-linear parameters
  
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  
  arma::mat mBeta = log(mPhi/(1.0 - mPhi));
  
  arma::vec parvec = join_cols(join_cols(vectorise(mDelta),vectorise(cGamma)),vectorise(mBeta));
  
  arma::mat mBeta_Score=zeros(iN,iK*iT);
  int iroll = 0;
  for(t = 0; t < iT; t++){
    for(k = 0; k < iK; k++){
      for(m = 0; m < iM; m++){
        mBeta_Score.col(iroll) += cPMX.slice(m).col(t)%(mY.col(k) - mPhi(k,t));
      }
      iroll += 1;
    }
  }
  
  arma::mat mDelta_Score_out(iN,(iM-1)*iPh);
  foo= 0;
  for(j = 0; j < iJ; j++){
    mDelta_Score_out.rows(foo,foo + vNj(j)-1)= repmat(mDelta_Score.row(j),vNj(j),1); 
    foo = foo + vNj(j);
  }
  
  arma::mat mScore = join_rows(join_rows(mDelta_Score_out,mGamma_Score),mBeta_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  // asymptotic SEs correction
  int uncondLatpars   = (iM-1) + (iT-1)*iM;
  int parsfree        = (iT - 1)*iP*iM + (iM - 1)*iPh;
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + (iT*iK)-1,uncondLatpars + (iT*iK)-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16); 
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + (iT*iK)-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixedpars==1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  List EMout;
  EMout["mPhi"]        = mPhi;
  EMout["cGamma"]      = cGamma;
  EMout["cGamma_Info"] = cGamma_Info;
  EMout["mGamma_Score"] = mGamma_Score;
  EMout["mDelta"]      = mDelta;
  EMout["cDelta_Info"] = cDelta_Info;
  EMout["mDelta_Score"] = mDelta_Score;
  EMout["cPi"]       = cPi;
  EMout["mOmega"]    = mOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIChigh"]       = BIChigh;
  EMout["BIClow"]       = BIClow;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["AIC"]       = AIC;  
  EMout["mScore"]    = mScore;
  EMout["Infomat"]   = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs_unc"]   = SEs_unc;
  EMout["SEs_cor"]   = SEs_cor;
  EMout["mV2"]       = mV2;
  EMout["mQ"]        = mQ;
  EMout["mVar_corr"] = mVar_corr;
  
  return EMout;
  
}

// [[Rcpp::export]]
List MLTLCA_cov(arma::mat mY, arma::mat mZ, arma::vec vNj, arma::vec vOmega_start, arma::cube cGamma_start, arma::mat mPhi_start, arma::mat mStep1Var, int maxIter = 1e3, double tol = 1e-8, int fixedpars = 0, double NRtol = 1e-6, int NRmaxit = 100){
  // mY is iJ*sum(nj) x K
  // mZ is iJ*sum(nj) x iP, where iP is the number of iP-1 covariates (can even be 1) +
  // a column vector of ones to include the intercept
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // cPi is iN x iT x iM
  // cGamma is iT-1 x iP x iM
  // first Class is taken as reference
  // 
  // fixedpars = 0 for one-step estimator
  // fixedpars = 1 for "pure" two-step (mPhi fixed, estimated using LCAfit or MLTLCA)
  // fixedpars = 2 for two-stage (fixing vOmega and mPhi, estimated using MLTLCA)
  
  int n,j,k,m,t;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int iT = cGamma_start.n_rows + 1 ;
  int iM = cGamma_start.n_slices;
  arma::cube cPi     = ones(iN,iT,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma_start.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPXag    = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPXag   = mPXag;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::mat mPMX = zeros(iN,iT*iM);
  arma::cube cLogdY = zeros(iN,iT,iK);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::vec vOmega = vOmega_start;
  arma::cube cGamma = cGamma_start;
  arma::mat mPhi = mPhi_start;
  arma::vec vLogOmega = log(vOmega);
  arma::cube cGamma_Next = cGamma;
  arma::mat mGamma = zeros((iT-1)*iM,iP);
  arma::cube cGamma_Info(iP,iP,iM*(iT-1));
  arma::mat mGamma_Score(iN,(iT-1)*iP*iM);
  arma::vec vOmega_Next = vOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  arma::vec foovec = zeros(iN,1);
  // 
  if(fixedpars>0){
    // compute log-densities
    for(n = 0; n < iN; n++){
      for(t = 0; t < iT; t++){
        for(k = 0; k < iK; k++){
          cLogdY(n,t,k) = Rf_dbinom(mY(n,k), size, mPhi(k,t), isize);
          mLogdKY(n,t) +=cLogdY(n,t,k);
        }
      }
    }
  }
  
  List NR_step;
  while(eps > tol && iter<maxIter){
    if(fixedpars==0){
      // compute log-densities
      mLogdKY.zeros();
      for(n = 0; n < iN; n++){
        for(t = 0; t < iT; t++){
          for(k = 0; k < iK; k++){
            cLogdY(n,t,k) = Rf_dbinom(mY(n,k), size, mPhi(k,t), isize);
            mLogdKY(n,t) +=cLogdY(n,t,k);
          }
        }
      }
    }
    //
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPXag(j,m) = accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
        mlogPW(j,m) = vLogOmega(m) + mlogPXag(j,m);
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1);
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX  = exp(clogPMX);
    mPXag = exp(mlogPXag);
    // M step 
    //
    int iFoo1 = 0;
    int iFoo2 = 0;
    for(m = 0; m < iM; m++){
      NR_step = NR_step_covIT_wei(mZ, cGamma.slice(m), cPX.slice(m), mPW_N.col(m), NRtol, NRmaxit);
      arma::mat mGamma_foo = NR_step["beta"];
      arma::cube cGammaInfo_foo = NR_step["ibeta"];
      arma::mat mGammaScore_foo = NR_step["mSbeta"];
      arma::mat mPi_foo = NR_step["w_i"];
      cPi_foo.slice(m) = mPi_foo;
      cGamma_Next.slice(m) = mGamma_foo;
      cGamma_Info.slices(iFoo1,iFoo1 + iT-2) = cGammaInfo_foo;
      mGamma_Score.cols(iFoo2, iFoo2 + ((iT-1)*iP) -1) = mGammaScore_foo;
      iFoo1 = iFoo1 + iT-1;
      iFoo2 = iFoo2 + ((iT-1)*iP);
    }
    
    if(fixedpars != 2){
      vOmega_Next = sum(mPW,0).t();
      vOmega_Next = vOmega_Next/accu(vOmega_Next);
    }
    if(fixedpars == 0){
      mPhi_Next.zeros();
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
            mPhi_Next(k,t) += accu(cPMX.slice(m).col(t)%mY.col(k));
          }
          mPhi_Next(k,t) = probcheck(mPhi_Next(k,t)/foo);
          mPMsumX.col(t) = foovec;
        }
      }
    }
    if(fixedpars > 0){
      for(k = 0; k < iK; k++){
        for(t = 0; t < iT; t++){
          foo = 0.0;
          foovec.zeros();
          for(m = 0; m < iM; m++){
            foo += accu(cPMX.slice(m).col(t));
            foovec = foovec + cPMX.slice(m).col(t);
          }
          mPMsumX.col(t) = foovec;
        }
      }
    }
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    cGamma    = cGamma_Next;
    cPi = cPi_foo;
    cLogPi    = log(cPi);
    vOmega    = vOmega_Next;
    vLogOmega = log(vOmega);
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  double BIClow;
  double BIChigh;
  double AIC;
  BIClow  = -2.0*LLKSeries(iter-1) + log(iN)*(iK*iT + (iT - 1.0)*iP*iM + (iM - 1.0));
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*(iK*iT + (iT - 1.0)*iP*iM + (iM - 1.0));
  AIC     = -2.0*LLKSeries(iter-1) + 2.0*(iK*iT + (iT - 1.0)*iP*iM + (iM - 1.0));
  
  // computing log-linear parameters
  arma::vec vDeltafoo(iM);
  vDeltafoo = log(vOmega/vOmega(0));
  arma::vec vDelta(iM-1);
  vDelta = vDeltafoo.subvec(1,iM-1);
  
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  
  arma::mat mBeta = log(mPhi/(1.0 - mPhi));
  
  arma::vec parvec = join_cols(join_cols(vDelta,vectorise(cGamma)),vectorise(mBeta));
  
  arma::mat mOmega_Score(iN,iM-1);
  for(m =1; m< iM; m++){
    mOmega_Score.col(m-1) = mPW_N.col(m)*(1.0 - vOmega(m));
  }
  arma::mat mBeta_Score=zeros(iN,iK*iT);
  int iroll = 0;
  for(t = 0; t < iT; t++){
    for(k = 0; k < iK; k++){
      for(m = 0; m < iM; m++){
        mBeta_Score.col(iroll) += cPMX.slice(m).col(t)%(mY.col(k) - mPhi(k,t));
      }
      iroll += 1;
    }
  }
  arma::mat mScore = join_rows(join_rows(mOmega_Score,mGamma_Score),mBeta_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  // asymptotic SEs correction
  int uncondLatpars   = (iM-1) + (iT-1)*iM;
  int parsfree        = (iT - 1)*iP*iM + (iM - 1);
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars,uncondLatpars,uncondLatpars + (iT*iK)-1,uncondLatpars + (iT*iK)-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = psinv(mJmat,1000,2.220446e-16); 
  arma::mat mH        = Infomat.submat(0,parsfree,parsfree-1,parsfree + (iT*iK)-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixedpars==1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  List EMout;
  EMout["mPhi"]        = mPhi;
  EMout["cGamma"]      = cGamma;
  EMout["cGamma_Info"] = cGamma_Info;
  EMout["mGamma_Score"] = mGamma_Score;
  EMout["cPi"]       = cPi;
  EMout["vOmega"]    = vOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["vDelta"] = vDelta;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIChigh"]       = BIChigh;
  EMout["AIC"]       = AIC;
  EMout["BIClow"]       = BIClow;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  EMout["mScore"]    = mScore;
  EMout["Infomat"]   = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs_unc"]   = SEs_unc;
  EMout["SEs_cor"]   = SEs_cor;
  EMout["mV2"]       = mV2;
  EMout["mQ"]        = mQ;
  EMout["mVar_corr"] = mVar_corr;
  
  return EMout;
  
}

// [[Rcpp::export]]
List MLTLCA(arma::mat mY, arma::vec vNj, arma::vec vOmega, arma::mat mPi, arma::mat mPhi, arma::vec reord_user, arma::vec reord_user_high, int maxIter = 1e3, double tol = 1e-8, int reord = 1){
  // mY is iJ*sum(nj) x K
  // iN = iJ*sum(nj)
  // iK is the number of items
  // iT is the number of individual classes
  // iM is the number of group classes
  // vGrId is the vector of group id's (should start from zero)
  // mPhi is iK x iT
  // vW is iM x 1
  // mPi is iT x iM
  
  int n,j,k,m,t;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = 1;
  
  int iT = mPi.n_rows;
  int iM = mPi.n_cols;
  arma::mat mPW      = zeros(iJ,iM);
  arma::mat mPW_N    = zeros(iN,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mPMsumX  = zeros(iN,iT);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::cube cLogdY = zeros(iN,iT,iK);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::mat mLogPi = log(mPi);
  arma::vec vLogOmega = log(vOmega);
  arma::mat mPi_Next = mPi;
  arma::vec vOmega_Next = vOmega;
  arma::mat mPhi_Next = mPhi;
  arma::vec vLLK = zeros(iJ,1);
  arma::vec LLKSeries(maxIter);
  double eps = 1.0;
  double iter = 0.0;
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  arma::vec foovec = zeros(iN,1);
  
  while(eps > tol && iter<maxIter){
    // compute log-densities
    mLogdKY.zeros();
    for(n = 0; n < iN; n++){
      for(t = 0; t < iT; t++){
        for(k = 0; k < iK; k++){
          cLogdY(n,t,k) = Rf_dbinom(mY(n,k), size, mPhi(k,t), isize);
          mLogdKY(n,t) +=cLogdY(n,t,k);
        }
      }
    }
    // 
    // E step
    // (working with log-probabilities to avoid numerical over/underflow)
    for(n = 0; n < iN; n++){
      for(m = 0; m < iM; m++){
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = mLogPi(t,m) + mLogdKY(n,t);
        }
        mSumPX(n,m) = MixtDensityScale(mPi.col(m),mLogdKY.row(n).t(), iT);
        for(t = 0; t < iT; t++){
          clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
          cPX(n,t,m) = exp(clogPX(n,t,m));
        }
      }
    }
    foo = 0.0;
    for(j = 0; j < iJ; j++){
      dLK =0.0;
      for(m = 0; m < iM; m++){
        mlogPW(j,m) = vLogOmega(m) + accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
      }
      foomax = max(mlogPW.row(j));
      mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
      
      for (m = 0; m < iM; m++) {
        dLK += exp(mlogPW_sc(j,m));
      }
      vLLK(j) = foomax + log(dLK);
      
      mlogPW.row(j) = mlogPW.row(j) - vLLK(j);
      mPW.row(j) = exp(mlogPW.row(j));
      for(m = 0; m <iM; m++){
        for(t = 0; t < iT; t++){
          clogPMX.slice(m).col(t).subvec(foo,foo + vNj[j]-1) = mlogPW(j,m) + clogPX.slice(m).col(t).subvec(foo,foo + vNj[j]-1); 
        }
        mPW_N.col(m).subvec(foo,foo + vNj[j]-1).fill(mPW(j,m)); 
      }
      foo = foo + vNj[j];
    }
    cPMX = exp(clogPMX);
    // 
    // M step
    //
    mPhi_Next.zeros();
    mPi_Next = sum(cPMX,0);
    vOmega_Next = sum(mPW,0).t();
    vOmega_Next = vOmega_Next/accu(vOmega_Next);
    for(m = 0; m < iM; m++){
      mPi_Next.col(m) = mPi_Next.col(m)/accu(mPi_Next.col(m)); 
    }
    
    mPhi_Next.zeros();
    for(k = 0; k < iK; k++){
      for(t = 0; t < iT; t++){
        foo = 0.0;
        foovec.zeros();
        for(m = 0; m < iM; m++){
          foo += accu(cPMX.slice(m).col(t));
          foovec = foovec + cPMX.slice(m).col(t);
          mPhi_Next(k,t) += accu(cPMX.slice(m).col(t)%mY.col(k));
        }
        mPhi_Next(k,t) = probcheck(mPhi_Next(k,t)/foo);
        mPMsumX.col(t) = foovec;
      }
    }
    for(m = 0; m < iM; m++){
      mPi_Next.col(m) = OmegaCheck(mPi_Next.col(m), iT);
    }
    vOmega_Next = OmegaCheck(vOmega_Next, iM);
    
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
    
    // Update parameters
    
    mPhi      = mPhi_Next;
    mPi       = mPi_Next;
    vOmega    = vOmega_Next;
    mLogPi    = log(mPi);
    vLogOmega = log(vOmega);
    
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::cube cPhi_foo = zeros(iK,iT,iM);
    // int fp  = iK;
    int fp  = iT;
    arma::mat mPhisum = zeros(fp,iM);
    for(k = 0; k < iK; k++){
      for(m = 0; m < iM; m++){ 
        for(t = 0; t < iT; t++){
          cPhi_foo(k,t,m) = accu(cPMX.slice(m).col(t)%mY.col(k))/accu(cPMX.slice(m).col(t));
        }
      }
    }
    
    
    arma::umat mLow_order(iT,iM);
    for(m = 0; m < iM; m++){
      mPhisum.col(m) = sum(cPhi_foo.slice(m)).t();
      mLow_order.col(m) = sort_index(mPhisum.col(m),"descending"); 
    }
    arma::vec vPhisum = sum(mPhi).t();
    arma::uvec low_order = sort_index(vPhisum,"descending");
    arma::uvec high_order = sort_index(vOmega,"descending");
    int ireord_check = 0;
    for(t = 0; t < iT; t++){
      if(reord_user(t) != t){
        ireord_check = 1;
      }
    }
    if(ireord_check == 1){
      arma::umat mLow_order_foo = mLow_order;
      arma::uvec low_order_foo = low_order;
      int ifoo_reord_user;
      for(t = 0; t < iT; t++){
        ifoo_reord_user = reord_user(t);
        low_order_foo(t) = low_order(ifoo_reord_user);
        for(m = 0; m < iM; m++){
          mLow_order_foo(t,m) = mLow_order(ifoo_reord_user,m);
        }
      }
      low_order = low_order_foo;
      mLow_order = mLow_order_foo;
    }
    //
    int ireord_check_high = 0;
    for(m = 0; m < iM; m++){
      if(reord_user_high(m) != m){
        ireord_check_high = 1;
      }
    }
    if(ireord_check_high == 1){
      arma::umat mLow_order_foo = mLow_order;
      arma::uvec high_order_foo = high_order;
      int ifoo_reord_user_high;
      for(m = 0; m < iM; m++){
        ifoo_reord_user_high = reord_user_high(m);
        high_order_foo(m) = high_order(ifoo_reord_user_high);
        for(t = 0; t < iT; t++){
          mLow_order_foo(t,m) = mLow_order(t,ifoo_reord_user_high);
        }
      }
      high_order = high_order_foo;
      mLow_order = mLow_order_foo;
    }
    //
    int ifoo = 0;
    int ifoo_high = 0;
    arma::vec vOmega_sorted = vOmega;
    arma::mat mPhi_sorted = mPhi;
    arma::mat mPi_sorted = mPi;
    arma::mat mPW_sorted = mPW;
    arma::mat mPW_N_sorted = mPW_N;
    arma::cube cPX_sorted = cPX;
    arma::cube cPMX_sorted = cPMX;
    arma::cube clogPX_sorted = cPX;
    arma::cube clogPMX_sorted = cPMX;
    arma::mat mSumPX_sorted = mSumPX;
    arma::mat mPMsumX_sorted = mPMsumX;
    arma::cube cLogdY_sorted = cLogdY;
    arma::mat mLogdKY_sorted = mLogdKY;
    for(t=0; t< iT; t++){
      ifoo               = low_order(t);
      mPhi_sorted.col(t) = mPhi.col(ifoo);
      mPMsumX_sorted.col(t)  = mPMsumX.col(ifoo);
      mLogdKY_sorted.col(t)  = mLogdKY.col(ifoo);
      for(m = 0; m < iM; m++){
        ifoo =  mLow_order(t,m);
        mPi_sorted(t,m)    = mPi(ifoo,m);
        cPX_sorted.slice(m).col(t)  = cPX.slice(m).col(ifoo);
        clogPX_sorted.slice(m).col(t)  = clogPX.slice(m).col(ifoo);
        cPMX_sorted.slice(m).col(t)  = cPMX.slice(m).col(ifoo);
        clogPMX_sorted.slice(m).col(t)  = clogPMX.slice(m).col(ifoo);
      }
    }
    for(m = 0; m < iM; m++){
      ifoo_high = high_order(m);
      vOmega_sorted(m) = vOmega(ifoo_high);
      mPi_sorted.col(m)  = mPi.col(ifoo_high);
      cPX_sorted.slice(m)  = cPX.slice(ifoo_high);
      clogPX_sorted.slice(m)  = clogPX.slice(ifoo_high);
      cPMX_sorted.slice(m)  = cPMX.slice(ifoo_high);
      clogPMX_sorted.slice(m)  = clogPMX.slice(ifoo_high);
      mSumPX_sorted.col(m)  = mSumPX.col(ifoo_high);
    }
    vOmega = vOmega_sorted;
    mPhi = mPhi_sorted;
    mPi = mPi_sorted;
    cPX = cPX_sorted;
    clogPX = clogPX_sorted;
    cPMX = cPMX_sorted;
    clogPMX = clogPMX_sorted;
    mPMsumX = mPMsumX_sorted;
    mLogdKY = mLogdKY_sorted;
    mSumPX = mSumPX_sorted;
  }
  
  
  
  double BIClow;
  double BIChigh;
  double AIC;
  BIClow = -2.0*LLKSeries(iter-1) + log(iN)*(iK*iT + (iT - 1.0)*iM + (iM - 1.0));
  BIChigh = -2.0*LLKSeries(iter-1) + log(iJ)*(iK*iT + (iT - 1.0)*iM + (iM - 1.0));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iK*iT + (iT - 1.0)*iM + (iM - 1.0));
  
  // computing log-linear parameters
  arma::vec vDeltafoo(iM);
  vDeltafoo = log(vOmega/vOmega(0));
  arma::vec vDelta = vDeltafoo.subvec(1,iM-1);
  arma::mat mGammafoo(iT,iM);
  arma::mat mGamma(iT-1,iM);
  for(m=0; m < iM; m++){
    mGammafoo.col(m) = log(mPi.col(m)/mPi(0,m));
    mGamma.col(m) = mGammafoo.col(m).subvec(1,iT-1);
  }
  // 
  arma::vec vPosthigh(iM);
  for(j = 0; j < iJ; j++){
    vPosthigh = PostCheck(mPW.row(j).t(),iM);
    mPW.row(j) = vPosthigh.t();
  }
  mlogPW = log(mPW);
  arma::vec vPW = mean(mPW).t();
  double Terr_high = iJ*accu(-vPW%log(vPW));
  mlogPW.elem( find_nonfinite(mlogPW) ).zeros();
  mlogPW.elem( find_nan(mlogPW) ).zeros();
  double Perr_high = accu(-mPW%mlogPW);
  double R2entr_high = 1.0-(Perr_high/Terr_high);
  arma::mat mPXmarg = sum(cPMX,2);
  arma::vec vFooPos(iT);
  for(n = 0; n< iN; n++){
    vFooPos = PostCheck( mPXmarg.row(n).t(),iT);
    mPXmarg.row(n) = vFooPos.t();
  }
  arma::vec vPimarg = mean(mPXmarg).t();
  double Terr_low = iN*accu(-vPimarg%log(vPimarg));
  arma::mat mlogPXmarg = log(mPXmarg);
  mlogPXmarg.elem(find_nonfinite(mlogPXmarg) ).zeros();
  mlogPXmarg.elem(find_nan(mlogPXmarg) ).zeros();
  double Perr_low = accu(-mPXmarg%mlogPXmarg);
  double R2entr_low = (Terr_low - Perr_low)/Terr_low;
  
  
  double ICL_BIClow;
  double ICL_BIChigh;
  ICL_BIClow = BIClow + 2.0*Perr_low;
  ICL_BIChigh = BIChigh + 2.0*Perr_high;
  
  arma::mat mBeta = log(mPhi/(1.0 - mPhi));
  
  arma::vec parvec = join_cols(join_cols(vDelta,vectorise(mGamma)),vectorise(mBeta));
  
  // Computing the score
  
  arma::mat mOmega_Score(iN,iM-1);
  for(m =1; m< iM; m++){
    mOmega_Score.col(m-1) = mPW_N.col(m)*(1.0 - vOmega(m));
  }
  arma::mat mGamma_Score(iN,(iT-1)*iP*iM);
  int iFoo2 = 0;
  arma::mat mGammaScore_foo = zeros(iN,(iT-1)*iP);
  for(m = 0; m < iM; m++){
    mGammaScore_foo.fill(0.0);
    for(t = 1; t < iT; t++){
      mGammaScore_foo.col(t-1) = (cPX.slice(m).col(t) - mPi(t,m))%(mPW_N.col(m));
    }
    
    mGamma_Score.cols(iFoo2, iFoo2 + ((iT-1)*iP) -1) = mGammaScore_foo;
    iFoo2 = iFoo2 + ((iT-1)*iP);
  }
  iFoo2 = 0;
  
  arma::mat mBeta_Score=zeros(iN,iK*iT);
  int iroll = 0;
  for(t = 0; t < iT; t++){
    for(k = 0; k < iK; k++){
      for(m = 0; m < iM; m++){
        mBeta_Score.col(iroll) += cPMX.slice(m).col(t)%(mY.col(k) - mPhi(k,t));
      }
      iroll += 1;
    }
  }
  
  arma::mat mScore = join_rows(join_rows(mOmega_Score,mGamma_Score),mBeta_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = psinv(Infomat,1000,2.220446e-16)/iN;
  arma::vec SEs =  sqrt(Varmat.diag());
  
  List EMout;
  EMout["mScore"]    = mScore;
  EMout["Infomat"]    = Infomat;
  EMout["Varmat"]    = Varmat;
  EMout["SEs"]    = SEs;
  EMout["mPhi"]      = mPhi;
  EMout["mPi"]       = mPi;
  EMout["vOmega"]    = vOmega;
  EMout["cPMX"]      = cPMX;
  EMout["cLogPMX"]   = clogPMX;
  EMout["mPXmarg"]   = mPXmarg;
  EMout["cPX"]       = cPX;
  EMout["cLogPX"]    = clogPX;
  EMout["mSumPX"]    = mSumPX;
  EMout["mPW"]       = mPW;
  EMout["mPW_N"]     = mPW_N;
  EMout["mlogPW"]    = mlogPW;
  EMout["mPMsumX"]   = mPMsumX;
  EMout["vDelta"] = vDelta;
  EMout["mGamma"] = mGamma;
  EMout["mBeta"]  = mBeta;
  EMout["parvec"] = parvec;
  EMout["LLKSeries"] = LLKSeries;
  EMout["vLLK"]      = vLLK;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIClow"]       = BIClow;
  EMout["BIChigh"]       = BIChigh;
  EMout["AIC"]       = AIC;
  EMout["ICL_BIClow"]    = ICL_BIClow;
  EMout["ICL_BIChigh"]    = ICL_BIChigh;
  EMout["R2entr_high"]       = R2entr_high;
  EMout["R2entr_low"]       = R2entr_low;
  
  return EMout;
}

// [[Rcpp::export]]
double MLTLCA_cov_LLK(arma::vec parvec,arma::mat mY, arma::mat mZ, arma::vec vNj, int iM, int iT){
  int n,j,k,m,t,p;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int foopar = 0;
  arma::cube cGamma(iT-1,iP,iM);
  arma::mat mPhi(iK,iT);
  arma::vec vOmega = ones(iM,1);
  vOmega.subvec(1,iM-1) = exp(parvec.subvec(0,iM-2));
  vOmega = vOmega/accu(vOmega);
  foopar = foopar + iM-1;
  for(m = 0; m < iM; m++){
    for(p = 0; p < iP; p++){
      for(t =1; t< iT; t++){
        cGamma(t-1,p,m) = parvec(foopar);
        foopar = foopar +1;
      }
    }
  }
  mPhi = reshape(conv_to<mat>::from(parvec(span(foopar,foopar + iK*iT -1))),iK,iT);
  mPhi = exp(mPhi)/(1.0 + exp(mPhi));
  arma::cube cPi     = ones(iN,iT,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::cube cLogdY = zeros(iN,iT,iK);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::vec vLogOmega = log(vOmega);
  arma::vec vLLK = zeros(iJ,1);
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  
  
  mLogdKY.zeros();
  for(n = 0; n < iN; n++){
    for(t = 0; t < iT; t++){
      for(k = 0; k < iK; k++){
        cLogdY(n,t,k) = Rf_dbinom(mY(n,k), size, mPhi(k,t), isize);
        mLogdKY(n,t) +=cLogdY(n,t,k);
      }
    }
  }
  for(n = 0; n < iN; n++){
    for(m = 0; m < iM; m++){
      for(t = 0; t < iT; t++){
        clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
      }
      mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
      for(t = 0; t < iT; t++){
        clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
        cPX(n,t,m) = exp(clogPX(n,t,m));
      }
    }
  }
  foo = 0.0;
  for(j = 0; j < iJ; j++){
    dLK =0.0;
    for(m = 0; m < iM; m++){
      mlogPW(j,m) = vLogOmega(m) + accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
    }
    foomax = max(mlogPW.row(j));
    mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
    
    for (m = 0; m < iM; m++) {
      dLK += exp(mlogPW_sc(j,m));
    }
    
    vLLK(j) = foomax + log(dLK);
    foo = foo + vNj[j];
  }
  
  double log_like = accu(vLLK);
  
  return(log_like);
  
}

// [[Rcpp::export]]
arma::vec MLTLCA_cov_LLK_j(arma::vec parvec,arma::mat mY, arma::mat mZ, arma::vec vNj, int iM, int iT){
  int n,j,k,m,t,p;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int iP = mZ.n_cols;
  int foopar = 0;
  arma::cube cGamma(iT-1,iP,iM);
  arma::mat mPhi(iK,iT);
  arma::vec vOmega = ones(iM,1);
  vOmega.subvec(1,iM-1) = exp(parvec.subvec(0,iM-2));
  vOmega = vOmega/accu(vOmega);
  foopar = foopar + iM-1;
  for(m = 0; m < iM; m++){
    for(p = 0; p < iP; p++){
      for(t =1; t< iT; t++){
        cGamma(t-1,p,m) = parvec(foopar);
        foopar = foopar +1;
      }
    }
  }
  mPhi = reshape(conv_to<mat>::from(parvec(span(foopar,foopar + iK*iT -1))),iK,iT);
  mPhi = exp(mPhi)/(1.0 + exp(mPhi));
  arma::cube cPi     = ones(iN,iT,iM);
  for(m = 0; m < iM; m++){
    for(n = 0; n< iN; n++){
      for(t = 1; t < iT; t++){
        cPi(n,t,m) = exp(accu(mZ.row(n)%cGamma.slice(m).row(t-1)));
      }
      cPi.slice(m).row(n) = cPi.slice(m).row(n)/accu(cPi.slice(m).row(n));
    }
  }
  arma::cube cPi_foo = cPi;
  arma::cube cLogPi  = log(cPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::cube cLogdY = zeros(iN,iT,iK);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::vec vLogOmega = log(vOmega);
  arma::vec vLLK = zeros(iJ,1);
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  
  
  mLogdKY.zeros();
  for(n = 0; n < iN; n++){
    for(t = 0; t < iT; t++){
      for(k = 0; k < iK; k++){
        cLogdY(n,t,k) = Rf_dbinom(mY(n,k), size, mPhi(k,t), isize);
        mLogdKY(n,t) +=cLogdY(n,t,k);
      }
    }
  }
  for(n = 0; n < iN; n++){
    for(m = 0; m < iM; m++){
      for(t = 0; t < iT; t++){
        clogPX(n,t,m) = cLogPi(n,t,m) + mLogdKY(n,t);
      }
      mSumPX(n,m) = MixtDensityScale(cPi.slice(m).row(n).t(),mLogdKY.row(n).t(), iT);
      for(t = 0; t < iT; t++){
        clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
        cPX(n,t,m) = exp(clogPX(n,t,m));
      }
    }
  }
  foo = 0.0;
  for(j = 0; j < iJ; j++){
    dLK =0.0;
    for(m = 0; m < iM; m++){
      mlogPW(j,m) = vLogOmega(m) + accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
    }
    foomax = max(mlogPW.row(j));
    mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
    
    for (m = 0; m < iM; m++) {
      dLK += exp(mlogPW_sc(j,m));
    }
    
    vLLK(j) = foomax + log(dLK);
    foo = foo + vNj[j];
  }
  
  return(vLLK);
  
}

// [[Rcpp::export]]
double MLTLCA_LLK(arma::vec parvec,arma::mat mY, arma::vec vNj, int iM, int iT){
  int n,j,k,m,t;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int foopar = 0;
  arma::mat mGamma(iT-1,iM);
  arma::mat mPhi(iK,iT);
  arma::vec vOmega = ones(iM,1);
  vOmega.subvec(1,iM-1) = exp(parvec.subvec(0,iM-2));
  vOmega = vOmega/accu(vOmega);
  foopar = foopar + iM-1;
  for(m = 0; m < iM; m++){
    for(t =1; t< iT; t++){
      mGamma(t-1,m) = parvec(foopar);
      foopar = foopar +1;
    }
  }
  mPhi = reshape(conv_to<mat>::from(parvec(span(foopar,foopar + iK*iT -1))),iK,iT);
  mPhi = exp(mPhi)/(1.0 + exp(mPhi));
  arma::mat mPi     = ones(iT,iM);
  for(m = 0; m < iM; m++){
    for(t = 1; t < iT; t++){
      mPi(t,m) = exp(mGamma(t-1,m));
    }
    mPi.col(m) = mPi.col(m)/accu(mPi.col(m));
  }
  
  arma::mat mPi_foo = mPi;
  arma::mat mLogPi  = log(mPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::cube cLogdY = zeros(iN,iT,iK);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::vec vLogOmega = log(vOmega);
  arma::vec vLLK = zeros(iJ,1);
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  
  
  mLogdKY.zeros();
  for(n = 0; n < iN; n++){
    for(t = 0; t < iT; t++){
      for(k = 0; k < iK; k++){
        cLogdY(n,t,k) = Rf_dbinom(mY(n,k), size, mPhi(k,t), isize);
        mLogdKY(n,t) +=cLogdY(n,t,k);
      }
    }
  }
  for(n = 0; n < iN; n++){
    for(m = 0; m < iM; m++){
      for(t = 0; t < iT; t++){
        clogPX(n,t,m) = mLogPi(t,m) + mLogdKY(n,t);
      }
      mSumPX(n,m) = MixtDensityScale(mPi.col(m),mLogdKY.row(n).t(), iT);
      for(t = 0; t < iT; t++){
        clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
        cPX(n,t,m) = exp(clogPX(n,t,m));
      }
    }
  }
  foo = 0.0;
  for(j = 0; j < iJ; j++){
    dLK =0.0;
    for(m = 0; m < iM; m++){
      mlogPW(j,m) = vLogOmega(m) + accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
    }
    foomax = max(mlogPW.row(j));
    mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
    
    for (m = 0; m < iM; m++) {
      dLK += exp(mlogPW_sc(j,m));
    }
    vLLK(j) = foomax + log(dLK);
    
    foo = foo + vNj[j];
  }
  
  double log_like = accu(vLLK);
  
  return(log_like);
  
}

// [[Rcpp::export]]
arma::vec MLTLCA_LLK_j(arma::vec parvec,arma::mat mY, arma::vec vNj, int iM, int iT){
  int n,j,k,m,t;
  int isize = 1;
  double size = 1.0;
  int iN = mY.n_rows;
  int iK = mY.n_cols;
  int iJ = vNj.n_elem;
  int foopar = 0;
  arma::mat mGamma(iT-1,iM);
  arma::mat mPhi(iK,iT);
  arma::vec vOmega = ones(iM,1);
  vOmega.subvec(1,iM-1) = exp(parvec.subvec(0,iM-2));
  vOmega = vOmega/accu(vOmega);
  foopar = foopar + iM-1;
  for(m = 0; m < iM; m++){
    for(t =1; t< iT; t++){
      mGamma(t-1,m) = parvec(foopar);
      foopar = foopar +1;
    }
  }
  mPhi = reshape(conv_to<mat>::from(parvec(span(foopar,foopar + iK*iT -1))),iK,iT);
  mPhi = exp(mPhi)/(1.0 + exp(mPhi));
  arma::mat mPi     = ones(iT,iM);
  for(m = 0; m < iM; m++){
    for(t = 1; t < iT; t++){
      mPi(t,m) = exp(mGamma(t-1,m));
    }
    mPi.col(m) = mPi.col(m)/accu(mPi.col(m));
  }
  
  arma::mat mPi_foo = mPi;
  arma::mat mLogPi  = log(mPi);
  arma::mat mPW      = zeros(iJ,iM);
  arma::cube cPX     = zeros(iN,iT,iM);
  arma::cube cPMX    = cPX;
  arma::mat mSumPX = zeros(iN,iM);
  arma::mat mlogPW   = mPW;
  arma::mat mlogPW_sc = mlogPW;
  arma::cube clogPX  = cPX;
  arma::cube clogPMX = cPMX;
  arma::cube cLogdY = zeros(iN,iT,iK);
  arma::mat mLogdKY = zeros(iN,iT);
  arma::vec vLogOmega = log(vOmega);
  arma::vec vLLK = zeros(iJ,1);
  double foo = 0.0;
  double foomax = 0.0;
  double dLK =0.0;
  
  
  mLogdKY.zeros();
  for(n = 0; n < iN; n++){
    for(t = 0; t < iT; t++){
      for(k = 0; k < iK; k++){
        cLogdY(n,t,k) = Rf_dbinom(mY(n,k), size, mPhi(k,t), isize);
        mLogdKY(n,t) +=cLogdY(n,t,k);
      }
    }
  }
  for(n = 0; n < iN; n++){
    for(m = 0; m < iM; m++){
      for(t = 0; t < iT; t++){
        clogPX(n,t,m) = mLogPi(t,m) + mLogdKY(n,t);
      }
      mSumPX(n,m) = MixtDensityScale(mPi.col(m),mLogdKY.row(n).t(), iT);
      for(t = 0; t < iT; t++){
        clogPX(n,t,m) = clogPX(n,t,m) - mSumPX(n,m);
        cPX(n,t,m) = exp(clogPX(n,t,m));
      }
    }
  }
  foo = 0.0;
  for(j = 0; j < iJ; j++){
    dLK =0.0;
    for(m = 0; m < iM; m++){
      mlogPW(j,m) = vLogOmega(m) + accu(mSumPX.col(m).subvec(foo, foo + vNj[j] -1));
    }
    foomax = max(mlogPW.row(j));
    mlogPW_sc.row(j) = mlogPW.row(j) - foomax;
    
    for (m = 0; m < iM; m++) {
      dLK += exp(mlogPW_sc(j,m));
    }
    vLLK(j) = foomax + log(dLK);
    
    foo = foo + vNj[j];
  }
  
  return(vLLK);
  
}
