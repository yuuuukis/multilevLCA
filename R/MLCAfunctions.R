.onAttach <- function(libname, pkgname){
  packageStartupMessage("multilevLCA: Estimates and plots single- and multilevel latent class models.")
}

multiLCA = function(data, Y, iT,
                    id_high = NULL, iM = NULL,
                    Z = NULL, Zh = NULL,
                    incomplete = FALSE, fixedslopes = FALSE,
                    startval = NULL, kmea = TRUE,
                    extout = FALSE, dataout = TRUE,
                    sequential = TRUE, numFreeCores = 2,
                    maxIter = 1e3, tol = 1e-8,
                    reord = TRUE,
                    reord_user = NULL, reord_user_high = NULL,
                    fixedpars = 1,
                    NRmaxit = 100, NRtol = 1e-6,
                    verbose = TRUE){
  check_inputs1(data,Y,iT,id_high,iM,Z,Zh,startval)
  data = data[,c(Y,id_high,Z,Zh,startval)]
  nrow_data = nrow(data)
  if(incomplete){
    if(nrow(na.omit(data[,Y]))==nrow_data){
      incomplete = FALSE
      mDesign = NULL
    } else{
      all_Y_missing = apply(data[,Y],1,function(x){sum(as.integer(!is.na(x)))==0})
      if(any(all_Y_missing)){
        if(verbose){
          pb=txtProgressBar(char="Discarding rows with missing values for all indicators...",width=1)
          setTxtProgressBar(pb,1)
          close(pb)
        }
      }
      data = data[!all_Y_missing,]
      mDesign = matrix(1,nrow(data[,Y]),ncol(data[,Y]))
      mDesign[is.na(data[,Y])] = 0
      colnames(mDesign) = Y
    }
  } else{
    data = data[complete.cases(data[,Y]),]
    if(nrow(data)<nrow_data){
      if(verbose){
        pb=txtProgressBar(char="Discarding rows with missing values for any indicator...",width=1)
        setTxtProgressBar(pb,1)
        close(pb)
      }
    }
    mDesign = NULL
  }
  approach  = check_inputs2(data,Y,iT,id_high,iM,Z,Zh,startval)
  
  if(is_tibble_manual(data)){
    data = as.data.frame(data)
  }
  
  if(length(iT)==1&(!is.null(reord_user))){
    if(any(reord_user!=c(1:iT))){
      reord = TRUE
    }
    reord_user = reord_user-1
  } else if(length(iT)==1&is.null(reord_user)){
    reord_user = 1:iT
    reord_user = reord_user-1
  }
  if(!is.null(iM)){
    if(length(iM)==1&(!is.null(reord_user_high))){
      if(any(reord_user_high!=c(1:iM))){
        reord = TRUE
      }
      reord_user_high = reord_user_high-1
    } else if(length(iM)==1&is.null(reord_user_high)){
      reord_user_high = 1:iM
      reord_user_high = reord_user_high-1
    }
  }
  reord = as.numeric(reord)
  fixed = fixedpars
  if((fixed==2)&is.null(iM)) fixed=fixedpars=1
  
  if(!is.null(id_high)){
    id_high_levs    = sort(unique(data[,id_high]))
    data[,id_high]  = as.numeric(factor(data[,id_high]))
    data            = data[order(data[,id_high]),]
    id_high_name    = id_high
    id_high         = data[,id_high]
    vNj             = table(id_high)
  }
  mY        = as.matrix(data[,Y])
  ivItemcat = apply(mY,2,function(x){length(na.omit(unique(x)))})
  itemchar  = any(apply(data[,Y],2,function(x){!is.numeric(x)}))
  if(any(ivItemcat>2)&!itemchar){
    mY  = update_YmY(mY,Y,ivItemcat,mDesign)
    Y   = mY$Y
    mY  = mY$mY
  } else if(itemchar){
    mY  = update_YmY_nonnum(mY,Y,ivItemcat,mDesign)
    Y   = mY$Y
    mY  = mY$mY
  } else if(incomplete){
    mY  = update_YmY(mY,Y,ivItemcat,mDesign)
    Y   = mY$Y
    mY  = mY$mY
  }
  if(any(ivItemcat>2)&incomplete){
    mDesign = matrix(1,nrow(mY),ncol(mY))
    mDesign[is.na(mY)] = 0
  }
  if(!is.null(Z)){
    mZ  = clean_cov(data,Z)
    Z   = mZ$covnam
    mZ  = mZ$cov
    }
  if(!is.null(Zh)){
    mZh = clean_cov(data,Zh)
    Zh  = mZh$covnam
    mZh = mZh$cov
  }
  if(!is.null(startval)){
    startval = data[,startval]
  }
  
  if(approach=="direct"){
    #
    if(verbose){
      if(fixed==0){
        pb=txtProgressBar(char="Fitting LC model...",width=1)
      } else{
        pb=txtProgressBar(char="Fitting measurement model...",width=1)
      }
      setTxtProgressBar(pb,1)
      close(pb)
    }
    #
    if(any(ivItemcat>2)|incomplete|(fixedslopes&!is.null(id_high)&!is.null(Z))){
      #
      if(is.null(id_high)&is.null(Z)&is.null(Zh)){
        out = LCA_fast_init_poly(mY,iT,ivItemcat,incomplete,kmea,maxIter,tol,reord,reord_user,startval)
        out = clean_output1(out,Y,iT,length(Y),extout,dataout,ivItemcat)
      } else if(is.null(id_high)&!is.null(Z)&is.null(Zh)){
        out = LCA_fast_init_wcov_poly(mY,mZ,iT,ivItemcat,incomplete,kmea,maxIter,tol,fixed,reord,reord_user,
                                      NRtol,NRmaxit,verbose,startval)
        out = clean_output2(out,Y,iT,c("Intercept",Z),length(Y),length(Z)+1,extout,dataout,ivItemcat)
      } else if(!is.null(id_high)&is.null(Z)&is.null(Zh)){
        init  = meas_Init_poly(mY,id_high,vNj,iM,iT,ivItemcat,incomplete,kmea,maxIter,tol,reord,reord_user,startval)
        if(!incomplete){
          out   = MLTLCA_poly(mY,vNj,init$vOmega_start,init$mPi_start,init$mPhi_start,
                              ivItemcat,init$first_poly,reord_user,reord_user_high,maxIter,tol,reord)
        } else{
          mY[is.na(mY)] = max(na.omit(mY))+1
          out   = MLTLCA_poly_includeall(mY,mDesign,vNj,init$vOmega_start,init$mPi_start,init$mPhi_start,
                                         ivItemcat,init$first_poly,reord_user,reord_user_high,maxIter,tol,reord)
        }
        out = clean_output3(out,Y,iT,iM,mY,id_high,length(Y),id_high_levs,id_high_name,extout,dataout,ivItemcat)
      } else if(!is.null(id_high)&!is.null(Z)&is.null(Zh)){
        init1 = meas_Init_poly(mY,id_high,vNj,iM,iT,ivItemcat,incomplete,kmea,maxIter,tol,reord,reord_user,startval)
        if(!incomplete){
          init2 = MLTLCA_poly(mY,vNj,init1$vOmega_start,init1$mPi_start,init1$mPhi_start,
                              ivItemcat,init1$first_poly,reord_user,reord_user_high,maxIter,tol,reord)
        } else{
          mY[is.na(mY)] = max(na.omit(mY))+1
          init2 = MLTLCA_poly_includeall(mY,mDesign,vNj,init1$vOmega_start,init1$mPi_start,init1$mPhi_start,
                                         ivItemcat,init1$first_poly,reord_user,reord_user_high,maxIter,tol,reord)
        }
        P             = ncol(mZ)
        vOmega_start  = init2$vOmega
        cGamma_start  = array(c(rbind(init2$mGamma,matrix(0,(iT-1)*(P-1),iM))),c(iT-1,P,iM))
        mPhi_start    = init2$mPhi
        mStep1Var     = init2$Varmat
        #
        mY      = mY[complete.cases(mZ),]
        mDesign = mDesign[complete.cases(mZ),]
        id_high = id_high[complete.cases(mZ)]
        vNj     = table(id_high)
        mZ      = mZ[complete.cases(mZ),]
        #
        if(verbose){
          if(fixed==1|fixed==2){
            pb=txtProgressBar(char="Fitting structural model...",width=1)
            setTxtProgressBar(pb,1)
            close(pb)
          }
        }
        #
        if(!incomplete){
          if(!fixedslopes){
            out = MLTLCA_cov_poly(mY,mZ,vNj,vOmega_start,cGamma_start,mPhi_start,mStep1Var,
                                  ivItemcat,maxIter,tol,fixedpars,1,NRtol,NRmaxit)
          } else{
            out = MLTLCA_covWfixed_poly(mY,mZ,vNj,vOmega_start,cGamma_start,mPhi_start,mStep1Var,
                                        ivItemcat,maxIter,tol,fixedpars,1,NRtol,NRmaxit)
          }
        } else{
          if(!fixedslopes){
            out = MLTLCA_cov_poly_includeall(mY,mDesign,mZ,vNj,vOmega_start,cGamma_start,mPhi_start,mStep1Var,
                                             ivItemcat,maxIter,tol,fixedpars,1,NRtol,NRmaxit)
          } else{
            out = MLTLCA_covWfixed_poly_includeall(mY,mDesign,mZ,vNj,vOmega_start,cGamma_start,mPhi_start,mStep1Var,
                                                   ivItemcat,maxIter,tol,fixedpars,1,NRtol,NRmaxit)
          }
        }
        out = clean_output4(out,Y,iT,iM,c("Intercept",Z),mY,mZ,id_high,length(Y),P,id_high_levs,id_high_name,extout,dataout,ivItemcat)
      } else if(!is.null(id_high)&!is.null(Z)&!is.null(Zh)){
        init1 = meas_Init_poly(mY,id_high,vNj,iM,iT,ivItemcat,incomplete,kmea,maxIter,tol,reord,reord_user,startval)
        if(!incomplete){
          init2 = MLTLCA_poly(mY,vNj,init1$vOmega_start,init1$mPi_start,init1$mPhi_start,
                              ivItemcat,init1$first_poly,reord_user,reord_user_high,maxIter,tol,reord)
        } else{
          mY[is.na(mY)] = max(na.omit(mY))+1
          init2 = MLTLCA_poly_includeall(mY,mDesign,vNj,init1$vOmega_start,init1$mPi_start,init1$mPhi_start,
                                         ivItemcat,init1$first_poly,reord_user,reord_user_high,maxIter,tol,reord)
        }
        P                 = ncol(mZ)
        P_high            = ncol(mZh)
        cGamma_start      = array(c(rbind(init2$mGamma,matrix(0,(iT-1)*(P-1),iM))),c(iT-1,P,iM))
        mPhi_start        = init2$mPhi
        mStep1Var         = init2$Varmat
        mDelta_start      = matrix(0,iM-1,P_high)
        mDelta_start[,1]  = init2$vDelta
        #
        nomissing = complete.cases(mZ)&complete.cases(mZh)
        mY        = mY[nomissing,]
        mDesign   = mDesign[nomissing,]
        id_high   = id_high[nomissing]
        vNj       = table(id_high)
        mZ        = mZ[nomissing,]
        mZh       = mZh[nomissing,]
        mZh       = mZh[!duplicated(id_high),]
        #
        if(verbose){
          if(fixed==1|fixed==2){
            pb=txtProgressBar(char="Fitting structural model...",width=1)
            setTxtProgressBar(pb,1)
            close(pb)
          }
        }
        #
        if(!incomplete){
          if(!fixedslopes){
            out = MLTLCA_covlowhigh_poly(mY,mZ,mZh,vNj,mDelta_start,cGamma_start,mPhi_start,
                                         mStep1Var,ivItemcat,maxIter,tol,fixedpars,1,NRtol,NRmaxit)
          } else{
            out = MLTLCA_covWfixedlowhigh_poly(mY,mZ,mZh,vNj,mDelta_start,cGamma_start,mPhi_start,
                                               mStep1Var,ivItemcat,maxIter,tol,fixedpars,1,NRtol,NRmaxit)
          }
        } else{
          if(!fixedslopes){
            out = MLTLCA_covlowhigh_poly_includeall(mY,mDesign,mZ,mZh,vNj,mDelta_start,cGamma_start,mPhi_start,
                                                    mStep1Var,ivItemcat,maxIter,tol,fixedpars,1,NRtol,NRmaxit)
          } else{
            out = MLTLCA_covWfixedlowhigh_poly_includeall(mY,mDesign,mZ,mZh,vNj,mDelta_start,cGamma_start,mPhi_start,
                                                          mStep1Var,ivItemcat,maxIter,tol,fixedpars,1,NRtol,NRmaxit)
          }
        }
        out = clean_output5(out,Y,iT,iM,c("Intercept",Z),c("Intercept",Zh),mY,mZ,mZh,id_high,length(Y),P,P_high,id_high_levs,id_high_name,extout,dataout,ivItemcat)
      }
    } else{
      if(is.null(id_high)&is.null(Z)&is.null(Zh)){
        out = LCA_fast_init(mY,iT,kmea,maxIter,tol,reord,reord_user,startval)
        out = clean_output1(out,Y,iT,length(Y),extout,dataout,ivItemcat)
      } else if(is.null(id_high)&!is.null(Z)&is.null(Zh)){
        out = LCA_fast_init_wcov(mY,mZ,iT,kmea,maxIter,tol,fixed,reord,reord_user,
                                 NRtol,NRmaxit,verbose,startval)
        out = clean_output2(out,Y,iT,c("Intercept",Z),length(Y),length(Z)+1,extout,dataout,ivItemcat)
      } else if(!is.null(id_high)&is.null(Z)&is.null(Zh)){
        init  = meas_Init(mY,id_high,vNj,iM,iT,kmea,maxIter,tol,reord,reord_user,startval)
        out   = MLTLCA(mY,vNj,init$vOmega_start,init$mPi_start,init$mPhi_start,reord_user,reord_user_high,
                       maxIter,tol,reord)
        out = clean_output3(out,Y,iT,iM,mY,id_high,length(Y),id_high_levs,id_high_name,extout,dataout,ivItemcat)
      } else if(!is.null(id_high)&!is.null(Z)&is.null(Zh)){
        init1 = meas_Init(mY,id_high,vNj,iM,iT,kmea,maxIter,tol,reord,reord_user,startval)
        init2 = MLTLCA(mY,vNj,init1$vOmega_start,init1$mPi_start,init1$mPhi_start,reord_user,reord_user_high,
                       maxIter,tol,reord)
        P             = ncol(mZ)
        vOmega_start  = init2$vOmega
        cGamma_start  = array(c(rbind(init2$mGamma,matrix(0,(iT-1)*(P-1),iM))),c(iT-1,P,iM))
        mPhi_start    = init2$mPhi
        mStep1Var     = init2$Varmat
        #
        mY      = mY[complete.cases(mZ),]
        mDesign = mDesign[complete.cases(mZ),]
        id_high = id_high[complete.cases(mZ)]
        vNj     = table(id_high)
        mZ      = mZ[complete.cases(mZ),]
        #
        if(verbose){
          if(fixed==1|fixed==2){
            pb=txtProgressBar(char="Fitting structural model...",width=1)
            setTxtProgressBar(pb,1)
            close(pb)
          }
        }
        #
        out = MLTLCA_cov(mY,mZ,vNj,vOmega_start,cGamma_start,mPhi_start,mStep1Var,
                         maxIter,tol,fixedpars,NRtol,NRmaxit)
        out = clean_output4(out,Y,iT,iM,c("Intercept",Z),mY,mZ,id_high,length(Y),P,id_high_levs,id_high_name,extout,dataout,ivItemcat)
      } else if(!is.null(id_high)&!is.null(Z)&!is.null(Zh)){
        init1   = meas_Init(mY,id_high,vNj,iM,iT,kmea,maxIter,tol,reord,reord_user,startval)
        init2   = MLTLCA(mY,vNj,init1$vOmega_start,init1$mPi_start,init1$mPhi_start,reord_user,reord_user_high,
                         maxIter,tol,reord)
        P                 = ncol(mZ)
        P_high            = ncol(mZh)
        cGamma_start      = array(c(rbind(init2$mGamma,matrix(0,(iT-1)*(P-1),iM))),c(iT-1,P,iM))
        mPhi_start        = init2$mPhi
        mStep1Var         = init2$Varmat
        mDelta_start      = matrix(0,iM-1,P_high)
        mDelta_start[,1]  = init2$vDelta
        #
        nomissing = complete.cases(mZ)&complete.cases(mZh)
        mY        = mY[nomissing,]
        mDesign   = mDesign[nomissing,]
        id_high   = id_high[nomissing]
        vNj       = table(id_high)
        mZ        = mZ[nomissing,]
        mZh       = mZh[nomissing,]
        mZh       = mZh[!duplicated(id_high),]
        #
        if(verbose){
          if(fixed==1|fixed==2){
            pb=txtProgressBar(char="Fitting structural model...",width=1)
            setTxtProgressBar(pb,1)
            close(pb)
          }
        }
        #
        out = MLTLCA_covlowhigh(mY,mZ,mZh,vNj,mDelta_start,cGamma_start,mPhi_start,
                                mStep1Var,maxIter,tol,fixedpars,NRtol,NRmaxit)
        out = clean_output5(out,Y,iT,iM,c("Intercept",Z),c("Intercept",Zh),mY,mZ,mZh,id_high,length(Y),P,P_high,id_high_levs,id_high_name,extout,dataout,ivItemcat)
      }
    }
  } else if(approach=="model selection on low"){
    if(any(ivItemcat>2)|incomplete){
      out = sel_other_poly(mY,mDesign,id_high,iT,iM,ivItemcat,approach,verbose,kmea,maxIter,tol,reord)
      out_mat = matrix(round(unlist(out[3:7]),2),length(iT),5,
                       dimnames=list(paste0("iT=",iT),
                                     c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
      out_mat[out_mat=="Inf"|out_mat=="-Inf"|is.na(out_mat)] = "-"
      optimal = matrix(out$iT_best,dimnames=list("iT=",""))
      list_sel = list(model_selection=out_mat,optimal=optimal)
      if(optimal==1){
        if(verbose)print(noquote(list_sel))
        stop("iT=1 optimal.",call.=FALSE)
      }
      out = clean_output1(out$outmuLCA,Y,out$iT_best,length(Y),extout,dataout,ivItemcat)
      out$model_selection = list_sel
      if(verbose)print(noquote(list_sel))
    } else{
      out = sel_other(mY,id_high,iT,iM,approach,verbose,kmea,maxIter,tol,reord)
      out_mat = matrix(round(unlist(out[3:7]),2),length(iT),5,
                       dimnames=list(paste0("iT=",iT),
                                     c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
      out_mat[out_mat=="Inf"|out_mat=="-Inf"|is.na(out_mat)] = "-"
      optimal = matrix(out$iT_best,dimnames=list("iT=",""))
      list_sel = list(model_selection=out_mat,optimal=optimal)
      if(optimal==1){
        if(verbose)print(noquote(list_sel))
        stop("iT=1 optimal.",call.=FALSE)
      }
      out = clean_output1(out$outmuLCA,Y,out$iT_best,length(Y),extout,dataout,ivItemcat)
      out$model_selection = list_sel
      if(verbose)print(noquote(list_sel))
    }
  } else if(approach=="model selection on low with high"){
    if(any(ivItemcat>2)|incomplete){
      out = sel_other_poly(mY,mDesign,id_high,iT,iM,ivItemcat,approach,verbose,kmea,maxIter,tol,reord)
      out_mat = matrix(round(unlist(out[3:7]),2),length(iT),5,
                       dimnames=list(paste0("iT=",iT,",iM=",iM),
                                     c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
      out_mat[out_mat=="Inf"|out_mat=="-Inf"|is.na(out_mat)] = "-"
      optimal = matrix(out$iT_best,dimnames=list("iT=",""))
      list_sel = list(model_selection=out_mat,optimal=optimal)
      if(optimal==1){
        if(verbose)print(noquote(list_sel))
        stop("iT=1 optimal.",call.=FALSE)
      }
      out = clean_output3(out$outmuLCA,Y,out$iT_best,iM,mY,id_high,length(Y),id_high_levs,id_high_name,extout,dataout,ivItemcat)
      out$model_selection = list_sel
      if(verbose)print(noquote(list_sel))
    } else{
      out = sel_other(mY,id_high,iT,iM,approach,verbose,kmea,maxIter,tol,reord)
      out_mat = matrix(round(unlist(out[3:7]),2),length(iT),5,
                       dimnames=list(paste0("iT=",iT,",iM=",iM),
                                     c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
      out_mat[out_mat=="Inf"|out_mat=="-Inf"|is.na(out_mat)] = "-"
      optimal = matrix(out$iT_best,dimnames=list("iT=",""))
      list_sel = list(model_selection=out_mat,optimal=optimal)
      if(optimal==1){
        if(verbose)print(noquote(list_sel))
        stop("iT=1 optimal.",call.=FALSE)
      }
      out = clean_output3(out$outmuLCA,Y,out$iT_best,iM,mY,id_high,length(Y),id_high_levs,id_high_name,extout,dataout,ivItemcat)
      out$model_selection = list_sel
      if(verbose)print(noquote(list_sel))
    }
  } else if(approach=="model selection on high"){
    if(any(ivItemcat>2)|incomplete){
      out = sel_other_poly(mY,mDesign,id_high,iT,iM,ivItemcat,approach,verbose,kmea,maxIter,tol,reord)
      out_mat = matrix(round(unlist(out[3:7]),2),length(iM),5,
                       dimnames=list(paste0("iT=",iT,",iM=",iM),
                                     c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
      out_mat[out_mat=="Inf"|out_mat=="-Inf"|is.na(out_mat)] = "-"
      optimal = matrix(out$iM_best,dimnames=list("iM=",""))
      list_sel = list(model_selection=out_mat,optimal=optimal)
      if(out$iM_best>1){
        out = clean_output3(out$outmuLCA,Y,iT,out$iM_best,mY,id_high,length(Y),id_high_levs,id_high_name,extout,dataout,ivItemcat)
      } else if(out$iM_best==1){
        out = clean_output1(out$outmuLCA,Y,iT,length(Y),extout,dataout,ivItemcat)
      }
      out$model_selection = list_sel
      if(verbose)print(noquote(list_sel))
    } else{
      out = sel_other(mY,id_high,iT,iM,approach,verbose,kmea,maxIter,tol,reord)
      out_mat = matrix(round(unlist(out[3:7]),2),length(iM),5,
                       dimnames=list(paste0("iT=",iT,",iM=",iM),
                                     c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
      out_mat[out_mat=="Inf"|out_mat=="-Inf"|is.na(out_mat)] = "-"
      optimal = matrix(out$iM_best,dimnames=list("iM=",""))
      list_sel = list(model_selection=out_mat,optimal=optimal)
      if(out$iM_best>1){
        out = clean_output3(out$outmuLCA,Y,iT,out$iM_best,mY,id_high,length(Y),id_high_levs,id_high_name,extout,dataout,ivItemcat)
      } else if(out$iM_best==1){
        out = clean_output1(out$outmuLCA,Y,iT,length(Y),extout,dataout,ivItemcat)
      }
      out$model_selection = list_sel
      if(verbose)print(noquote(list_sel))
    }
  } else if(approach=="model selection on low and high"){
    if(sequential){
      if(any(ivItemcat>2)|incomplete){
        out = lukosel_fun_poly(mY,mDesign,id_high,iT,iM,ivItemcat,verbose,kmea,maxIter,tol,reord)
        step1mat = matrix(as.character(round(unlist(out[4:8]),2)),length(iT),5,
                          dimnames=list(paste0("iT=",iT),
                                        c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
        step1mat[step1mat=="Inf"|step1mat=="-Inf"|is.na(step1mat)] = "-"
        step2mat = matrix(as.character(round(unlist(out[9:13]),2)),length(iM),5,
                          dimnames=list(paste0("iT*,","iM=",iM),
                                        c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
        step2mat[step2mat=="Inf"|step2mat=="-Inf"|is.na(step2mat)] = "-"
        step3mat = matrix(as.character(round(unlist(out[14:18]),2)),length(iT),5,
                          dimnames=list(paste0("iT=",iT,",iM*"),
                                        c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
        step3mat[step3mat=="Inf"|step3mat=="-Inf"|is.na(step3mat)] = "-"
        optimal = matrix(c(out$iT_opt,out$iM_opt),dimnames=list(c("iT=","iM="),""))
        list_luk = list(step1 = step1mat,step2 = step2mat,step3 = step3mat, optimal = optimal)
        if(out$iM_opt>1){
          out = clean_output3(out$outmuLCA_step3,Y,out$iT_opt,out$iM_opt,mY,id_high,length(Y),
                              id_high_levs,id_high_name,extout,dataout,ivItemcat)
        } else if(out$iM_opt==1){
          out = clean_output1(out$outmuLCA_step3,Y,out$iT_opt,length(Y),extout,dataout,ivItemcat)
        }
        out$model_selection = list_luk
        if(verbose)print(noquote(list_luk))
      } else{
        out = lukosel_fun(mY,id_high,iT,iM,verbose,kmea,maxIter,tol,reord)
        step1mat = matrix(as.character(round(unlist(out[4:8]),2)),length(iT),5,
                          dimnames=list(paste0("iT=",iT),
                                        c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
        step1mat[step1mat=="Inf"|step1mat=="-Inf"|is.na(step1mat)] = "-"
        step2mat = matrix(as.character(round(unlist(out[9:13]),2)),length(iM),5,
                          dimnames=list(paste0("iT*,","iM=",iM),
                                        c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
        step2mat[step2mat=="Inf"|step2mat=="-Inf"|is.na(step2mat)] = "-"
        step3mat = matrix(as.character(round(unlist(out[14:18]),2)),length(iT),5,
                          dimnames=list(paste0("iT=",iT,",iM*"),
                                        c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
        step3mat[step3mat=="Inf"|step3mat=="-Inf"|is.na(step3mat)] = "-"
        optimal = matrix(c(out$iT_opt,out$iM_opt),dimnames=list(c("iT=","iM="),""))
        list_luk = list(step1 = step1mat,step2 = step2mat,step3 = step3mat, optimal = optimal)
        if(out$iM_opt>1){
          out = clean_output3(out$outmuLCA_step3,Y,out$iT_opt,out$iM_opt,mY,id_high,length(Y),
                              id_high_levs,id_high_name,extout,dataout,ivItemcat)
        } else if(out$iM_opt==1){
          out = clean_output1(out$outmuLCA_step3,Y,out$iT_opt,length(Y),extout,dataout,ivItemcat)
        }
        out$model_selection = list_luk
        if(verbose)print(noquote(list_luk))
      }
    } else{
      #
      if(verbose){
        pb=txtProgressBar(char="Fitting measurement models simultaneously...",width=1)
        setTxtProgressBar(pb,1)
        close(pb)
      }
      #
      if(any(ivItemcat>2)|incomplete){
        mYfoo       = mY
        mDesignfoo  = mDesign
        id_highfoo  = id_high
        ivItemcatfoo = ivItemcat
        select_mat  = as.matrix(tidyr::expand_grid(iT,iM))
        nmod        = nrow(select_mat)
        iV          = parallel::detectCores()-numFreeCores
        cluster     = parallel::makeCluster(iV)
        parallel::clusterExport(cluster,c("mYfoo","mDesignfoo","id_highfoo","ivItemcatfoo","select_mat","nmod"), 
                                envir = environment())
        parallel::clusterExport(cluster, 
                                unclass(lsf.str(envir = asNamespace("multilevLCA"), 
                                                all = T)),
                                envir = as.environment(asNamespace("multilevLCA"))
        )
        parallel::clusterEvalQ(cluster, {
          library(clustMixType)
          library(multilevLCA)
          library(mclust)
          library(tictoc)
          library(numDeriv)
          library(klaR)
          library(tidyverse)
          library(MASS)})
        simultaneous_out = parallel::parLapply(cluster,1:nmod,
                                               function(x){
                                                 iFoo       = select_mat[x,]
                                                 iT_curr    = iFoo[1]
                                                 iM_curr    = iFoo[2]
                                                 out_simult = simultsel_fun_poly(mYfoo,mDesignfoo,id_highfoo,iT_curr,iM_curr,ivItemcatfoo,kmea,maxIter,tol,reord)
                                                 BIClow     = out_simult$BIClow
                                                 BIChigh    = out_simult$BIChigh
                                                 AIC        = out_simult$AIC
                                                 ICLlow     = out_simult$ICL_BIClow
                                                 ICLhigh    = out_simult$ICL_BIChigh
                                                 nout       = c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")
                                                 out        = c(BIClow,BIChigh,AIC,ICLlow,ICLhigh)
                                                 names(out) = nout
                                                 return(list(out=out,fit=out_simult$modFIT))
                                               }
        )
        parallel::stopCluster(cluster)
        mod_sel = apply(t(sapply(simultaneous_out,function(x){x$out})),
                        2,
                        function(x){round(as.numeric(x),2)})
        rownames(mod_sel) = paste0("iT=",select_mat[,1],",iM=",select_mat[,2])
        r_optimal = which.min(mod_sel[,1])
        optimal   = matrix(c(select_mat[r_optimal,1],select_mat[r_optimal,2]),dimnames=list(c("iT=","iM="),""))
        mod_sel[mod_sel==Inf|mod_sel==-Inf|is.na(mod_sel)] = "-"
        list_simul = list(model_selection=mod_sel,optimal=optimal)
        if(select_mat[r_optimal,1]==1){
          if(verbose)print(noquote(list_simul))
          stop("iT=1 optimal.",call.=FALSE)
        } else if(select_mat[r_optimal,2]>1){
          out = clean_output3(simultaneous_out[[r_optimal]]$fit,
                              Y,select_mat[r_optimal,1],select_mat[r_optimal,2],
                              mY,id_high,length(Y),id_high_levs,id_high_name,extout,dataout,ivItemcat)
        } else if(select_mat[r_optimal,2]==1){
          out = clean_output1(simultaneous_out[[r_optimal]]$fit,
                              Y,select_mat[r_optimal,1],length(Y),extout,dataout,ivItemcat)
        }
        out$model_selection = list_simul
        if(verbose)print(noquote(list_simul))
      } else{
        mYfoo       = mY
        id_highfoo  = id_high
        select_mat  = as.matrix(tidyr::expand_grid(iT,iM))
        nmod        = nrow(select_mat)
        iV          = parallel::detectCores()-numFreeCores
        cluster     = parallel::makeCluster(iV)
        parallel::clusterExport(cluster,c("mYfoo","id_highfoo","select_mat","nmod"), envir = environment())
        parallel::clusterExport(cluster, 
                                unclass(lsf.str(envir = asNamespace("multilevLCA"), 
                                                all = T)),
                                envir = as.environment(asNamespace("multilevLCA"))
        )
        parallel::clusterEvalQ(cluster, {
          library(clustMixType)
          library(multilevLCA)
          library(mclust)
          library(tictoc)
          library(numDeriv)
          library(klaR)
          library(tidyverse)
          library(MASS)})
        simultaneous_out = parallel::parLapply(cluster,1:nmod,
                                               function(x){
                                                 iFoo       = select_mat[x,]
                                                 iT_curr    = iFoo[1]
                                                 iM_curr    = iFoo[2]
                                                 out_simult = simultsel_fun(mYfoo,id_highfoo,iT_curr,iM_curr,kmea,maxIter,tol,reord)
                                                 BIClow     = out_simult$BIClow
                                                 BIChigh    = out_simult$BIChigh
                                                 AIC        = out_simult$AIC
                                                 ICLlow     = out_simult$ICL_BIClow
                                                 ICLhigh    = out_simult$ICL_BIChigh
                                                 nout       = c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")
                                                 out        = c(BIClow,BIChigh,AIC,ICLlow,ICLhigh)
                                                 names(out) = nout
                                                 return(list(out=out,fit=out_simult$modFIT))
                                               }
        )
        parallel::stopCluster(cluster)
        mod_sel = apply(t(sapply(simultaneous_out,function(x){x$out})),
                        2,
                        function(x){round(as.numeric(x),2)})
        rownames(mod_sel) = paste0("iT=",select_mat[,1],",iM=",select_mat[,2])
        r_optimal = which.min(mod_sel[,1])
        optimal   = matrix(c(select_mat[r_optimal,1],select_mat[r_optimal,2]),dimnames=list(c("iT=","iM="),""))
        mod_sel[mod_sel==Inf|mod_sel==-Inf|is.na(mod_sel)] = "-"
        list_simul = list(model_selection=mod_sel,optimal=optimal)
        if(select_mat[r_optimal,1]==1){
          if(verbose)print(noquote(list_simul))
          stop("iT=1 optimal.",call.=FALSE)
        } else if(select_mat[r_optimal,2]>1){
          out = clean_output3(simultaneous_out[[r_optimal]]$fit,
                              Y,select_mat[r_optimal,1],select_mat[r_optimal,2],
                              mY,id_high,length(Y),id_high_levs,id_high_name,extout,dataout,ivItemcat)
        } else if(select_mat[r_optimal,2]==1){
          out = clean_output1(simultaneous_out[[r_optimal]]$fit,
                              Y,select_mat[r_optimal,1],length(Y),extout,dataout,ivItemcat)
        }
        out$model_selection = list_simul
        if(verbose)print(noquote(list_simul))
      }
    }
  }
  
  if(approach=="direct"&!is.null(Z)){
    if(fixed==0){
      out$estimator = as.matrix("One-step")
      rownames(out$estimator) = colnames(out$estimator) = ""
    } else if(fixed==1){
      out$estimator = as.matrix("Two-step")
      rownames(out$estimator) = colnames(out$estimator) = ""
    } else if(fixed==2){
      out$estimator = as.matrix("Two-stage")
      rownames(out$estimator) = colnames(out$estimator) = ""
    }
  }
  if(nrow(data)==nrow_data){
    out$missing_values = as.matrix("None")
    rownames(out$missing_values) = colnames(out$missing_values) = ""
  } else if(!incomplete){
    out$missing_values = as.matrix("Row-wise deletion")
    rownames(out$missing_values) = colnames(out$missing_values) = ""
  } else if(incomplete){
    out$missing_values = as.matrix("Kept via FIML estimation")
    rownames(out$missing_values) = colnames(out$missing_values) = ""
  }
  if(approach=="direct"&!is.null(Z)){
    if(is.null(Zh)){
      out$sample_size = as.matrix(c(nrow(data),nrow(na.omit(mZ))))
      rownames(out$sample_size) = c("Measurement model:","Structural model:")
      colnames(out$sample_size) = ""
    } else if(!is.null(id_high)&!is.null(Zh)){
      out$sample_size = as.matrix(c(nrow(data),nrow(mZ),nrow(mZh)))
      rownames(out$sample_size) = c("Measurement model:",
                                    "Lower-level structural model:",
                                    "Higher-level structural model:")
      colnames(out$sample_size) = ""
    }
  } else{
    out$sample_size = as.matrix(nrow(data))
    rownames(out$sample_size) = colnames(out$sample_size) = ""
  }
  if(out$spec == "Single-level LC model with covariates"){
    gammas  = as.matrix(out$mGamma)
    SE      = as.matrix(out$SEs_cor_gamma)
    Zscore  = as.matrix(out$mGamma)/as.matrix(out$SEs_cor_gamma)
    pval    = 2*(1 - pnorm(abs(as.matrix(out$mGamma)/as.matrix(out$SEs_cor_gamma))))
    stru_inference = array(NA,c(nrow(gammas),4,ncol(gammas)),
                           list(paste0(substr(rownames(as.matrix(out$mGamma)), 1, nchar(rownames(as.matrix(out$mGamma))) - 1), ")"),c("Gamma", "S.E.", "Z-score", "p-value"),paste0("C",1+1:ncol(gammas))))
    for (i in 1:ncol(gammas)){
      stru_inference[,,i] = cbind(gammas[,i], SE[,i], Zscore[,i], pval[,i])
    }
    out$stru_inference = stru_inference
  } else if(out$spec == "Multilevel LC model with lower-level covariates"){
    if(is.null(out$mGamma_fixslope)){
      stru_inference = list()
      for (i in 1:dim(out$cGamma)[3]){
        gammas  = as.matrix(out$cGamma[,,i])
        SE      = as.matrix(out$SEs_cor_gamma[,,i])
        Zscore  = as.matrix(out$cGamma[,,i])/as.matrix(out$SEs_cor_gamma[,,i])
        pval    = 2*(1 - pnorm(abs(as.matrix(out$cGamma[,,i])/as.matrix(out$SEs_cor_gamma[,,i]))))
        stru_inference[[paste0("G",i)]] = array(NA,c(nrow(gammas),4,ncol(gammas)),
                                                list(paste0(substr(rownames(as.matrix(out$cGamma[,,i])), 1, nchar(rownames(as.matrix(out$cGamma[,,i]))) - 1), ",G", i, ")"),c("Gamma", "S.E.", "Z-score", "p-value"),paste0("C",1+1:ncol(gammas))))
        for (j in 1:ncol(gammas)){
          stru_inference[[paste0("G",i)]][,,j] = cbind(gammas[,j], SE[,j], Zscore[,j], pval[,j])
        }
      }
    } else{
      gammas  = as.matrix(out$mGamma_fixslope)
      SE      = as.matrix(out$SEs_cor_gamma)
      Zscore  = as.matrix(out$mGamma_fixslope)/as.matrix(out$SEs_cor_gamma)
      pval    = 2*(1 - pnorm(abs(as.matrix(out$mGamma_fixslope)/as.matrix(out$SEs_cor_gamma))))
      iM = nrow(out$vOmega)
      stru_inference = array(NA,c(nrow(gammas),4,ncol(gammas)),
                             list(c(paste0(substr(rownames(as.matrix(out$mGamma_fixslope))[1:iM], 1, 17), substr(rownames(as.matrix(out$mGamma_fixslope))[1:iM], 18, nchar(rownames(as.matrix(out$mGamma_fixslope))[1:iM]))),
                                    paste0(substr(rownames(as.matrix(out$mGamma_fixslope))[-c(1:iM)], 1, nchar(rownames(as.matrix(out$mGamma_fixslope))[-c(1:iM)]) - 1), ")")),
                                  c("Gamma", "S.E.", "Z-score", "p-value"),
                                  paste0("C",1+1:ncol(gammas))))
      for (i in 1:ncol(gammas)){
        stru_inference[,,i] = cbind(gammas[,i], SE[,i], Zscore[,i], pval[,i])
      }
    }
    out$stru_inference = stru_inference
  } else if(out$spec == "Multilevel LC model with lower- and higher-level covariates"){
    stru_inference = list()
    alphas  = as.matrix(out$mAlpha)
    SE      = as.matrix(out$SEs_cor_alpha)
    Zscore  = as.matrix(out$mAlpha)/as.matrix(out$SEs_cor_alpha)
    pval    = 2*(1 - pnorm(abs(as.matrix(out$mAlpha)/as.matrix(out$SEs_cor_alpha))))
    stru_inference$higher_level = array(NA,c(nrow(alphas),4,ncol(alphas)),
                                        list(paste0(substr(rownames(as.matrix(out$mAlpha)), 1, nchar(rownames(as.matrix(out$mAlpha))) - 1), ")"),c("Alpha", "S.E.", "Z-score", "p-value"),paste0("G",1+1:ncol(alphas))))
    for (i in 1:ncol(alphas)){
      stru_inference$higher_level[,,i] = noquote(cbind(alphas[,i], SE[,i], Zscore[,i], pval[,i]))
    }
    if(is.null(out$mGamma_fixslope)){
      stru_inference$lower_level = list()
      for (i in 1:dim(out$cGamma)[3]){
        gammas  = as.matrix(out$cGamma[,,i])
        SE      = as.matrix(out$SEs_cor_gamma[,,i])
        Zscore  = as.matrix(out$cGamma[,,i])/as.matrix(out$SEs_cor_gamma[,,i])
        pval    = 2*(1 - pnorm(abs(as.matrix(out$cGamma[,,i])/as.matrix(out$SEs_cor_gamma[,,i]))))
        stru_inference$lower_level[[paste0("G",i)]] = array(NA,c(nrow(gammas),4,ncol(gammas)),
                                                            list(paste0(substr(rownames(as.matrix(out$cGamma[,,i])), 1, nchar(rownames(as.matrix(out$cGamma[,,i]))) - 1), ",G", i, ")"),c("Gamma", "S.E.", "Z-score", "p-value"),paste0("C",1+1:ncol(gammas))))
        for (j in 1:ncol(gammas)){
          stru_inference$lower_level[[paste0("G",i)]][,,j] = cbind(gammas[,j], SE[,j], Zscore[,j], pval[,j])
        }
      }
    } else{
      gammas  = as.matrix(out$mGamma_fixslope)
      SE      = as.matrix(out$SEs_cor_gamma)
      Zscore  = as.matrix(out$mGamma_fixslope)/as.matrix(out$SEs_cor_gamma)
      pval    = 2*(1 - pnorm(abs(as.matrix(out$mGamma_fixslope)/as.matrix(out$SEs_cor_gamma))))
      iM = nrow(out$vOmega)
      stru_inference$lower_level = array(NA,c(nrow(gammas),4,ncol(gammas)),
                                         list(c(paste0(substr(rownames(as.matrix(out$mGamma_fixslope))[1:iM], 1, 17), substr(rownames(as.matrix(out$mGamma_fixslope))[1:iM], 18, nchar(rownames(as.matrix(out$mGamma_fixslope))[1:iM]))),
                                                paste0(substr(rownames(as.matrix(out$mGamma_fixslope))[-c(1:iM)], 1, nchar(rownames(as.matrix(out$mGamma_fixslope))[-c(1:iM)]) - 1), ")")),
                                              c("Gamma", "S.E.", "Z-score", "p-value"),
                                              paste0("C",1+1:ncol(gammas))))
      for (i in 1:ncol(gammas)){
        stru_inference$lower_level[,,i] = cbind(gammas[,i], SE[,i], Zscore[,i], pval[,i])
      }
    }
    out$stru_inference = stru_inference
  }

  out$call = match.call()
  class(out) = "multiLCA"
  return(out)
}
#
LCA_fast_init = function(mY, iT, kmea = T, maxIter = 1e3, tol = 1e-8, reord = 1, reord_user = 1:iT-1, startval = NULL){
  group_by_all  = NULL
  mY_df         = data.frame(mY)
  mY_aggr       = as.matrix(mY_df%>%group_by_all%>%count)
  iHf           = dim(mY_aggr)[2]
  freq          = mY_aggr[,iHf]
  mY_unique     = mY_aggr[,-iHf]
  if(is.null(startval)&kmea==FALSE){
    clusfoo     = klaR::kmodes(mY_unique,modes=iT,iter.max=100)$cluster
    mU          = vecTomatClass(clusfoo)
  } else if(is.null(startval)&kmea==TRUE){
    prscores    = prcomp(mY_unique)
    num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
    prscores    = prscores$x[,1:num_pc]
    spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
    clusfoo     = spectclust$cluster
    mU          = vecTomatClass(clusfoo)
  } else{
    clusfoo     = vecTomatClass(factor(startval))
    nclusfoo    = ncol(clusfoo)
    colnames(clusfoo) = paste0("clusfoo",1:nclusfoo)
    Y           = names(mY_df)
    mY_df       = data.frame(cbind(mY,clusfoo))
    clusfoo     = as.matrix(mY_df%>%group_by(across(all_of(Y)))%>%summarise_at(vars(paste0("clusfoo",1:nclusfoo)),mean))
    clusfoo     = clusfoo[,paste0("clusfoo",1:nclusfoo)]
    mU          = clusfoo
  }
  out           = LCA_fast(mY_unique,freq,iT,mU,reord_user,maxIter,tol,reord)
  out$mY_unique = mY_unique
  out$freq      = freq
  return(out)
}
#
LCA_fast_init_wcov = function(mY, mZ, iT, kmea = T, maxIter = 1e3, tol = 1e-8, fixed = 0, reord = 1, reord_user = 1:iT-1,
                              NRtol = 1e-6, NRmaxit = 100, verbose, startval){
  # mZ must include a column of ones!
  group_by_all    = NULL
  mY_df           = data.frame(mY)
  mY_aggr         = as.matrix(mY_df%>%group_by_all%>%count)
  iHf             = dim(mY_aggr)[2]
  freq            = mY_aggr[,iHf]
  mY_unique       = mY_aggr[,-iHf]
  if(is.null(startval)&kmea==FALSE){
    clusfoo     = klaR::kmodes(mY_unique,modes=iT,iter.max=100)$cluster
    mU          = vecTomatClass(clusfoo)
  } else if(is.null(startval)&kmea==TRUE){
    prscores    = prcomp(mY_unique)
    num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
    prscores    = prscores$x[,1:num_pc]
    spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
    clusfoo     = spectclust$cluster
    mU          = vecTomatClass(clusfoo)
  } else{
    clusfoo     = vecTomatClass(factor(startval))
    nclusfoo    = ncol(clusfoo)
    colnames(clusfoo) = paste0("clusfoo",1:nclusfoo)
    Y           = names(mY_df)
    mY_df       = data.frame(cbind(mY,clusfoo))
    clusfoo     = as.matrix(mY_df%>%group_by(across(all_of(Y)))%>%summarise_at(vars(paste0("clusfoo",1:nclusfoo)),mean))
    clusfoo     = clusfoo[,paste0("clusfoo",1:nclusfoo)]
    mU          = clusfoo
  }
  out             = LCA_fast(mY_unique,freq,iT,mU,reord_user,maxIter,tol,reord)
  P               = ncol(mZ)
  mBeta_init      = matrix(0,P,iT-1)
  mBeta_init[1,]  = out$alphas
  Step1Var        = out$Varmat
  mY              = mY[complete.cases(mZ),]
  mZ              = mZ[complete.cases(mZ),]
  #
  if(verbose){
    if(fixed==1|fixed==2){
      pb=txtProgressBar(char="Fitting structural model...",width=1)
      setTxtProgressBar(pb,1)
      close(pb)
    }
  }
  #
  outcov          = LCAcov(mY,mZ,iT,out$mPhi,mBeta_init,Step1Var,fixed,maxIter,
                           tol,NRtol,NRmaxit)
  return(list(out=out,outcov=outcov,mY=mY,mZ=mZ))
}
#
LCA_fast_init_whigh = function(mY, id_high, iT, kmea = T, maxIter = 1e3, tol = 1e-8, reord = 1, reord_user = 1:iT-1, startval){
  group_by_all = NULL
  mY_df     = data.frame(mY,id_high)
  mY_aggr   = as.matrix(mY_df%>%group_by_all%>%count)
  iHf       = dim(mY_aggr)[2]
  mY_aggr   = mY_aggr[order(mY_aggr[,iHf-1]),]
  freq      = mY_aggr[,iHf]
  mY_unique = mY_aggr[,-c(iHf-1,iHf)]
  if(is.null(startval)&kmea==FALSE){
    clusfoo     = klaR::kmodes(mY_unique,modes=iT,iter.max=100)$cluster
    mU          = vecTomatClass(clusfoo)
  } else if(is.null(startval)&kmea==TRUE){
    prscores    = prcomp(mY_unique)
    num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
    prscores    = prscores$x[,1:num_pc]
    spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
    clusfoo     = spectclust$cluster
    mU          = vecTomatClass(clusfoo)
  } else{
    clusfoo     = vecTomatClass(factor(startval))
    nclusfoo    = ncol(clusfoo)
    colnames(clusfoo) = paste0("clusfoo",1:nclusfoo)
    Y           = names(mY_df)
    mY_df       = data.frame(cbind(mY,id_high,clusfoo))
    clusfoo     = as.matrix(mY_df%>%group_by(across(all_of(c(Y,"id_high"))))%>%summarise_at(vars(paste0("clusfoo",1:nclusfoo)),mean))
    clusfoo     = clusfoo[,paste0("clusfoo",1:nclusfoo)]
    mU          = clusfoo
  }
  out       = LCA_fast(mY_unique,freq,iT,mU,reord_user,maxIter,tol,reord)
  mU        = out$mU[rep(1:nrow(mY_unique),freq),]
  return(list(out=out,mU=mU))
}
#
meas_Init = function(mY, id_high, vNj, iM, iT, kmea = T, maxIter = 1e3, tol = 1e-8, reord = 1, reord_user = 1:iT-1, startval = NULL){
  # Fixed number of low- (iT) and high- (iM) level classes
  iJ = length(vNj)
  iN = dim(mY)[1]
  iK = dim(mY)[2]
  # 
  # Working out starting values at the low level first
  # 
  out_LCA = LCA_fast_init_whigh(mY,id_high,iT,kmea,maxIter,tol,reord,reord_user,startval)
  # 
  # Now turning to the higher level
  # 
  lowlev_relclassprop = matrix(0,iJ,iT)
  foo = 0
  for(j in 1:iJ){
    lowlev_relclassprop[j,] = colMeans(out_LCA$mU[foo+(1:vNj[j]),,drop=FALSE])
    foo = foo +vNj[j]
  }
  high_out      = kmeans(lowlev_relclassprop,centers = iM,iter.max=100,nstart=100)
  vOmega_start  = high_out$size/iJ 
  Wmodal_mat    = vecTomatClass(high_out$cluster)
  # Reordering in decreasing order
  highreord         = order(vOmega_start,decreasing = T)
  Wmodal_mat_reord  = Wmodal_mat[,highreord]
  Wmodal            = apply(Wmodal_mat_reord,1,which.max)
  vOmega_start      = vOmega_start[highreord]
  index_indiv_high  = rep(Wmodal,times=vNj)
  # 
  mPhi_start  = out_LCA$out$mPhi
  mPi_fast    = table(index_indiv_high,rep(out_LCA$out$vModalAssnm+1,times=out_LCA$out$freq))
  mPi_start   = t(mPi_fast/rowSums(mPi_fast))
  if(nrow(mPi_start)!=iT){
    warning(paste0("Initialization suggests less than iT=",iT," lower-level classes."),call.=FALSE)
    mPi_start_replace                                   = matrix(NA,nrow=iT,ncol=iM)
    mPi_start_replace[as.numeric(rownames(mPi_start)),] = mPi_start
    mPi_start_replace[is.na(mPi_start_replace)]         = 1e-4
    mPi_start_replace                                   = apply(mPi_start_replace,2,function(x){x/sum(x)})
    mPi_start                                           = mPi_start_replace
  }
  # 
  return(list(vOmega_start=vOmega_start, mPi_start = mPi_start, mPhi_start = mPhi_start))
}
#
LCA_fast_init_poly = function(mY, iT, ivItemcat, incomplete = F, kmea = T, maxIter = 1e3, tol = 1e-8, reord = 1, reord_user = 1:iT-1, startval = NULL){
  # ivItemcat is the vector of number of categories for each item
  group_by_all  = NULL
  mY_df         = data.frame(mY)
  mY_aggr       = as.matrix(mY_df%>%group_by_all%>%count)
  iHf           = dim(mY_aggr)[2]
  freq          = mY_aggr[,iHf]
  mY_unique     = mY_aggr[,-iHf]
  if(incomplete){
    mDesign = matrix(1,nrow(mY_unique),ncol(mY_unique))
    mDesign[is.na(mY_unique)] = 0
  }
  if(is.null(startval)&kmea==FALSE){
    if(!incomplete){
      clusfoo     = klaR::kmodes(mY_unique,modes=iT,iter.max=100)$cluster
      mU          = vecTomatClass(clusfoo)
    } else{
      clusfoo     = klaR::kmodes(na.omit(mY_unique),modes=iT,iter.max=100)$cluster
      mU          = matrix(NA,nrow(mY_unique),iT)
      mU[complete.cases(mY_unique),] = vecTomatClass(clusfoo)
      mU[!complete.cases(mY_unique),] = 1/iT
      mY_unique[is.na(mY_unique)] = max(na.omit(mY_unique))+1
    }
  } else if(is.null(startval)&kmea==TRUE){
    if(!incomplete){
      prscores    = prcomp(mY_unique)
      num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
      prscores    = prscores$x[,1:num_pc]
      spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
      clusfoo     = spectclust$cluster
      mU          = vecTomatClass(clusfoo)
    } else{
      prscores    = prcomp(na.omit(mY_unique))
      num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
      prscores    = prscores$x[,1:num_pc]
      spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
      clusfoo     = spectclust$cluster
      mU          = matrix(NA,nrow(mY_unique),iT)
      mU[complete.cases(mY_unique),] = vecTomatClass(clusfoo)
      mU[!complete.cases(mY_unique),] = 1/iT
      mY_unique[is.na(mY_unique)] = max(na.omit(mY_unique))+1
    }
  } else{
    clusfoo     = vecTomatClass(factor(startval))
    nclusfoo    = ncol(clusfoo)
    colnames(clusfoo) = paste0("clusfoo",1:nclusfoo)
    Y           = names(mY_df)
    mY_df       = data.frame(cbind(mY,clusfoo))
    clusfoo     = as.matrix(mY_df%>%group_by(across(all_of(Y)))%>%summarise_at(vars(paste0("clusfoo",1:nclusfoo)),mean))
    clusfoo     = clusfoo[,paste0("clusfoo",1:nclusfoo)]
    mU          = clusfoo
  }
  first_poly    = sapply(ivItemcat,function(x){y=x;if(y==2)y=1;return(y)})
  first_poly    = cumsum(first_poly)-first_poly
  if(!incomplete){
    out       = LCA_fast_poly(mY_unique,freq,iT,mU,ivItemcat,first_poly,reord_user,maxIter,tol,reord)
  } else{
    out       = LCA_fast_poly_includeall(mY_unique,mDesign,freq,iT,mU,ivItemcat,first_poly,reord_user,maxIter,tol,reord)
  }
  mY_unique[mY_unique==max(mY_unique)] = NA
  out$mY_unique = mY_unique
  out$freq      = freq
  return(out)
}
#
LCA_fast_init_wcov_poly = function(mY, mZ, iT, ivItemcat, incomplete = F, kmea = T, maxIter = 1e3, tol = 1e-8, fixed = 0, reord = 1, reord_user = 1:iT-1,
                                   NRtol = 1e-6, NRmaxit = 100, verbose, startval){
  # mZ must include a column of ones!
  # ivItemcat is the vector of number of categories for each item
  group_by_all    = NULL
  mY_df           = data.frame(mY)
  mY_aggr         = as.matrix(mY_df%>%group_by_all%>%count)
  iHf             = dim(mY_aggr)[2]
  freq            = mY_aggr[,iHf]
  mY_unique       = mY_aggr[,-iHf]
  if(incomplete){
    mDesign = matrix(1,nrow(mY_unique),ncol(mY_unique))
    mDesign[is.na(mY_unique)] = 0
  }
  if(is.null(startval)&kmea==FALSE){
    if(!incomplete){
      clusfoo     = klaR::kmodes(mY_unique,modes=iT,iter.max=100)$cluster
      mU          = vecTomatClass(clusfoo)
    } else{
      clusfoo     = klaR::kmodes(na.omit(mY_unique),modes=iT,iter.max=100)$cluster
      mU          = matrix(NA,nrow(mY_unique),iT)
      mU[complete.cases(mY_unique),] = vecTomatClass(clusfoo)
      mU[!complete.cases(mY_unique),] = 1/iT
      mY_unique[is.na(mY_unique)] = max(na.omit(mY_unique))+1
    }
  } else if(is.null(startval)&kmea==TRUE){
    if(!incomplete){
      prscores    = prcomp(mY_unique)
      num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
      prscores    = prscores$x[,1:num_pc]
      spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
      clusfoo     = spectclust$cluster
      mU          = vecTomatClass(clusfoo)
    } else{
      prscores    = prcomp(na.omit(mY_unique))
      num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
      prscores    = prscores$x[,1:num_pc]
      spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
      clusfoo     = spectclust$cluster
      mU          = matrix(NA,nrow(mY_unique),iT)
      mU[complete.cases(mY_unique),] = vecTomatClass(clusfoo)
      mU[!complete.cases(mY_unique),] = 1/iT
      mY_unique[is.na(mY_unique)] = max(na.omit(mY_unique))+1
    }
  } else{
    clusfoo     = vecTomatClass(factor(startval))
    nclusfoo    = ncol(clusfoo)
    colnames(clusfoo) = paste0("clusfoo",1:nclusfoo)
    Y           = names(mY_df)
    mY_df       = data.frame(cbind(mY,clusfoo))
    clusfoo     = as.matrix(mY_df%>%group_by(across(all_of(Y)))%>%summarise_at(vars(paste0("clusfoo",1:nclusfoo)),mean))
    clusfoo     = clusfoo[,paste0("clusfoo",1:nclusfoo)]
    mU          = clusfoo
  }
  first_poly      = sapply(ivItemcat,function(x){y=x;if(y==2)y=1;return(y)})
  first_poly      = cumsum(first_poly)-first_poly
  if(!incomplete){
    out       = LCA_fast_poly(mY_unique,freq,iT,mU,ivItemcat,first_poly,reord_user,maxIter,tol,reord)
  } else{
    out       = LCA_fast_poly_includeall(mY_unique,mDesign,freq,iT,mU,ivItemcat,first_poly,reord_user,maxIter,tol,reord)
  }
  P               = ncol(mZ)
  mBeta_init      = matrix(0,P,iT-1)
  mBeta_init[1,]  = out$alphas
  Step1Var        = out$Varmat
  mY              = mY[complete.cases(mZ),]
  mZ              = mZ[complete.cases(mZ),]
  #
  if(verbose){
    if(fixed==1|fixed==2){
      pb=txtProgressBar(char="Fitting structural model...",width=1)
      setTxtProgressBar(pb,1)
      close(pb)
    }
  }
  #
  if(!incomplete){
    outcov          = LCAcov_poly(mY,mZ,iT,out$mPhi,mBeta_init,Step1Var,ivItemcat,fixed,maxIter,
                                  tol,NRtol,NRmaxit)
  } else{
    mDesign = matrix(1,nrow(mY),ncol(mY))
    mDesign[is.na(mY)] = 0
    mY[is.na(mY)] = max(na.omit(mY))+1
    outcov          = LCAcov_poly_includeall(mY,mDesign,mZ,iT,out$mPhi,mBeta_init,Step1Var,ivItemcat,fixed,maxIter,
                                             tol,NRtol,NRmaxit)
  }
  mY[mY==max(mY)] = NA
  return(list(out=out,outcov=outcov,mY=mY,mZ=mZ))
}
#
LCA_fast_init_whigh_poly = function(mY, id_high, iT, ivItemcat, incomplete = F, kmea = T, maxIter = 1e3, tol = 1e-8, reord = 1, reord_user = 1:iT-1, startval){
  # ivItemcat is the vector of number of categories for each item
  group_by_all = NULL
  mY_df     = data.frame(mY,id_high)
  mY_aggr   = as.matrix(mY_df%>%group_by_all%>%count)
  iHf       = dim(mY_aggr)[2]
  mY_aggr   = mY_aggr[order(mY_aggr[,iHf-1]),]
  freq      = mY_aggr[,iHf]
  mY_unique = mY_aggr[,-c(iHf-1,iHf)]
  if(incomplete){
    mDesign = matrix(1,nrow(mY_unique),ncol(mY_unique))
    mDesign[is.na(mY_unique)] = 0
  }
  if(is.null(startval)&kmea==FALSE){
    if(!incomplete){
      clusfoo     = klaR::kmodes(mY_unique,modes=iT,iter.max=100)$cluster
      mU          = vecTomatClass(clusfoo)
    } else{
      clusfoo     = klaR::kmodes(na.omit(mY_unique),modes=iT,iter.max=100)$cluster
      mU          = matrix(NA,nrow(mY_unique),iT)
      mU[complete.cases(mY_unique),] = vecTomatClass(clusfoo)
      mU[!complete.cases(mY_unique),] = 1/iT
      mY_unique[is.na(mY_unique)] = max(na.omit(mY_unique))+1
    }
  } else if(is.null(startval)&kmea==TRUE){
    if(!incomplete){
      prscores    = prcomp(mY_unique)
      num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
      prscores    = prscores$x[,1:num_pc]
      spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
      clusfoo     = spectclust$cluster
      mU          = vecTomatClass(clusfoo)
    } else{
      prscores    = prcomp(na.omit(mY_unique))
      num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
      prscores    = prscores$x[,1:num_pc]
      spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
      clusfoo     = spectclust$cluster
      mU          = matrix(NA,nrow(mY_unique),iT)
      mU[complete.cases(mY_unique),] = vecTomatClass(clusfoo)
      mU[!complete.cases(mY_unique),] = 1/iT
      mY_unique[is.na(mY_unique)] = max(na.omit(mY_unique))+1
    }
  } else{
    clusfoo     = vecTomatClass(factor(startval))
    nclusfoo    = ncol(clusfoo)
    colnames(clusfoo) = paste0("clusfoo",1:nclusfoo)
    Y           = names(mY_df)
    mY_df       = data.frame(cbind(mY,id_high,clusfoo))
    clusfoo     = as.matrix(mY_df%>%group_by(across(all_of(c(Y,"id_high"))))%>%summarise_at(vars(paste0("clusfoo",1:nclusfoo)),mean))
    clusfoo     = clusfoo[,paste0("clusfoo",1:nclusfoo)]
    mU          = clusfoo
  }
  first_poly  = sapply(ivItemcat,function(x){y=x;if(y==2)y=1;return(y)})
  first_poly  = cumsum(first_poly)-first_poly
  if(!incomplete){
    out       = LCA_fast_poly(mY_unique,freq,iT,mU,ivItemcat,first_poly,reord_user,maxIter,tol,reord)
  } else{
    out       = LCA_fast_poly_includeall(mY_unique,mDesign,freq,iT,mU,ivItemcat,first_poly,reord_user,maxIter,tol,reord)
  }
  mU          = out$mU[rep(1:nrow(mY_unique),freq),]
  return(list(out=out,mU=mU,first_poly=first_poly))
}
#
meas_Init_poly = function(mY, id_high, vNj, iM, iT, ivItemcat, incomplete = F, kmea = T, maxIter = 1e3, tol = 1e-8, reord = 1, reord_user = 1:iT-1, startval = NULL){
  # Fixed number of low- (iT) and high- (iM) level classes
  # ivItemcat is the vector of number of categories for each item
  iJ = length(vNj)
  iN = dim(mY)[1]
  iK = dim(mY)[2]
  # 
  # Working out starting values at the low level first
  # 
  out_LCA = LCA_fast_init_whigh_poly(mY,id_high,iT,ivItemcat,incomplete,kmea,maxIter,tol,reord,reord_user,startval)
  # 
  # now turning to higher level
  # 
  lowlev_relclassprop = matrix(0,iJ,iT)
  foo = 0
  for(j in 1:iJ){
    lowlev_relclassprop[j,] = colMeans(out_LCA$mU[foo+(1:vNj[j]),,drop=FALSE])
    foo = foo +vNj[j]
  }
  high_out      = kmeans(lowlev_relclassprop,centers = iM,iter.max=100,nstart=100)
  vOmega_start  = high_out$size/iJ 
  Wmodal_mat    = vecTomatClass(high_out$cluster)
  # Reordering in decreasing order
  highreord         = order(vOmega_start,decreasing = T)
  Wmodal_mat_reord  = Wmodal_mat[,highreord]
  Wmodal            = apply(Wmodal_mat_reord,1,which.max)
  vOmega_start      = vOmega_start[highreord]
  index_indiv_high  = rep(Wmodal,times=vNj)
  # 
  mPhi_start  = out_LCA$out$mPhi
  mPi_fast    = table(index_indiv_high,rep(out_LCA$out$vModalAssnm+1,times=out_LCA$out$freq))
  mPi_start   = t(mPi_fast/rowSums(mPi_fast))
  if(nrow(mPi_start)!=iT){
    warning(paste0("Initialization suggests less than iT=",iT," lower-level classes."),call.=FALSE)
    mPi_start_replace                                   = matrix(NA,nrow=iT,ncol=iM)
    mPi_start_replace[as.numeric(rownames(mPi_start)),] = mPi_start
    mPi_start_replace[is.na(mPi_start_replace)]         = 1e-4
    mPi_start_replace                                   = apply(mPi_start_replace,2,function(x){x/sum(x)})
    mPi_start                                           = mPi_start_replace
  }
  # 
  return(list(vOmega_start=vOmega_start, mPi_start = mPi_start, mPhi_start = mPhi_start, first_poly = out_LCA$first_poly))
}
#
#
#
simultsel_fun = function(mY,id_high,iT,iM,kmea,maxIter,tol,reord){
  LCAout  = NULL
  iH      = ncol(mY)
  K       = iH
  vNj     = table(id_high-1)
  iN      = length(vNj)
  if(iM ==1 & iT ==1){
    ll          = sum(apply(mY,2,function(x){dbinom(x,1,mean(x),log=T)}))
    BIClow      = -2*ll + K*log(sum(vNj))
    BIChigh     = -2*ll + K*log(iN)
    AIC         = -2*ll + 2*K
    ICL_BIClow  = BIClow
    ICL_BIChigh = BIChigh
    outmuLCA    = list()
  }else if(iM > 1 & iT == 1){
    ll          = -Inf
    BIClow      = Inf
    BIChigh     = Inf
    AIC         = Inf
    ICL_BIClow  = Inf
    ICL_BIChigh = Inf
    outmuLCA    = list()
  }else if(iM==1 & iT > 1){
    LCAout      = LCA_fast_init(mY,iT,kmea,maxIter,tol,reord)
    ll          = tail(LCAout$LLKSeries,1)
    npar        = iT*iH + iT - 1
    BIClow      = LCAout$BIC
    BIChigh     = -2*ll + npar*log(iN)
    AIC         = -2*ll + 2*npar
    ICL_BIClow  = BIClow + 2*sum(LCAout$freq*apply(-LCAout$mU*log(LCAout$mU),1,sum))
    ICL_BIChigh = Inf
    outmuLCA    = LCAout
  }else{
    # Actual multilevel LC fit
    start       = meas_Init(mY,id_high,vNj,iM,iT,kmea,maxIter,tol,reord)
    vOmegast    = start$vOmega_start
    mPi         = start$mPi_start
    outmuLCA    = MLTLCA(mY,vNj,vOmegast,mPi,start$mPhi_start,1:iT-1,1:iM-1,maxIter,tol,reord)
    ll          = tail(outmuLCA$LLKSeries,1) 
    BIClow      = outmuLCA$BIClow
    BIChigh     = outmuLCA$BIChigh
    AIC         = outmuLCA$AIC
    ICL_BIClow  = outmuLCA$ICL_BIClow
    ICL_BIChigh = outmuLCA$ICL_BIChigh
  }
  return(list(modFIT = outmuLCA,ll=ll, BIClow = BIClow, BIChigh = BIChigh, AIC = AIC,
              ICL_BIClow = ICL_BIClow, ICL_BIChigh = ICL_BIChigh))
}
#
lukosel_fun = function(mY,id_high,iT_range,iM_range,verbose,kmea,maxIter,tol,reord){
  iH      = ncol(mY)
  K       = iH
  vNj     = table(id_high-1)
  iN      = length(vNj)
  iT_max  = max(iT_range)
  iT_min  = min(iT_range)
  num_iT  = iT_max-iT_min+1
  iM_max  = max(iM_range)
  iM_min  = min(iM_range)
  num_iM  = iM_max-iM_min+1
  # Step 1 - lower level
  LCAout            = list()
  BIClow_step1      = rep(NA,num_iT)
  BIChigh_step1     = rep(NA,num_iT)
  AIC_step1         = rep(NA,num_iT)
  ICL_BIClow_step1  = rep(NA,num_iT)
  ICL_BIChigh_step1 = rep(NA,num_iT)
  for(i in iT_range){
    #
    if(verbose){
      pb=txtProgressBar(char="Fitting single-level measurement model (step 1)...",width=1)
      setTxtProgressBar(pb,1)
      close(pb)
    }
    #
    if(i==1){
      llfoo                 = sum(apply(mY,2,function(x){dbinom(x,1,mean(x),log=T)}))
      BIClow_step1[1]       = -2*llfoo + K*log(sum(vNj))
      BIChigh_step1[1]      = -2*llfoo + K*log(iN)
      AIC_step1[1]          = -2*llfoo + 2*K
      ICL_BIClow_step1[1]   = BIClow_step1[1] 
      ICL_BIChigh_step1[1]  = BIChigh_step1[1]
    } else{
      LCAout                        = LCA_fast_init(mY,i,kmea,maxIter,tol,reord)
      ll                            = tail(LCAout$LLKSeries,1)
      npar                          = i*iH + i - 1
      BIClow_step1[1+i-iT_min]      = LCAout$BIC
      BIChigh_step1[1+i-iT_min]     = -2*ll + npar*log(iN)
      AIC_step1[1+i-iT_min]         = -2*ll + 2*npar
      ICL_BIClow_step1[1+i-iT_min]  = BIClow_step1[i] + 2*sum(LCAout$freq*apply(-LCAout$mU*log(LCAout$mU),1,sum))
      ICL_BIChigh_step1[1+i-iT_min] = Inf
    }
  }
  iT_currbest = which.min(BIClow_step1)+iT_min-1
  if(iT_currbest==1){
    out = list(BIClow_step1 = BIClow_step1, BIChigh_step1 = BIChigh_step1, AIC_step1 = AIC_step1,
               ICL_BIClow_step1 = ICL_BIClow_step1, ICL_BIChigh_step1 = ICL_BIChigh_step1)
    step1mat = matrix(as.character(round(unlist(out),2)),length(iT_range),5,
                      dimnames=list(paste0("iT=",iT_range),
                                    c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
    step1mat[step1mat=="Inf"|step1mat=="-Inf"|is.na(step1mat)] = "-"
    if(verbose)print(noquote(step1mat))
    stop("iT=1 optimal.",call.=FALSE)
  }
  # Step 2 - higher level
  outmuLCA2         = list()
  BIClow_step2      = rep(NA,num_iM)
  BIChigh_step2     = rep(NA,num_iM)
  AIC_step2         = rep(NA,num_iM)
  ICL_BIClow_step2  = rep(NA,num_iM)
  ICL_BIChigh_step2 = rep(NA,num_iM)
  if(iT_currbest > 1){
    for(i in iM_range){
      #
      if(verbose){
        pb=txtProgressBar(char="Fitting multilevel measurement model (step 2)...",width=1)
        setTxtProgressBar(pb,1)
        close(pb)
      }
      #
      if(i==1){
        BIClow_step2[1]       = BIClow_step1[1+iT_currbest-iT_min]
        BIChigh_step2[1]      = BIChigh_step1[1+iT_currbest-iT_min]
        AIC_step2[1]          = AIC_step1[1+iT_currbest-iT_min]
        ICL_BIClow_step2[1]   = ICL_BIClow_step1[1+iT_currbest-iT_min]
        ICL_BIChigh_step2[1]  = ICL_BIChigh_step1[1+iT_currbest-iT_min]
      } else{
        start     = meas_Init(mY,id_high,vNj,i,iT_currbest,kmea,maxIter,tol,reord)
        vOmegast  = start$vOmega_start
        mPi       = start$mPi_start
        outmuLCA2[[1+i-iM_min]]         = MLTLCA(mY,vNj,vOmegast,mPi,start$mPhi_start,1:iT_currbest-1,1:i-1,maxIter,tol,reord)
        BIClow_step2[1+i-iM_min]        = outmuLCA2[[1+i-iM_min]]$BIClow
        BIChigh_step2[1+i-iM_min]       = outmuLCA2[[1+i-iM_min]]$BIChigh
        AIC_step2[1+i-iM_min]           = outmuLCA2[[1+i-iM_min]]$AIC
        ICL_BIClow_step2[1+i-iM_min]    = outmuLCA2[[1+i-iM_min]]$ICL_BIClow
        ICL_BIChigh_step2[1+i-iM_min]   = outmuLCA2[[1+i-iM_min]]$ICL_BIChigh
      }
    }
    iM_currbest = which.min(BIChigh_step2)+iM_min-1
  } else{
    iM_currbest = 1
  }
  # Step 3 - revisiting lower level
  outmuLCA3         = list()
  BIClow_step3      = rep(NA,num_iT)
  BIChigh_step3     = rep(NA,num_iT)
  AIC_step3         = rep(NA,num_iT)
  ICL_BIClow_step3  = rep(NA,num_iT)
  ICL_BIChigh_step3 = rep(NA,num_iT)
  if(iM_currbest > 1){
    for(i in iT_range){
      #
      if(verbose){
        pb=txtProgressBar(char="Fitting multilevel measurement model (step 3)...",width=1)
        setTxtProgressBar(pb,1)
        close(pb)
      }
      #
      if(i>1){
        start     = meas_Init(mY,id_high,vNj,iM_currbest,i,kmea,maxIter,tol,reord)
        vOmegast  = start$vOmega_start
        mPi       = start$mPi_start
        outmuLCA3[[1+i-iT_min]]       = MLTLCA(mY,vNj,vOmegast,mPi,start$mPhi_start,1:i-1,1:iM_currbest-1,maxIter,tol,reord)
        BIClow_step3[1+i-iT_min]      = outmuLCA3[[1+i-iT_min]]$BIClow
        BIChigh_step3[1+i-iT_min]     = outmuLCA3[[1+i-iT_min]]$BIChigh
        AIC_step3[1+i-iT_min]         = outmuLCA3[[1+i-iT_min]]$AIC
        ICL_BIClow_step3[1+i-iT_min]  = outmuLCA3[[1+i-iT_min]]$ICL_BIClow
        ICL_BIChigh_step3[1+i-iT_min] = outmuLCA3[[1+i-iT_min]]$ICL_BIChigh
      }
    }
    iT_currbest     = which(BIClow_step3==min(BIClow_step3,na.rm=T))+iT_min-1
    outmuLCA_step3  = outmuLCA3[[1+iT_currbest-iT_min]]
  } else{
    outmuLCA_step3 = LCA_fast_init(mY,iT_currbest,kmea,maxIter,tol,reord)
  }
  return(list(outmuLCA_step3=outmuLCA_step3,iT_opt=iT_currbest, iM_opt=iM_currbest, 
              BIClow_step1 = BIClow_step1, BIChigh_step1 = BIChigh_step1, AIC_step1 = AIC_step1,
              ICL_BIClow_step1 = ICL_BIClow_step1, ICL_BIChigh_step1 = ICL_BIChigh_step1,
              BIClow_step2 = BIClow_step2, BIChigh_step2 = BIChigh_step2, AIC_step2 = AIC_step2,
              ICL_BIClow_step2 = ICL_BIClow_step2, ICL_BIChigh_step2 = ICL_BIChigh_step2,
              BIClow_step3 = BIClow_step3, BIChigh_step3 = BIChigh_step3, AIC_step3 = AIC_step3,
              ICL_BIClow_step3 = ICL_BIClow_step3, ICL_BIChigh_step3 = ICL_BIChigh_step3))
}
#
sel_other = function(mY,id_high,iT_range,iM_range,approach,verbose,kmea,maxIter,tol,reord){
  iH      = ncol(mY)
  K       = iH
  vNj     = table(id_high-1)
  iN      = length(vNj)
  iT_max  = max(iT_range)
  iT_min  = min(iT_range)
  num_iT  = iT_max-iT_min+1
  if(!is.null(iM_range)){
    iM_max  = max(iM_range)
    iM_min  = min(iM_range)
    num_iM  = iM_max-iM_min+1
    if(num_iM==1){
      if(iM_range==1){
        iM_range=NULL
      }
    }
  }
  #
  if(approach=="model selection on low"|approach=="model selection on low with high"){
    BIClow      = rep(NA,num_iT)
    BIChigh     = rep(NA,num_iT)
    AIC         = rep(NA,num_iT)
    ICL_BIClow  = rep(NA,num_iT)
    ICL_BIChigh = rep(NA,num_iT)
    if(is.null(iM_range)){
      outmuLCA = list()
      for(i in iT_range){
        #
        if(verbose){
          pb=txtProgressBar(char="Fitting measurement model...",width=1)
          setTxtProgressBar(pb,1)
          close(pb)
        }
        #
        if(i==1){
          outmuLCA[[1]]   = NA
          llfoo           = sum(apply(mY,2,function(x){dbinom(x,1,mean(x),log=T)}))
          BIClow[1]       = -2*llfoo + K*log(sum(vNj))
          if(is.null(iM_range)) BIClow[1] = -2*llfoo + K*log(nrow(mY))
          if(!is.null(iM_range)) BIChigh[1] = -2*llfoo + K*log(iN)
          AIC[1]          = -2*llfoo + 2*K
          ICL_BIClow[1]   = BIClow[1]
          ICL_BIChigh[1]  = BIChigh[1]
        } else{
          outmuLCA[[1+i-iT_min]]  = LCA_fast_init(mY,i,kmea,maxIter,tol,reord)
          ll                      = tail(outmuLCA[[1+i-iT_min]]$LLKSeries,1)
          npar                    = i*iH + i - 1
          BIClow[1+i-iT_min]      = outmuLCA[[1+i-iT_min]]$BIC
          BIChigh[1+i-iT_min]     = -2*ll + npar*log(iN)
          AIC[1+i-iT_min]         = -2*ll + 2*npar
          ICL_BIClow[1+i-iT_min]  = BIClow[i] + 2*sum(outmuLCA[[1+i-iT_min]]$freq*apply(-outmuLCA[[1+i-iT_min]]$mU*log(outmuLCA[[1+i-iT_min]]$mU),1,sum))
          ICL_BIChigh[1+i-iT_min] = Inf
        }
      }
      iT_best   = which(BIClow==min(BIClow,na.rm=T))+iT_min-1
      outmuLCA  = outmuLCA[[1+iT_best-iT_min]]
    } else{
      outmuLCA  = list()
      for(i in iT_range){
        #
        if(verbose){
          pb=txtProgressBar(char="Fitting measurement model...",width=1)
          setTxtProgressBar(pb,1)
          close(pb)
        }
        #
        if(i>1){
          start     = meas_Init(mY,id_high,vNj,iM_range,i,kmea,maxIter,tol,reord)
          vOmegast  = start$vOmega_start
          mPi       = start$mPi_start
          outmuLCA[[1+i-iT_min]]  = MLTLCA(mY,vNj,vOmegast,mPi,start$mPhi_start,1:i-1,1:iM_range-1,maxIter,tol,reord)
          BIClow[1+i-iT_min]      = outmuLCA[[1+i-iT_min]]$BIClow
          BIChigh[1+i-iT_min]     = outmuLCA[[1+i-iT_min]]$BIChigh
          AIC[1+i-iT_min]         = outmuLCA[[1+i-iT_min]]$AIC
          ICL_BIClow[1+i-iT_min]  = outmuLCA[[1+i-iT_min]]$ICL_BIClow
          ICL_BIChigh[1+i-iT_min] = outmuLCA[[1+i-iT_min]]$ICL_BIChigh
        }
      }
      iT_best   = which(BIClow==min(BIClow,na.rm=T))+iT_min-1
      outmuLCA  = outmuLCA[[1+iT_best-iT_min]]
    }
    if(iT_best==1) stop("iT=1 optimal.",call.=FALSE)
    return(list(outmuLCA=outmuLCA,iT_best=iT_best,
                BIClow=BIClow,BIChigh=BIChigh,AIC=AIC,
                ICL_BIClow=ICL_BIClow,ICL_BIChigh=ICL_BIChigh))
  } else if(approach=="model selection on high"){
    BIClow      = rep(NA,num_iM)
    BIChigh     = rep(NA,num_iM)
    AIC         = rep(NA,num_iM)
    ICL_BIClow  = rep(NA,num_iM)
    ICL_BIChigh = rep(NA,num_iM)
    if(iT_range > 1){
      outmuLCA = list()
      for(i in iM_range){
        #
        if(verbose){
          pb=txtProgressBar(char="Fitting measurement model...",width=1)
          setTxtProgressBar(pb,1)
          close(pb)
        }
        #
        if(i==1){
          outmuLCA[[1]]   = LCA_fast_init(mY,iT_range,kmea,maxIter,tol,reord)
          ll              = tail(outmuLCA[[1]]$LLKSeries,1)
          npar            = iT_range*iH + iT_range - 1
          BIClow[1]       = outmuLCA[[1]]$BIC
          BIChigh[1]      = -2*ll + npar*log(iN)
          AIC[1]          = -2*ll + 2*npar
          ICL_BIClow[1]   = BIClow[1] + 2*sum(outmuLCA[[1]]$freq*apply(-outmuLCA[[1]]$mU*log(outmuLCA[[1]]$mU),1,sum))
          ICL_BIChigh[1]  = Inf
        } else{
          start     = meas_Init(mY,id_high,vNj,i,iT_range,kmea,maxIter,tol,reord)
          vOmegast  = start$vOmega_start
          mPi       = start$mPi_start
          outmuLCA[[1+i-iM_min]]    = MLTLCA(mY,vNj,vOmegast,mPi,start$mPhi_start,1:iT_range-1,1:i-1,maxIter,tol,reord)
          BIClow[1+i-iM_min]        = outmuLCA[[1+i-iM_min]]$BIClow
          BIChigh[1+i-iM_min]       = outmuLCA[[1+i-iM_min]]$BIChigh
          AIC[1+i-iM_min]           = outmuLCA[[1+i-iM_min]]$AIC
          ICL_BIClow[1+i-iM_min]    = outmuLCA[[1+i-iM_min]]$ICL_BIClow
          ICL_BIChigh[1+i-iM_min]   = outmuLCA[[1+i-iM_min]]$ICL_BIChigh
        }
      }
      iM_best   = which.min(BIChigh)+iM_min-1
      outmuLCA  = outmuLCA[[1+iM_best-iM_min]]
    } else{
      stop("iM=1 optimal.",call.=FALSE)
    }
    return(list(outmuLCA=outmuLCA,iM_best=iM_best,
                BIClow=BIClow,BIChigh=BIChigh,AIC=AIC,
                ICL_BIClow=ICL_BIClow,ICL_BIChigh=ICL_BIChigh))
  }
}
#
simultsel_fun_poly = function(mY,mDesign,id_high,iT,iM,ivItemcat,kmea,maxIter,tol,reord){
  LCAout  = NULL
  iH      = sum(ivItemcat-1)
  K       = iH
  vNj     = table(id_high-1)
  iN      = length(vNj)
  incomplete  = !is.null(mDesign)
  if(iM ==1 & iT ==1){
    ll          = sum(apply(mY,2,function(x){dbinom(x,1,mean(x,na.rm=T),log=T)}),na.rm=T)
    BIClow      = -2*ll + K*log(sum(vNj))
    BIChigh     = -2*ll + K*log(iN)
    AIC         = -2*ll + 2*K
    ICL_BIClow  = BIClow
    ICL_BIChigh = BIChigh
    outmuLCA    = list()
  }else if(iM > 1 & iT == 1){
    ll          = -Inf
    BIClow      = Inf
    BIChigh     = Inf
    AIC         = Inf
    ICL_BIClow  = Inf
    ICL_BIChigh = Inf
    outmuLCA    = list()
  }else if(iM==1 & iT > 1){
    LCAout      = LCA_fast_init_poly(mY,iT,ivItemcat,incomplete,kmea,maxIter,tol,reord)
    ll          = tail(LCAout$LLKSeries,1)
    npar        = iT*iH + iT - 1
    BIClow      = LCAout$BIC
    BIChigh     = -2*ll + npar*log(iN)
    AIC         = -2*ll + 2*npar
    ICL_BIClow  = BIClow + 2*sum(LCAout$freq*apply(-LCAout$mU*log(LCAout$mU),1,sum))
    ICL_BIChigh = Inf
    outmuLCA    = LCAout
  }else{
    # Actual multilevel LC fit
    start       = meas_Init_poly(mY,id_high,vNj,iM,iT,ivItemcat,incomplete,kmea,maxIter,tol,reord)
    vOmegast    = start$vOmega_start
    mPi         = start$mPi_start
    if(!incomplete){
      outmuLCA    = MLTLCA_poly(mY,vNj,vOmegast,mPi,
                                start$mPhi_start,ivItemcat,start$first_poly,1:iT-1,1:iM-1,maxIter,tol,reord)
    } else{
      mY[is.na(mY)] = max(na.omit(mY))+1
      outmuLCA    = MLTLCA_poly_includeall(mY,mDesign,vNj,vOmegast,mPi,
                                           start$mPhi_start,ivItemcat,start$first_poly,1:iT-1,1:iM-1,maxIter,tol,reord)
    }
    ll          = tail(outmuLCA$LLKSeries,1)
    BIClow      = outmuLCA$BIClow
    BIChigh     = outmuLCA$BIChigh
    AIC         = outmuLCA$AIC
    ICL_BIClow  = outmuLCA$ICL_BIClow
    ICL_BIChigh = outmuLCA$ICL_BIChigh
  }
  return(list(modFIT = outmuLCA,ll=ll, BIClow = BIClow, BIChigh = BIChigh, AIC = AIC,
              ICL_BIClow = ICL_BIClow, ICL_BIChigh = ICL_BIChigh))
}
#
lukosel_fun_poly = function(mY,mDesign,id_high,iT_range,iM_range,ivItemcat,verbose,kmea,maxIter,tol,reord){
  iH      = sum(ivItemcat-1)
  K       = iH
  vNj     = table(id_high-1)
  iN      = length(vNj)
  iT_max  = max(iT_range)
  iT_min  = min(iT_range)
  num_iT  = iT_max-iT_min+1
  iM_max  = max(iM_range)
  iM_min  = min(iM_range)
  num_iM  = iM_max-iM_min+1
  incomplete  = !is.null(mDesign)
  # Step 1 - lower level
  LCAout            = list()
  BIClow_step1      = rep(NA,num_iT)
  BIChigh_step1     = rep(NA,num_iT)
  AIC_step1         = rep(NA,num_iT)
  ICL_BIClow_step1  = rep(NA,num_iT)
  ICL_BIChigh_step1 = rep(NA,num_iT)
  for(i in iT_range){
    #
    if(verbose){
      pb=txtProgressBar(char="Fitting single-level measurement model (step 1)...",width=1)
      setTxtProgressBar(pb,1)
      close(pb)
    }
    #
    if(i==1){
      llfoo                 = sum(apply(mY,2,function(x){dbinom(x,1,mean(x,na.rm=T),log=T)}),na.rm=T)
      BIClow_step1[1]       = -2*llfoo + K*log(sum(vNj))
      BIChigh_step1[1]      = -2*llfoo + K*log(iN)
      AIC_step1[1]          = -2*llfoo + 2*K
      ICL_BIClow_step1[1]   = BIClow_step1[1]
      ICL_BIChigh_step1[1]  = BIChigh_step1[1]
    } else{
      LCAout                        = LCA_fast_init_poly(mY,i,ivItemcat,incomplete,kmea,maxIter,tol,reord)
      ll                            = tail(LCAout$LLKSeries,1)
      npar                          = i*iH + i - 1
      BIClow_step1[1+i-iT_min]      = LCAout$BIC
      BIChigh_step1[1+i-iT_min]     = -2*ll + npar*log(iN)
      AIC_step1[1+i-iT_min]         = -2*ll + 2*npar
      ICL_BIClow_step1[1+i-iT_min]  = BIClow_step1[i] + 2*sum(LCAout$freq*apply(-LCAout$mU*log(LCAout$mU),1,sum))
      ICL_BIChigh_step1[1+i-iT_min] = Inf
    }
  }
  iT_currbest = which.min(BIClow_step1)+iT_min-1
  if(iT_currbest==1){
    out = list(BIClow_step1 = BIClow_step1, BIChigh_step1 = BIChigh_step1, AIC_step1 = AIC_step1,
               ICL_BIClow_step1 = ICL_BIClow_step1, ICL_BIChigh_step1 = ICL_BIChigh_step1)
    step1mat = matrix(as.character(round(unlist(out),2)),length(iT_range),5,
                      dimnames=list(paste0("iT=",iT_range),
                                    c("BIClow","BIChigh","AIC","ICL_BIClow","ICL_BIChigh")))
    step1mat[step1mat=="Inf"|step1mat=="-Inf"|is.na(step1mat)] = "-"
    if(verbose)print(noquote(step1mat))
    stop("iT=1 optimal.",call.=FALSE)
  }
  # Step 2 - higher level
  outmuLCA2         = list()
  BIClow_step2      = rep(NA,num_iM)
  BIChigh_step2     = rep(NA,num_iM)
  AIC_step2         = rep(NA,num_iM)
  ICL_BIClow_step2  = rep(NA,num_iM)
  ICL_BIChigh_step2 = rep(NA,num_iM)
  if(iT_currbest > 1){
    for(i in iM_range){
      #
      if(verbose){
        pb=txtProgressBar(char="Fitting multilevel measurement model (step 2)...",width=1)
        setTxtProgressBar(pb,1)
        close(pb)
      }
      #
      if(i==1){
        BIClow_step2[1]       = BIClow_step1[1+iT_currbest-iT_min]
        BIChigh_step2[1]      = BIChigh_step1[1+iT_currbest-iT_min]
        AIC_step2[1]          = AIC_step1[1+iT_currbest-iT_min]
        ICL_BIClow_step2[1]   = ICL_BIClow_step1[1+iT_currbest-iT_min]
        ICL_BIChigh_step2[1]  = ICL_BIChigh_step1[1+iT_currbest-iT_min]
      } else{
        start     = meas_Init_poly(mY,id_high,vNj,i,iT_currbest,ivItemcat,incomplete,kmea,maxIter,tol,reord)
        vOmegast  = start$vOmega_start
        mPi       = start$mPi_start
        if(!incomplete){
          outmuLCA2[[1+i-iM_min]]         = MLTLCA_poly(mY,vNj,vOmegast,mPi,start$mPhi_start,ivItemcat,start$first_poly,1:iT_currbest-1,1:i-1,maxIter,tol,reord)
        } else{
          mY[is.na(mY)] = max(na.omit(mY))+1
          outmuLCA2[[1+i-iM_min]]         = MLTLCA_poly_includeall(mY,mDesign,vNj,vOmegast,mPi,start$mPhi_start,ivItemcat,start$first_poly,1:iT_currbest-1,1:i-1,maxIter,tol,reord)
        }
        BIClow_step2[1+i-iM_min]        = outmuLCA2[[1+i-iM_min]]$BIClow
        BIChigh_step2[1+i-iM_min]       = outmuLCA2[[1+i-iM_min]]$BIChigh
        AIC_step2[1+i-iM_min]           = outmuLCA2[[1+i-iM_min]]$AIC
        ICL_BIClow_step2[1+i-iM_min]    = outmuLCA2[[1+i-iM_min]]$ICL_BIClow
        ICL_BIChigh_step2[1+i-iM_min]   = outmuLCA2[[1+i-iM_min]]$ICL_BIChigh
      }
    }
    iM_currbest = which.min(BIChigh_step2)+iM_min-1
  } else{
    iM_currbest = 1
  }
  # Step 3 - revisiting lower level
  outmuLCA3         = list()
  BIClow_step3      = rep(NA,num_iT)
  BIChigh_step3     = rep(NA,num_iT)
  AIC_step3         = rep(NA,num_iT)
  ICL_BIClow_step3  = rep(NA,num_iT)
  ICL_BIChigh_step3 = rep(NA,num_iT)
  if(iM_currbest > 1){
    for(i in iT_range){
      #
      if(verbose){
        pb=txtProgressBar(char="Fitting multilevel measurement model (step 3)...",width=1)
        setTxtProgressBar(pb,1)
        close(pb)
      }
      #
      if(i>1){
        start     = meas_Init_poly(mY,id_high,vNj,iM_currbest,i,ivItemcat,incomplete,kmea,maxIter,tol,reord)
        vOmegast  = start$vOmega_start
        mPi       = start$mPi_start
        if(!incomplete){
          outmuLCA3[[1+i-iT_min]]       = MLTLCA_poly(mY,vNj,vOmegast,mPi,start$mPhi_start,ivItemcat,start$first_poly,1:i-1,1:iM_currbest-1,maxIter,tol,reord)
        } else{
          mY[is.na(mY)] = max(na.omit(mY))+1
          outmuLCA3[[1+i-iT_min]]       = MLTLCA_poly_includeall(mY,mDesign,vNj,vOmegast,mPi,start$mPhi_start,ivItemcat,start$first_poly,1:i-1,1:iM_currbest-1,maxIter,tol,reord)
        }
        BIClow_step3[1+i-iT_min]      = outmuLCA3[[1+i-iT_min]]$BIClow
        BIChigh_step3[1+i-iT_min]     = outmuLCA3[[1+i-iT_min]]$BIChigh
        AIC_step3[1+i-iT_min]         = outmuLCA3[[1+i-iT_min]]$AIC
        ICL_BIClow_step3[1+i-iT_min]  = outmuLCA3[[1+i-iT_min]]$ICL_BIClow
        ICL_BIChigh_step3[1+i-iT_min] = outmuLCA3[[1+i-iT_min]]$ICL_BIChigh
      }
    }
    iT_currbest     = which(BIClow_step3==min(BIClow_step3,na.rm=T))+iT_min-1
    outmuLCA_step3  = outmuLCA3[[1+iT_currbest-iT_min]]
  } else{
    outmuLCA_step3 = LCA_fast_init_poly(mY,iT_currbest,ivItemcat,incomplete,kmea,maxIter,tol,reord)
  }
  return(list(outmuLCA_step3=outmuLCA_step3,iT_opt=iT_currbest, iM_opt=iM_currbest, 
              BIClow_step1 = BIClow_step1, BIChigh_step1 = BIChigh_step1, AIC_step1 = AIC_step1,
              ICL_BIClow_step1 = ICL_BIClow_step1, ICL_BIChigh_step1 = ICL_BIChigh_step1,
              BIClow_step2 = BIClow_step2, BIChigh_step2 = BIChigh_step2, AIC_step2 = AIC_step2,
              ICL_BIClow_step2 = ICL_BIClow_step2, ICL_BIChigh_step2 = ICL_BIChigh_step2,
              BIClow_step3 = BIClow_step3, BIChigh_step3 = BIChigh_step3, AIC_step3 = AIC_step3,
              ICL_BIClow_step3 = ICL_BIClow_step3, ICL_BIChigh_step3 = ICL_BIChigh_step3))
}
#
sel_other_poly = function(mY,mDesign,id_high,iT_range,iM_range,ivItemcat,approach,verbose,kmea,maxIter,tol,reord){
  iH      = sum(ivItemcat-1)
  K       = iH
  vNj     = table(id_high-1)
  iN      = length(vNj)
  iT_max  = max(iT_range)
  iT_min  = min(iT_range)
  num_iT  = iT_max-iT_min+1
  if(!is.null(iM_range)){
    iM_max  = max(iM_range)
    iM_min  = min(iM_range)
    num_iM  = iM_max-iM_min+1
    if(num_iM==1){
      if(iM_range==1){
        iM_range=NULL
      }
    }
  }
  incomplete = !is.null(mDesign)
  #
  if(approach=="model selection on low"|approach=="model selection on low with high"){
    BIClow      = rep(NA,num_iT)
    BIChigh     = rep(NA,num_iT)
    AIC         = rep(NA,num_iT)
    ICL_BIClow  = rep(NA,num_iT)
    ICL_BIChigh = rep(NA,num_iT)
    if(is.null(iM_range)){
      outmuLCA = list()
      for(i in iT_range){
        #
        if(verbose){
          pb=txtProgressBar(char="Fitting measurement model...",width=1)
          setTxtProgressBar(pb,1)
          close(pb)
        }
        #
        if(i==1){
          outmuLCA[[1]]   = NA
          llfoo           = sum(apply(mY,2,function(x){dbinom(x,1,mean(x,na.rm=T),log=T)}),na.rm=T)
          BIClow[1]       = -2*llfoo + K*log(sum(vNj))
          if(is.null(iM_range)) BIClow[1] = -2*llfoo + K*log(nrow(mY))
          if(!is.null(iM_range)) BIChigh[1] = -2*llfoo + K*log(iN)
          AIC[1]          = -2*llfoo + 2*K
          ICL_BIClow[1]   = BIClow[1]
          ICL_BIChigh[1]  = BIChigh[1]
        } else{
          outmuLCA[[1+i-iT_min]]  = LCA_fast_init_poly(mY,i,ivItemcat,incomplete,kmea,maxIter,tol,reord)
          ll                      = tail(outmuLCA[[1+i-iT_min]]$LLKSeries,1)
          npar                    = i*iH + i - 1
          BIClow[1+i-iT_min]      = outmuLCA[[1+i-iT_min]]$BIC
          BIChigh[1+i-iT_min]     = -2*ll + npar*log(iN)
          AIC[1+i-iT_min]         = -2*ll + 2*npar
          ICL_BIClow[1+i-iT_min]  = BIClow[i] + 2*sum(outmuLCA[[1+i-iT_min]]$freq*apply(-outmuLCA[[1+i-iT_min]]$mU*log(outmuLCA[[1+i-iT_min]]$mU),1,sum))
          ICL_BIChigh[1+i-iT_min] = Inf
        }
      }
      iT_best   = which(BIClow==min(BIClow,na.rm=T))+iT_min-1
      outmuLCA  = outmuLCA[[1+iT_best-iT_min]]
    } else{
      outmuLCA  = list()
      for(i in iT_range){
        #
        if(verbose){
          pb=txtProgressBar(char="Fitting measurement model...",width=1)
          setTxtProgressBar(pb,1)
          close(pb)
        }
        #
        if(i>1){
          start     = meas_Init_poly(mY,id_high,vNj,iM_range,i,ivItemcat,incomplete,kmea,maxIter,tol,reord)
          vOmegast  = start$vOmega_start
          mPi       = start$mPi_start
          if(!incomplete){
            outmuLCA[[1+i-iT_min]]  = MLTLCA_poly(mY,vNj,vOmegast,mPi,start$mPhi_start,ivItemcat,start$first_poly,1:i-1,1:iM_range-1,maxIter,tol,reord)
          } else{
            mY[is.na(mY)] = max(na.omit(mY))+1
            outmuLCA[[1+i-iT_min]]  = MLTLCA_poly_includeall(mY,mDesign,vNj,vOmegast,mPi,start$mPhi_start,ivItemcat,start$first_poly,1:i-1,1:iM_range-1,maxIter,tol,reord)
          }
          BIClow[1+i-iT_min]      = outmuLCA[[1+i-iT_min]]$BIClow
          BIChigh[1+i-iT_min]     = outmuLCA[[1+i-iT_min]]$BIChigh
          AIC[1+i-iT_min]         = outmuLCA[[1+i-iT_min]]$AIC
          ICL_BIClow[1+i-iT_min]  = outmuLCA[[1+i-iT_min]]$ICL_BIClow
          ICL_BIChigh[1+i-iT_min] = outmuLCA[[1+i-iT_min]]$ICL_BIChigh
        }
      }
      iT_best   = which(BIClow==min(BIClow,na.rm=T))+iT_min-1
      outmuLCA  = outmuLCA[[1+iT_best-iT_min]]
    }
    if(iT_best==1) stop("iT=1 optimal.",call.=FALSE)
    return(list(outmuLCA=outmuLCA,iT_best=iT_best,
                BIClow=BIClow,BIChigh=BIChigh,AIC=AIC,
                ICL_BIClow=ICL_BIClow,ICL_BIChigh=ICL_BIChigh))
  } else if(approach=="model selection on high"){
    BIClow      = rep(NA,num_iM)
    BIChigh     = rep(NA,num_iM)
    AIC         = rep(NA,num_iM)
    ICL_BIClow  = rep(NA,num_iM)
    ICL_BIChigh = rep(NA,num_iM)
    if(iT_range > 1){
      outmuLCA = list()
      for(i in iM_range){
        #
        if(verbose){
          pb=txtProgressBar(char="Fitting measurement model...",width=1)
          setTxtProgressBar(pb,1)
          close(pb)
        }
        #
        if(i==1){
          outmuLCA[[1]]   = LCA_fast_init_poly(mY,iT_range,ivItemcat,incomplete,kmea,maxIter,tol,reord)
          ll              = tail(outmuLCA[[1]]$LLKSeries,1)
          npar            = i*iH + i - 1
          BIClow[1]       = outmuLCA[[1]]$BIC
          BIChigh[1]      = -2*ll + npar*log(iN)
          AIC[1]          = -2*ll + 2*npar
          ICL_BIClow[1]   = BIClow[i] + 2*sum(outmuLCA[[1]]$freq*apply(-outmuLCA[[1]]$mU*log(outmuLCA[[1]]$mU),1,sum))
          ICL_BIChigh[1]  = Inf
        } else{
          start     = meas_Init_poly(mY,id_high,vNj,i,iT_range,ivItemcat,incomplete,kmea,maxIter,tol,reord)
          vOmegast  = start$vOmega_start
          mPi       = start$mPi_start
          if(!incomplete){
            outmuLCA[[1+i-iM_min]]    = MLTLCA_poly(mY,vNj,vOmegast,mPi,start$mPhi_start,ivItemcat,start$first_poly,1:iT_range-1,1:i-1,maxIter,tol,reord)
          } else{
            mY[is.na(mY)] = max(na.omit(mY))+1
            outmuLCA[[1+i-iM_min]]    = MLTLCA_poly_includeall(mY,mDesign,vNj,vOmegast,mPi,start$mPhi_start,ivItemcat,start$first_poly,1:iT_range-1,1:i-1,maxIter,tol,reord)
          }
          BIClow[1+i-iM_min]        = outmuLCA[[1+i-iM_min]]$BIClow
          BIChigh[1+i-iM_min]       = outmuLCA[[1+i-iM_min]]$BIChigh
          AIC[1+i-iM_min]           = outmuLCA[[1+i-iM_min]]$AIC
          ICL_BIClow[1+i-iM_min]    = outmuLCA[[1+i-iM_min]]$ICL_BIClow
          ICL_BIChigh[1+i-iM_min]   = outmuLCA[[1+i-iM_min]]$ICL_BIChigh
        }
      }
      iM_best   = which.min(BIChigh)+iM_min-1
      outmuLCA  = outmuLCA[[1+iM_best-iM_min]]
    } else{
      stop("iM=1 optimal.",call.=FALSE)
    }
    return(list(outmuLCA=outmuLCA,iM_best=iM_best,
                BIClow=BIClow,BIChigh=BIChigh,AIC=AIC,
                ICL_BIClow=ICL_BIClow,ICL_BIChigh=ICL_BIChigh))
  }
}
#
#
#
clean_output1 = function(output,Y,iT,iH,extout,dataout,ivItemcat){
  if(!extout){
    
    # Create empty list for output
    out = list()
    
    # Add vector of class proportions
    out$vPi           = output$pg
    rownames(out$vPi) = paste0("P(C",1:iT,")")
    colnames(out$vPi) = ""
    
    # Add matrix of conditional response probabilities
    out$mPhi            = output$mPhi
    rownames(out$mPhi)  = paste0("P(",Y,"|C)")
    colnames(out$mPhi)  = paste0("C",1:iT)
    
    # Add average proportion of classification errors
    out$AvgClassErrProb = output$dClassErr_tot
    
    # Add entropy R-sqr
    out$R2entr = output$R2entr
    
    # Add Bayesian information criterion
    out$BIC = output$BIC
    
    # Add Akaike information criterion
    out$AIC = output$AIC
    
    # Add number of iterations
    out$iter = output$iter
    
    # Add log-likelihood series
    out$LLKSeries = output$LLKSeries
    
    # Add specification
    out$spec            = as.matrix("Single-level LC model")
    rownames(out$spec)  = colnames(out$spec) = ""
    
  } else{
    
    Ylab_beta = unlist(lapply(ivItemcat,function(x){if(x==2){TRUE}else{c(FALSE,rep(TRUE,x-1))}}))
    
    # Create empty list for output
    out = list()
    
    # Add vector of class proportions
    out$vPi           = output$pg
    rownames(out$vPi) = paste0("P(C",1:iT,")")
    colnames(out$vPi) = ""
    
    # Add matrix of conditional response probabilities
    out$mPhi            = output$mPhi
    rownames(out$mPhi)  = paste0("P(",Y,"|C)")
    colnames(out$mPhi)  = paste0("C",1:iT)
    
    # Add matrix of posterior class membership probabilities
    out$mU            = output$mU[rep(1:dim(output$mU)[1],output$freq),]
    colnames(out$mU)  = colnames(out$mPhi)
    if(dataout){
      
      out$mU            = cbind(output$mY_unique[rep(1:dim(output$mY_unique)[1],output$freq),],out$mU)
      rownames(out$mU)  = NULL
      
    }
    
    # Add matrix of modal class assignment
    out$mU_modal            = output$mModalAssnm[rep(1:dim(output$mModalAssnm)[1],output$freq),]
    colnames(out$mU_modal)  = colnames(out$mPhi)
    if(dataout){
      
      out$mU_modal            = cbind(output$mY_unique[rep(1:dim(output$mY_unique)[1],output$freq),],out$mU_modal)
      rownames(out$mU_modal)  = NULL
      
    }
    
    # Add vector of modal class assignment
    out$vU_modal = output$vModalAssnm[rep(1:dim(output$vModalAssnm)[1],output$freq),,drop=FALSE]+1
    colnames(out$vU_modal) = "C"
    if(dataout){
      
      out$vU_modal            = cbind(output$mY_unique[rep(1:dim(output$mY_unique)[1],output$freq),],out$vU_modal)
      rownames(out$vU_modal)  = NULL
      
    }
    
    # Add matrix of classification errors
    out$mClassErr           = output$mClassErr
    rownames(out$mClassErr) = paste0("C",1:iT,"_true")
    colnames(out$mClassErr) = paste0("C",1:iT,"_pred")
    
    # Add matrix of proportion of classification errors
    out$mClassErrProb           = output$mClassErrProb
    rownames(out$mClassErrProb) = rownames(out$mClassErr)
    colnames(out$mClassErrProb) = colnames(out$mClassErr)
    
    # Add average proportion of classification errors
    out$AvgClassErrProb = output$dClassErr_tot
    
    # Add entropy R-sqr
    out$R2entr = output$R2entr
    
    # Add Bayesian information criterion
    out$BIC = output$BIC
    
    # Add Akaike information criterion
    out$AIC = output$AIC
    
    # Add gammas
    out$vGamma            = output$alphas
    rownames(out$vGamma)  = paste0("gamma(C",2:iT,")")
    colnames(out$vGamma)  = ""
    
    # Add betas
    out$mBeta           = output$gamma
    rownames(out$mBeta) = paste0("beta(",Y[Ylab_beta],"|C)")
    colnames(out$mBeta) = colnames(out$mPhi)
    
    # Add vector of model parameters
    out$parvec            = output$parvec
    rownames(out$parvec)  = c(rownames(out$vGamma),paste0(rep(substr(rownames(out$mBeta),1,nchar(rownames(out$mBeta))-2),iT),rep(colnames(out$mBeta),rep(nrow(out$mBeta),iT)),")"))
    colnames(out$parvec)  = ""
    
    # Add vector of standard errors
    out$SEs           = output$SEs
    rownames(out$SEs) = rownames(out$parvec)
    colnames(out$SEs) = ""
    
    # Add variance-covariance matrix
    out$Varmat            = output$Varmat
    rownames(out$Varmat)  = colnames(out$Varmat) = rownames(out$parvec)
    
    # Add number of iterations
    out$iter = output$iter
    
    # Add epsilon
    out$eps = output$eps
    
    # Add log-likelihood series
    out$LLKSeries = output$LLKSeries
    
    # Add matrix of model parameter contributions to log-likelihood score
    out$mScore            = output$mScore
    colnames(out$mScore)  = rownames(out$parvec)
    
    # Add specification
    out$spec            = as.matrix("Single-level LC model")
    rownames(out$spec)  = colnames(out$spec) = ""
    
  }
  return(out)
}
clean_output2 = function(output,Y,iT,Z,iH,P,extout,dataout,ivItemcat){
  mY      = output$mY
  mZ      = output$mZ
  output  = output$outcov
  if(!extout){
    
    # Create empty list for output
    out = list()
    
    # Add vector of sample means of class proportions
    out$vPi_avg           = as.matrix(apply(output$mPg,2,mean))
    rownames(out$vPi_avg) = paste0("P(C",1:iT,")")
    colnames(out$vPi_avg) = ""
    
    # Add matrix of conditional response probabilities
    out$mPhi            = output$mPhi
    rownames(out$mPhi)  = paste0("P(",Y,"|C)")
    colnames(out$mPhi)  = paste0("C",1:iT)
    
    # Add average proportion of classification errors
    out$AvgClassErrProb = output$dClassErr_tot
    
    # Add entropy R-sqr
    out$R2entr = output$R2entr
    
    # Add Bayesian information criterion
    out$BIC = output$BIC
    
    # Add Akaike information criterion
    out$AIC = output$AIC
    
    # Add gammas
    out$mGamma          = output$beta
    rownames(out$mGamma) = paste0("gamma(",Z,"|C)")
    colnames(out$mGamma) = paste0("C",2:iT)
    
    # Add corrected standard errors for gamma
    out$SEs_cor_gamma            = matrix(output$SEs_cor[1:((iT-1)*P),],P,iT-1)
    rownames(out$SEs_cor_gamma)  = rownames(out$mGamma)
    colnames(out$SEs_cor_gamma)  = colnames(out$mGamma)
    
    # Add number of iterations
    out$iter = output$iter
    
    # Add log-likelihood series
    out$LLKSeries = output$LLKSeries
    
    # Add specification
    out$spec            = as.matrix("Single-level LC model with covariates")
    rownames(out$spec)  = colnames(out$spec) = ""
    
  } else{
    
    Ylab_beta = unlist(lapply(ivItemcat,function(x){if(x==2){TRUE}else{c(FALSE,rep(TRUE,x-1))}}))
    
    # Create empty list for output
    out = list()
    
    # Add matrix of class proportions
    out$mPi           = output$mPg
    colnames(out$mPi) = paste0("P(C",1:iT,")")
    
    # Add vector of sample means of class proportions
    out$vPi_avg           = as.matrix(apply(output$mPg,2,mean))
    rownames(out$vPi_avg) = colnames(out$mPi)
    colnames(out$vPi_avg) = ""
    
    # Add matrix of conditional response probabilities
    out$mPhi            = output$mPhi
    rownames(out$mPhi)  = paste0("P(",Y,"|C)")
    colnames(out$mPhi)  = paste0("C",1:iT)
    
    # Add matrix of posterior class membership probabilities
    out$mU            = output$mU
    colnames(out$mU)  = colnames(out$mPhi)
    if(dataout){
      
      data              = cbind(mY,mZ)
      rownames(data)    = NULL
      colnames(data)    = c(Y,Z)
      data              = data[,-which(colnames(data)=="Intercept"),drop=FALSE]
      out$mU            = cbind(data,out$mU)
      rownames(out$mU)  = NULL
      
    }
    
    # Add matrix of classification errors
    out$mClassErr           = output$mClassErr
    rownames(out$mClassErr) = paste0("C",1:iT,"_true")
    colnames(out$mClassErr) = paste0("C",1:iT,"_pred")
    
    # Add matrix of proportion of classification errors
    out$mClassErrProb           = output$mClassErrProb
    rownames(out$mClassErrProb) = paste0("C",1:iT,"_true")
    colnames(out$mClassErrProb) = paste0("C",1:iT,"_pred")
    
    # Add average proportion of classification errors
    out$AvgClassErrProb = output$dClassErr_tot
    
    # Add entropy R-sqr
    out$R2entr = output$R2entr
    
    # Add Bayesian information criterion
    out$BIC = output$BIC
    
    # Add Akaike information criterion
    out$AIC = output$AIC
    
    # Add gammas
    out$mGamma          = output$beta
    rownames(out$mGamma) = paste0("gamma(",Z,"|C)")
    colnames(out$mGamma) = paste0("C",2:iT)
    
    # Add betas
    out$mBeta           = output$gamma
    rownames(out$mBeta) = paste0("beta(",Y[Ylab_beta],"|C)")
    colnames(out$mBeta) = colnames(out$mPhi)
    
    # Add vector of model parameters
    out$parvec            = output$parvec
    rownames(out$parvec)  = c(paste0(rep(substr(rownames(out$mGamma),1,nchar(rownames(out$mGamma))-2),iT-1),rep(colnames(out$mGamma),rep(P,iT-1)),")"),paste0(rep(substr(rownames(out$mBeta),1,nchar(rownames(out$mBeta))-2),iT),rep(colnames(out$mBeta),rep(nrow(out$mBeta),iT)),")"))
    colnames(out$parvec)  = ""
    
    # Add vector of uncorrected standard errors
    out$SEs_unc           = as.matrix(output$SEs_unc[1:((iT-1)*P),])
    rownames(out$SEs_unc) = rownames(out$mV2)
    colnames(out$SEs_unc) = ""
    
    # Add vector of corrected standard errors
    out$SEs_cor           = as.matrix(output$SEs_cor[1:((iT-1)*P),])
    rownames(out$SEs_cor) = rownames(out$mV2)
    colnames(out$SEs_cor) = ""
    
    # Add corrected standard errors for gamma
    out$SEs_cor_gamma            = matrix(output$SEs_cor[1:((iT-1)*P),],P,iT-1)
    rownames(out$SEs_cor_gamma)  = rownames(out$mGamma)
    colnames(out$SEs_cor_gamma)  = colnames(out$mGamma)
    
    # Add mQ
    out$mQ            = output$mQ
    rownames(out$mQ)  = colnames(out$mQ) = rownames(out$mV2)
    
    # Add uncorrected variance-covariance matrix
    out$Varmat_unc            = output$Varmat_unc[1:((iT-1)*P),1:((iT-1)*P)]
    rownames(out$Varmat_unc)  = colnames(out$Varmat_unc) = rownames(out$mV2)
    
    # Add corrected variance-covariance matrix
    out$Varmat_cor            = output$Varmat_cor
    rownames(out$Varmat_cor)  = colnames(out$Varmat_cor) = rownames(out$mV2)
    
    # Add inverse of the information matrix from the second step
    out$mV2           = output$mV2
    rownames(out$mV2) = colnames(out$mV2) = paste0(rep(substr(rownames(out$mGamma),1,nchar(rownames(out$mGamma))-2),iT-1),rep(colnames(out$mGamma),rep(P,iT-1)),")")
    
    # Add number of iterations
    out$iter = output$iter
    
    # Add epsilon
    out$eps = output$eps
    
    # Add log-likelihood series
    out$LLKSeries = output$LLKSeries
    
    # Add specification
    out$spec            = as.matrix("Single-level LC model with covariates")
    rownames(out$spec)  = colnames(out$spec) = ""
    
  }
  return(out)
}
clean_output3 = function(output,Y,iT,iM,mY,id_high,iH,id_high_levs,id_high_name,extout,dataout,ivItemcat){
  if(!extout){
    
    # Create empty list for output
    out = list()
    
    # Add matrix of high-level class proportions
    out$vOmega           = output$vOmega
    rownames(out$vOmega) = paste0("P(G",1:iM,")")
    colnames(out$vOmega) = ""
    
    # Add matrix of conditional low-level class proportions
    out$mPi           = output$mPi
    rownames(out$mPi) = paste0("P(C",1:iT,"|G)")
    colnames(out$mPi) = paste0("G",1:iM)
    
    # Add matrix of conditional response probabilities
    out$mPhi            = output$mPhi
    rownames(out$mPhi)  = paste0("P(",Y,"|C)")
    colnames(out$mPhi)  = paste0("C",1:iT)
    
    # Add low-level entropy R-sqr
    out$R2entr_low = output$R2entr_low
    
    # Add high-level entropy R-sqr
    out$R2entr_high = output$R2entr_high
    
    # Add low-level Bayesian information criterion
    out$BIClow = output$BIClow
    
    # Add high-level Bayesian information criterion
    out$BIChigh = output$BIChigh
    
    # Add low-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIClow = output$ICL_BIClow
    
    # Add high-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIChigh = output$ICL_BIChigh
    
    # Add Akaike information criterion
    out$AIC = output$AIC
    
    # Add number of iterations
    out$iter = output$iter
    
    # Add log-likelihood series
    out$LLKSeries = output$LLKSeries
    
    # Add specification
    out$spec            = as.matrix("Multilevel LC model")
    rownames(out$spec)  = colnames(out$spec) = ""
    
  } else{
    
    Ylab_beta = unlist(lapply(ivItemcat,function(x){if(x==2){TRUE}else{c(FALSE,rep(TRUE,x-1))}}))
    
    # Create empty list for output
    out = list()
    
    # Add matrix of high-level class proportions
    out$vOmega           = output$vOmega
    rownames(out$vOmega) = paste0("P(G",1:iM,")")
    colnames(out$vOmega) = ""
    
    # Add matrix of conditional low-level class proportions
    out$mPi           = output$mPi
    rownames(out$mPi) = paste0("P(C",1:iT,"|G)")
    colnames(out$mPi) = paste0("G",1:iM)
    
    # Add matrix of conditional response probabilities
    out$mPhi            = output$mPhi
    rownames(out$mPhi)  = paste0("P(",Y,"|C)")
    colnames(out$mPhi)  = paste0("C",1:iT)
    
    # Add cube of joint posterior class membership probabilities
    out$cPMX = array(output$cPMX,dim(output$cPMX),dimnames=list(NULL,paste0("C",1:iT,",G"),colnames(out$mPi)))
    if(dataout){
      
      data                = cbind(mY,id_high)
      rownames(data)      = NULL
      colnames(data)      = c(Y,id_high_name)
      for(i in 1:length(id_high_levs)){
        data[data[,id_high_name]==i,id_high_name] = id_high_levs[i]
      }
      out$cPMX            = array(unlist(lapply(1:dim(out$cPMX)[3],function(x){cbind(data,out$cPMX[,,x])})),
                                  c(dim(out$cPMX)[1],ncol(data)+ncol(out$cPMX),dim(out$cPMX)[3]),
                                  dimnames=list(dimnames(out$cPMX)[[1]],c(Y,id_high_name,dimnames(out$cPMX)[[2]]),dimnames(out$cPMX)[[3]]))
      
    }
    
    # Add cube of log of joint posterior class membership probabilities
    out$cLogPMX = array(output$cLogPMX,dim(output$cLogPMX),dimnames=list(NULL,paste0("C",1:iT,",G"),colnames(out$mPi)))
    if(dataout){
      
      out$cLogPMX = array(unlist(lapply(1:dim(out$cLogPMX)[3],function(x){cbind(data,out$cLogPMX[,,x])})),
                          c(dim(out$cLogPMX)[1],ncol(data)+ncol(out$cLogPMX),dim(out$cLogPMX)[3]),
                          dimnames=list(dimnames(out$cLogPMX)[[1]],c(Y,id_high_name,dimnames(out$cLogPMX)[[2]]),dimnames(out$cLogPMX)[[3]]))
      
    }
    
    # Add cube of conditional posterior low-level class membership probabilities
    out$cPX = array(output$cPX,dim(output$cPX),dimnames=list(NULL,paste0("C",1:iT,"|G"),colnames(out$mPi)))
    if(dataout){
      
      out$cPX = array(unlist(lapply(1:dim(out$cPX)[3],function(x){cbind(data,out$cPX[,,x])})),
                          c(dim(out$cPX)[1],ncol(data)+ncol(out$cPX),dim(out$cPX)[3]),
                          dimnames=list(dimnames(out$cPX)[[1]],c(Y,id_high_name,dimnames(out$cPX)[[2]]),dimnames(out$cPX)[[3]]))
      
    }
    
    # Add cube of log of conditional posterior low-level class membership probabilities
    out$cLogPX = array(output$cLogPX,dim(output$cLogPX),dimnames=list(NULL,paste0("C",1:iT,"|G"),colnames(out$mPi)))
    if(dataout){
      
      out$cLogPX = array(unlist(lapply(1:dim(out$cLogPX)[3],function(x){cbind(data,out$cLogPX[,,x])})),
                      c(dim(out$cLogPX)[1],ncol(data)+ncol(out$cLogPX),dim(out$cLogPX)[3]),
                      dimnames=list(dimnames(out$cLogPX)[[1]],c(Y,id_high_name,dimnames(out$cLogPX)[[2]]),dimnames(out$cLogPX)[[3]]))
      
    }
    
    # Add matrix of posterior high-level class membership probabilities for low-level units after marginalizing over low-level classes
    out$mSumPX = output$mSumPX
    colnames(out$mSumPX) = colnames(out$mPi)
    if(dataout){
      
      out$mSumPX            = cbind(data,out$mSumPX)
      rownames(out$mSumPX)  = NULL
      
    }
    
    # Add matrix of posterior high-level class membership probabilities for high-level units
    out$mPW           = output$mPW
    rownames(out$mPW) = id_high_levs
    colnames(out$mPW) = colnames(out$mPi)
    
    # Add matrix of log of posterior high-level class membership probabilities for high-level units
    out$mlogPW            = output$mlogPW
    rownames(out$mlogPW)  = id_high_levs
    colnames(out$mlogPW)  = colnames(out$mPi)
    
    # Add matrix of posterior high-level class membership probabilities for low-level units
    out$mPW_N           = output$mPW_N
    colnames(out$mPW_N) = colnames(out$mPi)
    if(dataout){
      
      out$mPW_N            = cbind(data,out$mPW_N)
      rownames(out$mPW_N)  = NULL
      
    }
    
    # Add matrix of posterior low-level class membership probabilities for low-level units after marginalizing over high-level classes
    out$mPMsumX           = output$mPMsumX
    colnames(out$mPMsumX) = colnames(out$mPhi)
    if(dataout){
      
      out$mPMsumX            = cbind(data,out$mPMsumX)
      rownames(out$mPMsumX)  = NULL
      
    }
    
    # Add low-level entropy R-sqr
    out$R2entr_low = output$R2entr_low
    
    # Add high-level entropy R-sqr
    out$R2entr_high = output$R2entr_high
    
    # Add low-level Bayesian information criterion
    out$BIClow = output$BIClow
    
    # Add high-level Bayesian information criterion
    out$BIChigh = output$BIChigh
    
    # Add low-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIClow = output$ICL_BIClow
    
    # Add high-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIChigh = output$ICL_BIChigh
    
    # Add Akaike information criterion
    out$AIC = output$AIC
    
    # Add alphas
    out$vAlpha            = output$vDelta
    rownames(out$vAlpha)  = paste0("alpha(G",2:iM,")")
    colnames(out$vAlpha)  = ""
    
    # Add gammas
    out$mGamma            = output$mGamma
    rownames(out$mGamma)  = paste0("gamma(C",2:iT,"|G)")
    colnames(out$mGamma)  = colnames(out$mPi)
    
    # Add betas
    out$mBeta           = output$mBeta
    rownames(out$mBeta) = paste0("beta(",Y[Ylab_beta],"|C)")
    colnames(out$mBeta) = colnames(out$mPhi)
    
    # Add vector of model parameters
    out$parvec            = output$parvec
    rownames(out$parvec)  = c(rownames(out$vAlpha),paste0(rep(substr(rownames(out$mGamma),1,nchar(rownames(out$mGamma))-2),iM),rep(colnames(out$mGamma),rep(iT-1,iM)),")"),paste0(rep(substr(rownames(out$mBeta),1,nchar(rownames(out$mBeta))-2),iT),rep(colnames(out$mBeta),rep(nrow(out$mBeta),iT)),")"))
    colnames(out$parvec)  = ""
    
    # Add vector of standard errors
    out$SEs           = output$SEs
    rownames(out$SEs) = rownames(out$parvec)
    colnames(out$SEs) = ""
    
    # Add variance-covariance matrix
    out$Varmat            = output$Varmat
    rownames(out$Varmat)  = colnames(out$Varmat) = rownames(out$parvec)
    
    # Add Fisher information matrix
    out$Infomat           = output$Infomat
    rownames(out$Infomat) = colnames(out$Infomat) = rownames(out$parvec)
    
    # Add number of iterations
    out$iter = output$iter
    
    # Add epsilon
    out$eps = output$eps
    
    # Add log-likelihood series
    out$LLKSeries = output$LLKSeries
    
    # Add current log-likelihood for high-level units
    out$vLLK            = output$vLLK
    rownames(out$vLLK)  = id_high_levs
    colnames(out$vLLK)  = ""
    
    # Add matrix of model parameter contributions to log-likelihood score
    out$mScore            = output$mScore
    colnames(out$mScore)  = rownames(out$parvec)
    
    # Add specification
    out$spec            = as.matrix("Multilevel LC model")
    rownames(out$spec)  = colnames(out$spec) = ""
    
  }
  return(out)
}
clean_output4 = function(output,Y,iT,iM,Z,mY,mZ,id_high,iH,P,id_high_levs,id_high_name,extout,dataout,ivItemcat){
  if(!extout){
    
    # Create empty list for output
    out = list()
    
    # Add matrix of high-level class proportions
    out$vOmega            = output$vOmega
    rownames(out$vOmega)  = paste0("P(G",1:iM,")")
    colnames(out$vOmega)  = ""
    
    # Add vector of sample means of conditional low-level class proportions
    out$mPi_avg           = apply(output$cPi,3,function(x){apply(x,2,mean)})
    rownames(out$mPi_avg) = paste0("P(C",1:iT,"|G)")
    colnames(out$mPi_avg) = paste0("G",1:iM)
    
    # Add matrix of conditional response probabilities
    out$mPhi            = output$mPhi
    rownames(out$mPhi)  = paste0("P(",Y,"|C)")
    colnames(out$mPhi)  = paste0("C",1:iT)
    
    # Add low-level entropy R-sqr
    out$R2entr_low = output$R2entr_low
    
    # Add high-level entropy R-sqr
    out$R2entr_high = output$R2entr_high
    
    # Add low-level Bayesian information criterion
    out$BIClow = output$BIClow
    
    # Add high-level Bayesian information criterion
    out$BIChigh = output$BIChigh
    
    # Add low-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIClow = output$ICL_BIClow
    
    # Add high-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIChigh = output$ICL_BIChigh
    
    # Add Akaike information criterion
    out$AIC = output$AIC
    
    # Add gammas
    if(is.null(output$mGamma_fixslope)){
      out$cGamma = array(apply(output$cGamma,3,t),c(P,iT-1,iM),dimnames=list(paste0("gamma(",Z,"|C)"),paste0("C",2:iT,",G"),colnames(out$mPi_avg)))
    } else{
      out$mGamma_fixslope = output$mGamma_fixslope
      rownames(out$mGamma_fixslope) = c(paste0("gamma(",Z[1],"_C)+gamma(",Z[1],"|G",1:iM,")"),paste0("gamma(",Z[-1],"|C)"))
      colnames(out$mGamma_fixslope) = paste0("C",2:iT)
    }
    
    # Add corrected standard errors for gamma
    if(is.null(output$mGamma_fixslope)){
      out$SEs_cor_gamma = array(apply(array(output$SEs_cor[iM:(iM-1+(P*(iT-1)*iM))],c(iT-1,P,iM)),3,t),dim(out$cGamma),dimnames=dimnames(out$cGamma))
    } else{
      out$SEs_cor_gamma = matrix(output$SEs_cor[iM:(iM-1+(iM+P-1)*(iT-1))],iM+P-1,(iT-1),dimnames=dimnames(out$mGamma_fixslope))
    }
    
    # Add number of iterations
    out$iter = output$iter
    
    # Add log-likelihood series
    out$LLKSeries = output$LLKSeries
    
    # Add specification
    out$spec            = as.matrix("Multilevel LC model with lower-level covariates")
    rownames(out$spec)  = colnames(out$spec) = ""
    
  } else{
    
    Ylab_beta = unlist(lapply(ivItemcat,function(x){if(x==2){TRUE}else{c(FALSE,rep(TRUE,x-1))}}))
    
    # Create empty list for output
    out = list()
    
    # Add matrix of high-level class proportions
    out$vOmega            = output$vOmega
    rownames(out$vOmega)  = paste0("P(G",1:iM,")")
    colnames(out$vOmega)  = ""
    
    # Add matrix of conditional low-level class proportions
    out$mPi           = matrix(output$cPi,nrow(output$cPi),iT*iM)
    colnames(out$mPi) = paste0(rep(paste0("P(C",1:iT,"|"), iM),rep(paste0("G",1:iM,")"),rep(iT,iM)))
    
    # Add vector of sample means of conditional low-level class proportions
    out$mPi_avg           = apply(output$cPi,3,function(x){apply(x,2,mean)})
    rownames(out$mPi_avg) = paste0("P(C",1:iT,"|G)")
    colnames(out$mPi_avg) = paste0("G",1:iM)
    
    # Add matrix of conditional response probabilities
    out$mPhi            = output$mPhi
    rownames(out$mPhi)  = paste0("P(",Y,"|C)")
    colnames(out$mPhi)  = paste0("C",1:iT)
    
    # Add cube of joint posterior class membership probabilities
    out$cPMX = array(output$cPMX,dim(output$cPMX),dimnames=list(NULL,paste0("C",1:iT,",G"),colnames(out$mPi_avg)))
    if(dataout){
      
      data            = cbind(mY,id_high,mZ)
      rownames(data)  = NULL
      colnames(data)  = c(Y,id_high_name,Z)
      for(i in 1:length(id_high_levs)){
        data[data[,id_high_name]==i,id_high_name] = id_high_levs[i]
      }
      data            = data[,-which(colnames(data)=="Intercept"),drop=FALSE]
      out$cPMX        = array(unlist(lapply(1:dim(out$cPMX)[3],function(x){cbind(data,out$cPMX[,,x])})),
                              c(dim(out$cPMX)[1],ncol(data)+ncol(out$cPMX),dim(out$cPMX)[3]),
                              dimnames=list(dimnames(out$cPMX)[[1]],c(Y,id_high_name,Z[-1],dimnames(out$cPMX)[[2]]),dimnames(out$cPMX)[[3]]))
      
    }
    
    # Add cube of log of joint posterior class membership probabilities
    out$cLogPMX = array(output$cLogPMX, dim(output$cLogPMX), dimnames=list(NULL,paste0("C",1:iT,",G"),colnames(out$mPi_avg)))
    if(dataout){
      
      out$cLogPMX = array(unlist(lapply(1:dim(out$cLogPMX)[3],function(x){cbind(data,out$cLogPMX[,,x])})),
                          c(dim(out$cLogPMX)[1],ncol(data)+ncol(out$cLogPMX),dim(out$cLogPMX)[3]),
                          dimnames=list(dimnames(out$cLogPMX)[[1]],c(Y,id_high_name,Z[-1],dimnames(out$cLogPMX)[[2]]),dimnames(out$cLogPMX)[[3]]))
      
    }
    
    # Add cube of conditional posterior low-level class membership probabilities
    out$cPX = array(output$cPX, dim(output$cPX),dimnames=list(NULL,paste0("C",1:iT,"|G"),colnames(out$mPi_avg)))
    if(dataout){
      
      out$cPX = array(unlist(lapply(1:dim(out$cPX)[3],function(x){cbind(data,out$cPX[,,x])})),
                      c(dim(out$cPX)[1],ncol(data)+ncol(out$cPX),dim(out$cPX)[3]),
                      dimnames=list(dimnames(out$cPX)[[1]],c(Y,id_high_name,Z[-1],dimnames(out$cPX)[[2]]),dimnames(out$cPX)[[3]]))
      
    }
    
    # Add cube of log of conditional posterior low-level class membership probabilities
    out$cLogPX = array(output$cLogPX,dim(output$cLogPX),dimnames=list(NULL,paste0("C",1:iT,"|G"),colnames(out$mPi_avg)))
    if(dataout){
      
      out$cLogPX = array(unlist(lapply(1:dim(out$cLogPX)[3],function(x){cbind(data,out$cLogPX[,,x])})),
                         c(dim(out$cLogPX)[1],ncol(data)+ncol(out$cLogPX),dim(out$cLogPX)[3]),
                         dimnames=list(dimnames(out$cLogPX)[[1]],c(Y,id_high_name,Z[-1],dimnames(out$cLogPX)[[2]]),dimnames(out$cLogPX)[[3]]))
      
    }
    
    # Add matrix of posterior high-level class membership probabilities for low-level units after marginalizing over low-level classes
    out$mSumPX            = output$mSumPX
    colnames(out$mSumPX)  = colnames(out$mPi_avg)
    if(dataout){
      
      out$mSumPX            = cbind(data,out$mSumPX)
      rownames(out$mSumPX)  = NULL
      
    }
    
    # Add matrix of posterior high-level class membership probabilities for high-level units
    out$mPW           = output$mPW
    rownames(out$mPW) = id_high_levs
    colnames(out$mPW) = colnames(out$mPi_avg)
    
    # Add matrix of log of posterior high-level class membership probabilities for high-level units
    out$mlogPW            = output$mlogPW
    rownames(out$mlogPW)  = id_high_levs
    colnames(out$mlogPW)  = colnames(out$mPi_avg)
    
    # Add matrix of posterior high-level class membership probabilities for low-level units
    out$mPW_N           = output$mPW_N
    colnames(out$mPW_N) = colnames(out$mPi_avg)
    if(dataout){
      
      out$mPW_N            = cbind(data,out$mPW_N)
      rownames(out$mPW_N)  = NULL
      
    }
    
    # Add matrix of posterior low-level class membership probabilities for low-level units after marginalizing over high-level classes
    out$mPMsumX           = output$mPMsumX
    colnames(out$mPMsumX) = colnames(out$mPhi)
    if(dataout){
      
      out$mPMsumX            = cbind(data,out$mPMsumX)
      rownames(out$mPMsumX)  = NULL
      
    }
    
    # Add low-level entropy R-sqr
    out$R2entr_low = output$R2entr_low
    
    # Add high-level entropy R-sqr
    out$R2entr_high = output$R2entr_high
    
    # Add low-level Bayesian information criterion
    out$BIClow = output$BIClow
    
    # Add high-level Bayesian information criterion
    out$BIChigh = output$BIChigh
    
    # Add low-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIClow = output$ICL_BIClow
    
    # Add high-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIChigh = output$ICL_BIChigh
    
    # Add Akaike information criterion
    out$AIC = output$AIC
    
    # Add alphas
    out$vAlpha            = output$vDelta
    rownames(out$vAlpha)  = paste0("alpha(G",2:iM,")")
    colnames(out$vAlpha)  = ""
    
    # Add gammas
    if(is.null(output$mGamma_fixslope)){
      out$cGamma = array(apply(output$cGamma,3,t),c(P,iT-1,iM),dimnames=list(paste0("gamma(",Z,"|C)"),paste0("C",2:iT,",G"),colnames(out$mPi_avg)))
    } else{
      out$mGamma_fixslope = output$mGamma_fixslope
      rownames(out$mGamma_fixslope) = c(paste0("gamma(",Z[1],"_C)+gamma(",Z[1],"|G",1:iM,")"),paste0("gamma(",Z[-1],"|C)"))
      colnames(out$mGamma_fixslope) = paste0("C",2:iT)
    }
    
    # Add betas
    out$mBeta           = output$mBeta
    rownames(out$mBeta) = paste0("beta(",Y[Ylab_beta],"|C)")
    colnames(out$mBeta) = colnames(out$mPhi)
    
    # Add vector of model parameters
    out$parvec            = output$parvec
    if(is.null(output$mGamma_fixslope)){
      rownames(out$parvec)  = c(rownames(out$vAlpha),paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),",")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")"),paste0(rep(substr(rownames(out$mBeta),1,nchar(rownames(out$mBeta))-2),iT),rep(colnames(out$mBeta),rep(nrow(out$mBeta),iT)),")"))
    } else{
      rownames(out$parvec)  = c(rownames(out$vAlpha),unlist(lapply(2:iT,function(x){c(paste0("gamma(",Z[1],"_C",x,")+gamma(",Z[1],"|G",1:iM,")"),paste0("gamma(",Z[-1],"|C",x,")"))})),paste0(rep(substr(rownames(out$mBeta),1,nchar(rownames(out$mBeta))-2),iT),rep(colnames(out$mBeta),rep(nrow(out$mBeta),iT)),")"))
    }
    colnames(out$parvec)  = ""
    
    # Add vector of uncorrected standard errors
    out$SEs_unc           = output$SEs_unc
    dimnames(out$SEs_unc) = dimnames(out$parvec)
    
    # Add vector of corrected standard errors
    out$SEs_cor           = output$SEs_cor
    dimnames(out$SEs_cor) = dimnames(out$parvec)
    
    # Add corrected standard errors for gamma
    if(is.null(output$mGamma_fixslope)){
      out$SEs_cor_gamma = array(apply(array(output$SEs_cor[iM:(iM-1+(P*(iT-1)*iM))],c(iT-1,P,iM)),3,t),dim(out$cGamma),dimnames=dimnames(out$cGamma))
    } else{
      out$SEs_cor_gamma = matrix(output$SEs_cor[iM:(iM-1+(iM+P-1)*(iT-1))],iM+P-1,(iT-1),dimnames=dimnames(out$mGamma_fixslope))
    }
    
    # Add mQ
    out$mQ = output$mQ
    if(is.null(output$mGamma_fixslope)){
      rownames(out$mQ)  = colnames(out$mQ) = c(rownames(out$vAlpha),paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),",")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")"))
    } else{
      rownames(out$mQ)  = colnames(out$mQ) = c(rownames(out$vAlpha),unlist(lapply(2:iT,function(x){c(paste0("gamma(",Z[1],"_C",x,")+gamma(",Z[1],"|G",1:iM,")"),paste0("gamma(",Z[-1],"|C",x,")"))})))
    }
    
    # Add uncorrected variance-covariance matrix
    out$Varmat_unc            = output$Varmat
    rownames(out$Varmat_unc)  = colnames(out$Varmat_unc) = rownames(out$parvec)
    
    # Add corrected variance-covariance matrix
    out$Varmat_cor            = output$mVar_corr
    dimnames(out$Varmat_cor) = dimnames(out$mQ)
    
    # Add Fisher information matrix
    out$Infomat           = output$Infomat
    rownames(out$Infomat) = colnames(out$Infomat) = rownames(out$parvec)
    
    # Add cube of Fisher information matrix for gamma
    if(is.null(output$mGamma_fixslope)){
      out$cGamma_Info = array(output$cGamma_Info,dim(output$cGamma_Info),dimnames=list(paste0("gamma(",Z,"|C,G)"),paste0("gamma(",Z,"|C,G)"),paste0("C", 2:iT,rep(paste0(",G", 1:iM),rep(iT-1,iM)))))
    }
    
    # Add inverse of the information matrix from the second step
    out$mV2           = output$mV2
    dimnames(out$mV2) = dimnames(out$mQ)
    
    # Add number of iterations
    out$iter = output$iter
    
    # Add epsilon
    out$eps = output$eps
    
    # Add log-likelihood series
    out$LLKSeries = output$LLKSeries
    
    # Add current log-likelihood for high-level units
    out$vLLK            = output$vLLK
    rownames(out$vLLK)  = id_high_levs
    colnames(out$vLLK)  = ""
    
    # Add matrix of model parameter contributions to log-likelihood score
    out$mScore            = output$mScore
    colnames(out$mScore)  = rownames(out$parvec)
    
    # Add subset of matrix of model parameter contributions to log-likelihood score for gamma
    out$mGamma_Score = output$mGamma_Score
    if(is.null(output$mGamma_fixslope)){
      colnames(out$mGamma_Score) = paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),",")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")")
    } else{
      colnames(out$mGamma_Score) = unlist(lapply(2:iT,function(x){c(paste0("gamma(",Z[1],"_C",x,")+gamma(",Z[1],"|G",1:iM,")"),paste0("gamma(",Z[-1],"|C",x,")"))}))
    }
    
    # Add specification
    out$spec            = as.matrix("Multilevel LC model with lower-level covariates")
    rownames(out$spec)  = colnames(out$spec) = ""
    
  }
  return(out)
}
clean_output5 = function(output,Y,iT,iM,Z,Zh,mY,mZ,mZh,id_high,iH,P,P_high,id_high_levs,id_high_name,extout,dataout,ivItemcat){
  if(!extout){
    
    # Create empty list for output
    out = list()
    
    # Add matrix of sample means of high-level class proportions
    out$vOmega_avg            = as.matrix(apply(output$mOmega,2,mean))
    rownames(out$vOmega_avg)  = paste0("P(G",1:iM,")")
    colnames(out$vOmega_avg)  = ""
    
    # Add matrix of sample means of conditional low-level class proportions
    out$mPi_avg           = apply(output$cPi,3,function(x){apply(x,2,mean)})
    rownames(out$mPi_avg) = paste0("P(C",1:iT,"|G)")
    colnames(out$mPi_avg) = paste0("G",1:iM)
    
    # Add matrix of conditional response probabilities
    out$mPhi            = output$mPhi
    rownames(out$mPhi)  = paste0("P(",Y,"|C)")
    colnames(out$mPhi)  = paste0("C",1:iT)
    
    # Add low-level entropy R-sqr
    out$R2entr_low = output$R2entr_low
    
    # Add high-level entropy R-sqr
    out$R2entr_high = output$R2entr_high
    
    # Add low-level Bayesian information criterion
    out$BIClow = output$BIClow
    
    # Add high-level Bayesian information criterion
    out$BIChigh = output$BIChigh
    
    # Add low-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIClow = output$ICL_BIClow
    
    # Add high-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIChigh = output$ICL_BIChigh
    
    # Add Akaike information criterion
    out$AIC = output$AIC
    
    # Add alphas
    out$mAlpha            = t(output$mDelta)
    rownames(out$mAlpha)  = paste0("alpha(",Zh,"|G)")
    colnames(out$mAlpha)  = paste0("G",2:iM)
    
    # Add gammas
    if(is.null(output$mGamma_fixslope)){
      out$cGamma = array(apply(output$cGamma,3,t),c(P,iT-1,iM),dimnames=list(paste0("gamma(",Z,"|C)"),paste0("C",2:iT,",G"),colnames(out$mPi_avg)))
    } else{
      out$mGamma_fixslope = output$mGamma_fixslope
      rownames(out$mGamma_fixslope) = c(paste0("gamma(",Z[1],"_C)+gamma(",Z[1],"|G",1:iM,")"),paste0("gamma(",Z[-1],"|C)"))
      colnames(out$mGamma_fixslope) = paste0("C",2:iT)
    }
    
    # Add corrected standard errors for alpha
    out$SEs_cor_alpha = t(matrix(output$SEs_cor[1:(P_high*(iM-1))],iM-1,P_high))
    rownames(out$SEs_cor_alpha) = rownames(out$mAlpha)
    colnames(out$SEs_cor_alpha) = colnames(out$mAlpha)
    
    # Add corrected standard errors for gamma
    if(is.null(output$mGamma_fixslope)){
      out$SEs_cor_gamma = array(apply(array(output$SEs_cor[(1+P_high*(iM-1)):(P_high*(iM-1)+P*(iT-1)*iM)],c(iT-1,P,iM)),3,t),dim(out$cGamma),dimnames=dimnames(out$cGamma))
    } else{
      out$SEs_cor_gamma = matrix(output$SEs_cor[(1+P_high*(iM-1)):(P_high*(iM-1)+(iM+P-1)*(iT-1))],iM+P-1,(iT-1),dimnames=dimnames(out$mGamma_fixslope))
    }
    
    # Add number of iterations
    out$iter = output$iter
    
    # Add log-likelihood series
    out$LLKSeries = output$LLKSeries
    
    # Add specification
    out$spec            = as.matrix("Multilevel LC model with lower- and higher-level covariates")
    rownames(out$spec)  = colnames(out$spec) = ""
    
  } else{
    
    Ylab_beta = unlist(lapply(ivItemcat,function(x){if(x==2){TRUE}else{c(FALSE,rep(TRUE,x-1))}}))
    
    # Create empty list for output
    out = list()
    
    # Add matrix of high-level class proportions
    out$mOmega            = output$mOmega
    colnames(out$mOmega)  = paste0("P(G",1:iM,")")
    
    # Add matrix of sample means of high-level class proportions
    out$vOmega_avg            = as.matrix(apply(output$mOmega,2,mean))
    rownames(out$vOmega_avg)  = colnames(out$mOmega)
    colnames(out$vOmega_avg)  = ""
    
    # Add matrix of conditional low-level class proportions
    out$mPi           = matrix(output$cPi,nrow(output$cPi),iT*iM)
    colnames(out$mPi) = paste0(rep(paste0("P(C",1:iT,"|"),iM),rep(paste0("G",1:iM,")"),rep(iT,iM)))
    
    # Add matrix of sample means of conditional low-level class proportions
    out$mPi_avg           = apply(output$cPi,3,function(x){apply(x,2,mean)})
    rownames(out$mPi_avg) = paste0("P(C",1:iT,"|G)")
    colnames(out$mPi_avg) = paste0("G",1:iM)
    
    # Add matrix of conditional response probabilities
    out$mPhi            = output$mPhi
    rownames(out$mPhi)  = paste0("P(",Y,"|C)")
    colnames(out$mPhi)  = paste0("C",1:iT)
    
    # Add cube of joint posterior class membership probabilities
    out$cPMX = array(output$cPMX,dim(output$cPMX),dimnames=list(NULL,paste0("C",1:iT,",G"),paste0("G", 1:iM)))
    if(dataout){
      
      data_low            = cbind(mY,id_high,mZ)
      rownames(data_low)  = NULL
      colnames(data_low)  = c(Y,id_high_name,Z)
      for(i in 1:length(id_high_levs)){
        data_low[data_low[,id_high_name]==i,id_high_name] = id_high_levs[i]
      }
      data_low            = data_low[,-which(colnames(data_low)=="Intercept"),drop=FALSE]
      out$cPMX            = array(unlist(lapply(1:dim(out$cPMX)[3],function(x){cbind(data_low,out$cPMX[,,x])})),
                                  c(dim(out$cPMX)[1],ncol(data_low)+ncol(out$cPMX),dim(out$cPMX)[3]),
                                  dimnames=list(dimnames(out$cPMX)[[1]],c(Y,id_high_name,Z[-1],dimnames(out$cPMX)[[2]]),dimnames(out$cPMX)[[3]]))
      
    }
    
    # Add cube of log of joint posterior class membership probabilities
    out$cLogPMX = array(output$cLogPMX,dim(output$cLogPMX),dimnames=list(NULL,paste0("C",1:iT,",G"),paste0("G", 1:iM)))
    if(dataout){
      
      out$cLogPMX = array(unlist(lapply(1:dim(out$cLogPMX)[3],function(x){cbind(data_low,out$cLogPMX[,,x])})),
                          c(dim(out$cLogPMX)[1],ncol(data_low)+ncol(out$cLogPMX),dim(out$cLogPMX)[3]),
                          dimnames=list(dimnames(out$cLogPMX)[[1]],c(Y,id_high_name,Z[-1],dimnames(out$cLogPMX)[[2]]),dimnames(out$cLogPMX)[[3]]))
      
    }
    
    # Add cube of conditional posterior low-level class membership probabilities
    out$cPX = array(output$cPX,dim(output$cPX),dimnames=list(NULL,paste0("C", 1:iT,"|G"),paste0("G", 1:iM)))
    if(dataout){
      
      out$cPX = array(unlist(lapply(1:dim(out$cPX)[3],function(x){cbind(data_low,out$cPX[,,x])})),
                      c(dim(out$cPX)[1],ncol(data_low)+ncol(out$cPX),dim(out$cPX)[3]),
                      dimnames=list(dimnames(out$cPX)[[1]],c(Y,id_high_name,Z[-1],dimnames(out$cPX)[[2]]),dimnames(out$cPX)[[3]]))
      
    }
    
    # Add cube of log of conditional posterior low-level class membership probabilities
    out$cLogPX = array(output$cLogPX,dim(output$cLogPX),dimnames=list(NULL,paste0("C", 1:iT,"|G"),paste0("G", 1:iM)))
    if(dataout){
      
      out$cLogPX = array(unlist(lapply(1:dim(out$cLogPX)[3],function(x){cbind(data_low,out$cLogPX[,,x])})),
                         c(dim(out$cLogPX)[1],ncol(data_low)+ncol(out$cLogPX),dim(out$cLogPX)[3]),
                         dimnames=list(dimnames(out$cLogPX)[[1]],c(Y,id_high_name,Z[-1],dimnames(out$cLogPX)[[2]]),dimnames(out$cLogPX)[[3]]))
      
    }
    
    # Add matrix of posterior high-level class membership probabilities for low-level units after marginalizing over low-level classes
    out$mSumPX            = output$mSumPX
    colnames(out$mSumPX)  = colnames(out$mPi_avg)
    if(dataout){
      
      out$mSumPX            = cbind(data_low,out$mSumPX)
      rownames(out$mSumPX)  = NULL
      
    }
    
    # Add matrix of posterior high-level class membership probabilities for high-level units
    out$mPW           = output$mPW
    rownames(out$mPW) = id_high_levs
    colnames(out$mPW) = colnames(out$mPi_avg)
    if(dataout){
      
      data_high           = mZh
      rownames(data_high) = id_high_levs
      colnames(data_high) = Zh
      data_high           = data_high[,-which(colnames(data_high)=="Intercept"),drop=FALSE]
      out$mPW             = cbind(data_high,out$mPW)
      
    }
    
    # Add matrix of log of posterior high-level class membership probabilities for high-level units
    out$mlogPW            = output$mlogPW
    rownames(out$mlogPW)  = id_high_levs
    colnames(out$mlogPW)  = colnames(out$mPi_avg)
    if(dataout){
      
      out$mlogPW = cbind(data_high,out$mlogPW)
      
    }
    
    # Add matrix of posterior high-level class membership probabilities for low-level units
    out$mPW_N           = output$mPW_N
    colnames(out$mPW_N) = colnames(out$mPi_avg)
    if(dataout){
      
      out$mPW_N            = cbind(data_low,out$mPW_N)
      rownames(out$mPW_N)  = NULL
      
    }
    
    # Add matrix of posterior low-level class membership probabilities for low-level units after marginalizing over high-level classes
    out$mPMsumX           = output$mPMsumX
    colnames(out$mPMsumX) = colnames(out$mPhi)
    if(dataout){
      
      out$mPMsumX            = cbind(data_low,out$mPMsumX)
      rownames(out$mPMsumX)  = NULL
      
    }
    
    # Add low-level entropy R-sqr
    out$R2entr_low = output$R2entr_low
    
    # Add high-level entropy R-sqr
    out$R2entr_high = output$R2entr_high
    
    # Add low-level Bayesian information criterion
    out$BIClow = output$BIClow
    
    # Add high-level Bayesian information criterion
    out$BIChigh = output$BIChigh
    
    # Add low-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIClow = output$ICL_BIClow
    
    # Add high-level integrated completed likelihood Bayesian information criterion
    out$ICL_BIChigh = output$ICL_BIChigh
    
    # Add Akaike information criterion
    out$AIC = output$AIC
    
    # Add alphas
    out$mAlpha            = t(output$mDelta)
    rownames(out$mAlpha)  = paste0("alpha(",Zh,"|G)")
    colnames(out$mAlpha)  = paste0("G",2:iM)
    
    # Add gammas
    if(is.null(output$mGamma_fixslope)){
      out$cGamma = array(apply(output$cGamma,3,t),c(P,iT-1,iM),dimnames=list(paste0("gamma(",Z,"|C)"),paste0("C",2:iT,",G"),colnames(out$mPi_avg)))
    } else{
      out$mGamma_fixslope = output$mGamma_fixslope
      rownames(out$mGamma_fixslope) = c(paste0("gamma(",Z[1],"_C)+gamma(",Z[1],"|G",1:iM,")"),paste0("gamma(",Z[-1],"|C)"))
      colnames(out$mGamma_fixslope) = paste0("C",2:iT)
    }
    
    # Add betas
    out$mBeta           = output$mBeta
    rownames(out$mBeta) = paste0("beta(",Y[Ylab_beta],"|C)")
    colnames(out$mBeta) = colnames(out$mPhi)
    
    # Add vector of model parameters
    out$parvec            = output$parvec
    if(is.null(output$mGamma_fixslope)){
      rownames(out$parvec)  = c(paste0(rep(substr(rownames(out$mAlpha),1,nchar(rownames(out$mAlpha))-2),rep(iM-1,P_high)),rep(colnames(out$mAlpha),P_high),")"),paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),",")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")"),paste0(rep(substr(rownames(out$mBeta),1,nchar(rownames(out$mBeta))-2),iT),rep(colnames(out$mBeta),rep(nrow(out$mBeta),iT)),")"))
    } else{
      rownames(out$parvec)  = c(paste0(rep(substr(rownames(out$mAlpha),1,nchar(rownames(out$mAlpha))-2),rep(iM-1,P_high)),rep(colnames(out$mAlpha),P_high),")"),unlist(lapply(2:iT,function(x){c(paste0("gamma(",Z[1],"_C",x,")+gamma(",Z[1],"|G",1:iM,")"),paste0("gamma(",Z[-1],"|C",x,")"))})),paste0(rep(substr(rownames(out$mBeta),1,nchar(rownames(out$mBeta))-2),iT),rep(colnames(out$mBeta),rep(nrow(out$mBeta),iT)),")"))
    }
    colnames(out$parvec)  = ""
    
    # Add vector of uncorrected standard errors
    out$SEs_unc           = output$SEs_unc
    dimnames(out$SEs_unc) = dimnames(out$parvec)
    
    # Add vector of corrected standard errors
    out$SEs_cor           = output$SEs_cor
    dimnames(out$SEs_cor) = dimnames(out$parvec)
    
    # Add corrected standard errors for alpha
    out$SEs_cor_alpha = t(matrix(output$SEs_cor[1:(P_high*(iM-1))],iM-1,P_high))
    rownames(out$SEs_cor_alpha) = rownames(out$mAlpha)
    colnames(out$SEs_cor_alpha) = colnames(out$mAlpha)
    
    # Add corrected standard errors for gamma
    if(is.null(output$mGamma_fixslope)){
      out$SEs_cor_gamma = array(apply(array(output$SEs_cor[(1+P_high*(iM-1)):(P_high*(iM-1)+P*(iT-1)*iM)],c(iT-1,P,iM)),3,t),dim(out$cGamma),dimnames=dimnames(out$cGamma))
    } else{
      out$SEs_cor_gamma = matrix(output$SEs_cor[(1+P_high*(iM-1)):(P_high*(iM-1)+(iM+P-1)*(iT-1))],iM+P-1,(iT-1),dimnames=dimnames(out$mGamma_fixslope))
    }
    
    # Add mQ
    out$mQ = output$mQ
    if(is.null(output$mGamma_fixslope)){
      rownames(out$mQ)  = colnames(out$mQ) = c(paste0(rep(substr(rownames(out$mAlpha),1,nchar(rownames(out$mAlpha))-2),rep(iM-1,P_high)),rep(colnames(out$mAlpha),P_high),")"),paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),",")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")"))
    } else{
      rownames(out$mQ)  = colnames(out$mQ) = c(paste0(rep(substr(rownames(out$mAlpha),1,nchar(rownames(out$mAlpha))-2),rep(iM-1,P_high)),rep(colnames(out$mAlpha),P_high),")"),unlist(lapply(2:iT,function(x){c(paste0("gamma(",Z[1],"_C",x,")+gamma(",Z[1],"|G",1:iM,")"),paste0("gamma(",Z[-1],"|C",x,")"))})))
    }
    
    # Add uncorrected variance-covariance matrix
    out$Varmat_unc            = output$Varmat
    rownames(out$Varmat_unc)  = colnames(out$Varmat_unc) = rownames(out$parvec)
    
    # Add corrected variance-covariance matrix
    out$Varmat_cor            = output$mVar_corr
    dimnames(out$Varmat_cor) = dimnames(out$mQ)
    
    # Add Fisher information matrix
    out$Infomat           = output$Infomat
    rownames(out$Infomat) = colnames(out$Infomat) = rownames(out$parvec)
    
    # Add cube of Fisher information matrix for alpha
    out$cAlpha_Info = array(output$cDelta_Info, dim(output$cDelta_Info), dimnames=list(paste0("alpha(",Zh,"|G)"),paste0("alpha(",Zh,"|G)"),colnames(out$mAlpha)))
    
    # Add cube of Fisher information matrix for gamma
    if(is.null(output$mGamma_fixslope)){
      out$cGamma_Info = array(output$cGamma_Info, dim(output$cGamma_Info), dimnames=list(paste0("gamma(",Z,"|C,G)"),paste0("gamma(",Z,"|C,G)"),paste0("C",2:iT,rep(paste0(",G",1:iM),rep(iT-1,iM)))))
    }
    
    # Add inverse of the information matrix from the second step
    out$mV2           = output$mV2
    dimnames(out$mV2) = dimnames(out$mQ)
    
    # Add number of iterations
    out$iter = output$iter
    
    # Add epsilon
    out$eps = output$eps
    
    # Add log-likelihood series
    out$LLKSeries = output$LLKSeries
    
    # Add current log-likelihood for high-level units
    out$vLLK            = output$vLLK
    rownames(out$vLLK)  = id_high_levs
    colnames(out$vLLK)  = ""
    
    # Add matrix of model parameter contributions to log-likelihood score
    out$mScore            = output$mScore
    colnames(out$mScore)  = rownames(out$parvec)
    
    # Add subset of matrix of model parameter contributions to log-likelihood score for alpha
    out$mAlpha_Score            = output$mDelta_Score
    rownames(out$mAlpha_Score)  = id_high_levs
    colnames(out$mAlpha_Score)  = paste0(rep(substr(rownames(out$mAlpha),1,nchar(rownames(out$mAlpha))-2),rep(iM-1,P_high)),rep(colnames(out$mAlpha),P_high),")")
    
    # Add subset of matrix of model parameter contributions to log-likelihood score for gamma
    out$mGamma_Score = output$mGamma_Score
    if(is.null(output$mGamma_fixslope)){
      colnames(out$mGamma_Score) = paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),",")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")")
    } else{
      colnames(out$mGamma_Score) = unlist(lapply(2:iT,function(x){c(paste0("gamma(",Z[1],"_C",x,")+gamma(",Z[1],"|G",1:iM,")"),paste0("gamma(",Z[-1],"|C",x,")"))}))
    }
    
    # Add specification
    out$spec            = as.matrix("Multilevel LC model with lower- and higher-level covariates")
    rownames(out$spec)  = colnames(out$spec) = ""
    
  }
  return(out)
}
#
clean_cov = function(data,covnam){
  cov = data[,covnam,drop=FALSE]
  if(any(!sapply(cov,is.numeric))){
    covnam_upd = c()
    for(i in covnam){
      if(!sapply(cov,is.numeric)[which(colnames(cov)==i)]){
        covnam_lev            = sort(unique(cov[,i]))
        cov_missing           = is.na(as.numeric(factor(cov[,i])))
        cov_upd               = vecTomatClass(as.numeric(factor(cov[,i])))
        cov_upd[cov_missing,] = NA
        colnames(cov_upd)     = paste0(i,".",covnam_lev)
        cov_upd               = cov_upd[,-1,drop=FALSE]
        covnam_upd            = c(covnam_upd,colnames(cov_upd))
      } else{
        cov_upd     = cov[,i,drop=FALSE]
        covnam_upd  = c(covnam_upd,i)
      }
      cov = cov[,-which(colnames(cov)==i),drop=FALSE]
      cov = cbind(cov,cov_upd)
    }
    covnam        = covnam_upd
    colnames(cov) = covnam
  }
  cov = as.matrix(cbind(1,cov))
  return(list(cov=cov,covnam=covnam))
}
#
update_YmY = function(mY,Y,ivItemcat,mDesign){
  Y_upd = c()
  for(i in Y){
    if(ivItemcat[i]>2){
      mY_upd            = vecTomatClass(mY[,i]+1)
      colnames(mY_upd)  = paste0(i,".",1:ncol(mY_upd)-1)
      Y_upd             = c(Y_upd,colnames(mY_upd))
    } else{
      mY_upd  = mY[,i,drop=FALSE]
      Y_upd   = c(Y_upd,i)
    }
    if(!is.null(mDesign)){
      row_missing = mDesign[,i]==0
      mY_upd[row_missing,] = NA
    }
    mY = mY[,-which(colnames(mY)==i)]
    mY = cbind(mY,mY_upd)
  }
  Y             = Y_upd
  colnames(mY)  = Y
  return(list(mY=mY,Y=Y))
}
#
update_YmY_nonnum = function(mY,Y,ivItemcat,mDesign){
  itemnam = lapply(data.frame(mY),function(x){levels(factor(x))})
  mY      = apply(mY,2,function(x){as.numeric(factor(x))-1})
  Y_upd   = c()
  for(i in Y){
    if(ivItemcat[i]>2){
      mY_upd            = vecTomatClass(mY[,i]+1)
      colnames(mY_upd)  = paste0(i,".",itemnam[[i]])
      Y_upd             = c(Y_upd,colnames(mY_upd))
    } else{
      mY_upd  = mY[,i,drop=FALSE]
      Y_upd   = c(Y_upd,paste0(i,".",itemnam[[i]][2]))
    }
    if(!is.null(mDesign)){
      row_missing = mDesign[,i]==0
      mY_upd[row_missing,] = NA
    }
    mY = mY[,-which(colnames(mY)==i)]
    mY = cbind(mY,mY_upd)
  }
  Y             = Y_upd
  colnames(mY)  = Y
  return(list(mY=mY,Y=Y))
}
#
check_inputs1 = function(data,Y,iT,id_high,iM,Z,Zh,startval){
  
  # <data> is matrix or dataframe
  if(!(is.matrix(data)|is.data.frame(data))){
    
    stop("data must be matrix or dataframe.",call.=FALSE)
    
  }
  
  # <Y>, <id_high>, <Z>, <Zh> and <startval> are column names in <data>
  if(any(!c(Y,id_high,Z,Zh,startval)%in%colnames(data))){
    
    if(any(!Y%in%colnames(data))){
      
      stop("At least one element of Y not found in data.",call.=FALSE)
      
    }
    if(!is.null(id_high)&any(!id_high%in%colnames(data))){
      
      stop("id_high not found in data.",call.=FALSE)
      
    }
    if(!is.null(Z)&any(!Z%in%colnames(data))){
      
      stop("At least one element of Z not found in data.",call.=FALSE)
      
    }
    if(!is.null(Zh)&any(!Zh%in%colnames(data))){
      
      stop("At least one element of Zh not found in data.",call.=FALSE)
      
    }
    if(!is.null(startval)&any(!startval%in%colnames(data))){
      
      stop("startval not found in data.",call.=FALSE)
      
    }
    
  }
  
}
#
check_inputs2 = function(data,Y,iT,id_high,iM,Z,Zh,startval){
  
  # Single <id_high>
  if(!is.null(id_high)){
    
    if(length(id_high)!=1){
      
      stop("Invalid id_high.",call.=FALSE)
      
    }
    
  }
  
  # Single <startval>
  if(!is.null(startval)){
    
    if(length(startval)!=1){
      
      stop("Invalid startval",call.=FALSE)
      
    }
    
  }
  
  # <iT> positive integer or sequential positive integers
  if(!length(iT)>=1){
    
    stop("Invalid iT.",call.=FALSE)
    
  } else{
    
    if(length(iT)==1){
      
      if(!is.numeric(iT)){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(iT!=round(iT)){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(iT==1){
        
        stop("Invalid iT: not a mixture model.",call.=FALSE)
        
      } else if(iT<1){
        
        stop("Invalid iT.",call.=FALSE)
        
      }
      
    } else{
      
      if(!is.vector(iT)|is.list(iT)){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(!is.numeric(iT)){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(any(!iT==round(iT))){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(min(iT)<1){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(length(iT)!=length(min(iT):max(iT))){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(any(!iT==min(iT):max(iT))){
        
        stop("Invalid iT.",call.=FALSE)
        
      }
      
    }
    
  }
  
  # <iM> positive integer or sequential positive integers
  if(!is.null(iM)&!length(iM)>=1){
    
    stop("Invalid iM.",call.=FALSE)
    
  } else if(!is.null(iM)){
    
    if(length(iM)==1){
      
      if(!is.numeric(iM)){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(iM!=round(iM)){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(iM<=1){
        
        stop("Invalid iM.",call.=FALSE)
        
      }
      
    } else{
      
      if(!is.vector(iM)|is.list(iM)){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(!is.numeric(iM)){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(any(!iM==round(iM))){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(min(iM)<1){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(length(iM)!=length(min(iM):max(iM))){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(any(!iM==min(iM):max(iM))){
        
        stop("Invalid iM.",call.=FALSE)
        
      }
      
    }
    
  }
  
  # <Y> consecutive integers from 0 or be non-numeric character
  if(any(sapply(as.data.frame(data[,Y]),function(x){(any(sort(na.omit(unique(x)))!=0:(length(na.omit(unique(x)))-1)))&!is.character(x)}))){  
    stop("Items must contain consecutive integers from 0 or be non-numeric character.",call.=FALSE)
    
  }
  
  # No missings in <id_high>
  if(any(is.na(data[,id_high]))){
    
    stop("id_high cannot contain missing values.",call.=FALSE)
    
  }
  
  # No missings in <startval>
  if(any(is.na(data[,startval]))){
    
    stop("startval cannot contain missing values.",call.=FALSE)
    
  }
  
  # <id_high> not labelled "Intercept"
  if(!is.null(id_high)){
    
    if(id_high=="Intercept"){
      
      stop("id_high may not be labelled 'Intercept'.",call.=FALSE)
      
    }
    
  }
  
  # No duplicate <Y>, <Z> or >Zh>
  if(any(duplicated(as.list(as.data.frame(data[,Y]))))){
    
    stop("Duplicate items.",call.=FALSE)
    
  }
  if(!is.null(Z)){
    
    if(length(Z)>=2){
      
      if(any(duplicated(as.list(as.data.frame(data[,Z]))))){
        
        stop("Duplicate lower-level covariates.",call.=FALSE)
        
      }
      
    }
    
  }
  if(!is.null(Zh)){
    
    if(length(Zh)>=2){
      
      if(any(duplicated(as.list(as.data.frame(data[,Zh]))))){
        
        stop("Duplicate higher-level covariates.",call.=FALSE)
        
      }
      
    }
    
  }
  
  # No constant <Y>, <Z> or <Zh>
  if(any(apply(data[,Y],2,function(x){length(na.omit(unique(x)))==1}))){  
    stop("Constant in indicators.",call.=FALSE)
    
  }
  if(!is.null(Z)){
    if(any(apply(data[,Z,drop=FALSE],2,function(x){length(na.omit(unique(x)))==1}))){
      
      stop("Constant in lower-level covariates.",call.=FALSE)
      
    }
    
  }
  if(!is.null(Zh)){
    if(any(apply(data[,Zh,drop=FALSE],2,function(x){length(na.omit(unique(x)))==1}))){  
      stop("Constant in higher-level covariates.",call.=FALSE)
      
    }
    
  }
  
  # Identify specification
  if(is.null(id_high)&is.null(iM)){
    
    specification = "single-level"
    
  } else if(!is.null(id_high)&!is.null(iM)){
    
    specification = "multilevel"
    
  } else{
    
    if(is.null(id_high)){
      
      stop("iM specified, id_high not.",call.=FALSE)
      
    }
    if(is.null(iM)){
      
      stop("id_high specified, iM not.",call.=FALSE)
      
    }
    
  }
  if(is.null(Z)){
    
    specification = paste(specification,"measurement")
    
    if(!is.null(Zh)){
      
      stop("Zh specified, Z not.",call.=FALSE)
      
    }
    
  } else if(!is.null(Z)){
    
    specification = paste(specification,"structural")
    
    if(specification=="single-level structural"&!is.null(Zh)){
      
      stop("Zh specified, iM and id_high not.",call.=FALSE)
      
    }
    
  }
  
  # Identify approach
  if(length(iT)==1&is.null(iM)){
    
    approach = "direct"
    
  } else if(length(iT)==1&length(iM)==1){
    
    approach = "direct"
    
  } else if(length(iT)>1&is.null(iM)){
    
    approach = "model selection on low"
    
  } else if(length(iT)>1&length(iM)==1){
    
    if(iM==1){
      
      approach = "model selection on low"
      
    } else{
      
      approach = "model selection on low with high"
      
    }
    
  } else if(length(iT)==1&length(iM)>1){
    
    approach = "model selection on high"
    
  } else if(length(iT)>1&length(iM)>1){
    
    approach = "model selection on low and high"
    
  }
  
  return(approach)
  
}
#
is_tibble_manual = function(x) {
  return(inherits(x, "tbl_df") && inherits(x, "tbl") && inherits(x, "data.frame"))
}
#
#
#
print.multiLCA = function(x,...){
  
  out = NULL
  
  if(x$spec == "Single-level LC model"){
    
    ######################################
    ### Single-level measurement model ###
    ######################################
    
    cat("\nCALL:\n")
    print(x$call)
    cat("\nSPECIFICATION:\n")
    print(noquote(x$spec))
    cat("\nMISSING VALUES ON THE INDICATORS:\n")
    print(noquote(x$missing_values))
    cat("\nFINAL SAMPLE SIZE:\n")
    print(x$sample_size)
    
    iter            = t(matrix(c(x$iter, x$LLKSeries[1], tail(x$LLKSeries, 1))))
    colnames(iter)  = c("EMiter", "LLfirst", "LLlast")
    rownames(iter)  = ""
    
    cat("\nESTIMATION DETAILS:\n\n")
    print(iter)
    cat("\n---------------------------\n")
    cat("\nCLASS PROPORTIONS:\n")
    print(round(x$vPi, 4))
    cat("\nRESPONSE PROBABILITIES:\n\n")
    print(round(x$mPhi, 4))
    cat("\n---------------------------\n")
    cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
    
    stat            = as.matrix(c(as.character(round(x$AvgClassErrProb, 4)), as.character(round(x$R2entr, 4)), as.character(round(x$BIC, 4)), as.character(round(x$AIC, 4))))
    rownames(stat)  = c("ClassErr", "EntR-sqr", "BIC", "AIC")
    colnames(stat)  = ""
    
    print(noquote(stat), right = TRUE)
    cat("\n")
    
  } else if(x$spec == "Single-level LC model with covariates"){
    
    #####################################
    ### Single-level structural model ###
    #####################################
    
    cat("\nCALL:\n")
    print(x$call)
    cat("\nSPECIFICATION:\n")
    print(noquote(x$spec))
    cat("\nMISSING VALUES ON THE INDICATORS:\n")
    print(noquote(x$missing_values))
    cat("\nFINAL SAMPLE SIZE:\n")
    print(x$sample_size)
    
    iter            = t(matrix(c(x$iter, x$LLKSeries[1], tail(x$LLKSeries, 1))))
    colnames(iter)  = c("EMiter", "LLfirst", "LLlast")
    rownames(iter)  = ""
    
    cat("\nESTIMATION DETAILS:\n\n")
    print(iter)
    cat("\n---------------------------\n")
    cat("\nCLASS PROPORTIONS (SAMPLE MEAN):\n")
    print(round(x$vPi_avg, 4))
    cat("\nRESPONSE PROBABILITIES:\n\n")
    print(round(x$mPhi, 4))
    cat("\n---------------------------\n")
    cat("\nMODEL AND CLASSIFICATION STATISTICS:\n\n")
    
    stat             = as.matrix(c(as.character(round(x$AvgClassErrProb, 4)), as.character(round(x$R2entr, 4)), as.character(round(x$BIC, 4)), as.character(round(x$AIC, 4))))
    rownames(stat)   = c("ClassErr", "EntR-sqr", "BIC", "AIC")
    colnames(stat)   = ""
    
    print(noquote(stat), right = TRUE)
    cat("\n")
    
    cat("\n---------------------------\n")
    cat("\nLOGISTIC MODEL FOR CLASS MEMBERSHIP:\n\n")
    cat("Estimator:",noquote(x$estimator),"\n\n")
    
    gammas  = format(round(as.matrix(x$mGamma), 4), nsmall = 4, scientific = FALSE)
    SE      = format(round(as.matrix(x$SEs_cor_gamma), 4), nsmall = 4, scientific = FALSE)
    Zscore  = format(round(as.matrix(x$mGamma)/as.matrix(x$SEs_cor_gamma), 4), nsmall = 4, scientific = FALSE)
    pval    = format(round(2*(1 - pnorm(abs(as.matrix(x$mGamma)/as.matrix(x$SEs_cor_gamma)))), 4), nsmall = 4, scientific = FALSE)
    psignf  = matrix("   ", nrow(gammas), ncol(gammas))
    psignf[pval < 0.1 & pval >= 0.05]   = "*  "
    psignf[pval < 0.05 & pval >= 0.01]  = "** "
    psignf[pval < 0.01]                 = "***"
    
    for (i in 1:ncol(gammas)){
      
      C = paste0("C", 1 + i)
      logit_params = noquote(cbind(gammas[,i], SE[,i], Zscore[,i], matrix(paste0(pval[,i], psignf[,i]))))
      colnames(logit_params) = c("Gamma", "S.E.", "Z-score", "p-value")
      rownames(logit_params) = paste0(substr(rownames(as.matrix(x$mGamma)), 1, nchar(rownames(as.matrix(x$mGamma))) - 1), 1 + i, ")")
      
      cat("\nMODEL FOR", C, "(BASE C1)\n\n")
      print(logit_params, right = TRUE)
      cat("\n")
      cat("\n", "*** p < 0.01, ** p < 0.05, * p < 0.1")
      cat("\n")
    }
    
  } else if(x$spec == "Multilevel LC model"){
    
    ####################################
    ### Multilevel measurement model ###
    ####################################
    
    cat("\nCALL:\n")
    print(x$call)
    cat("\nSPECIFICATION:\n")
    print(noquote(x$spec))
    cat("\nMISSING VALUES ON THE INDICATORS:\n")
    print(noquote(x$missing_values))
    cat("\nFINAL SAMPLE SIZE:\n")
    print(x$sample_size)
    
    iter            = t(matrix(c(x$iter, x$LLKSeries[1], tail(x$LLKSeries, 1))))
    colnames(iter)  = c("EMiter", "LLfirst", "LLlast")
    rownames(iter)  = ""
    
    cat("\nESTIMATION DETAILS:\n\n")
    print(iter)
    cat("\n---------------------------\n")
    cat("\nGROUP PROPORTIONS:\n")
    print(round(x$vOmega, 4))
    cat("\nCLASS PROPORTIONS:\n")
    print(round(x$mPi, 4))
    cat("\nRESPONSE PROBABILITIES:\n\n")
    print(round(x$mPhi, 4))
    cat("\n---------------------------\n")
    cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
    
    stat = as.matrix(c(as.character(round(x$R2entr_low, 4)), as.character(round(x$R2entr_high, 4)), as.character(round(x$BIClow, 4)), as.character(round(x$BIChigh, 4)), as.character(round(x$ICL_BIClow, 4)), as.character(round(x$ICL_BIChigh, 4)), as.character(round(x$AIC, 4))))
    rownames(stat) = c("R2entrlow", "R2entrhigh", "BIClow", "BIChigh", "ICLBIClow", "ICLBIChigh", "AIC")
    colnames(stat) = ""
    
    print(noquote(stat), right = TRUE)
    
  } else if(x$spec == "Multilevel LC model with lower-level covariates"){
    
    ###############################################################
    ### Multilevel structural model with lower-level covariates ###
    ###############################################################
    
    cat("\nCALL:\n")
    print(x$call)
    cat("\nSPECIFICATION:\n")
    print(noquote(x$spec))
    cat("\nMISSING VALUES ON THE INDICATORS:\n")
    print(noquote(x$missing_values))
    cat("\nFINAL SAMPLE SIZE:\n")
    print(x$sample_size)
    
    iter            = t(matrix(c(x$iter, x$LLKSeries[1], tail(x$LLKSeries, 1))))
    colnames(iter)  = c("EMiter", "LLfirst", "LLlast")
    rownames(iter)  = ""
    
    cat("\nESTIMATION DETAILS:\n\n")
    print(iter)
    cat("\n---------------------------\n")
    cat("\nGROUP PROPORTIONS:\n")
    print(round(x$vOmega, 4))
    cat("\nCLASS PROPORTIONS (SAMPLE MEAN):\n\n")
    print(round(x$mPi_avg, 4))
    cat("\nRESPONSE PROBABILITIES:\n\n")
    print(round(x$mPhi, 4))
    cat("\n---------------------------\n")
    cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
    
    stat = as.matrix(c(as.character(round(x$R2entr_low, 4)), as.character(round(x$R2entr_high, 4)), as.character(round(x$BIClow, 4)), as.character(round(x$BIChigh, 4)), as.character(round(x$ICL_BIClow, 4)), as.character(round(x$ICL_BIChigh, 4)), as.character(round(x$AIC, 4))))
    rownames(stat) = c("R2entrlow", "R2entrhigh", "BIClow", "BIChigh", "ICLBIClow", "ICLBIChigh", "AIC")
    colnames(stat) = ""
    
    print(noquote(stat), right = TRUE)
    cat("\n")
    
    cat("\n---------------------------\n")
    cat("\nLOGISTIC MODEL FOR LOWER-LEVEL CLASS MEMBERSHIP:\n\n")
    cat("Estimator:",noquote(x$estimator),"\n\n")
    if(!is.null(x$mGamma_fixslope)){
      cat("Fixed slopes across higher-level classes\n\n")
    }
    
    if(is.null(x$mGamma_fixslope)){
      for (i in 1:dim(x$cGamma)[3]){
        
        gammas  = format(round(as.matrix(x$cGamma[,,i]), 4), nsmall = 4, scientific = FALSE)
        SE      = format(round(as.matrix(x$SEs_cor_gamma[,,i]), 4), nsmall = 4, scientific = FALSE)
        Zscore  = format(round(as.matrix(x$cGamma[,,i])/as.matrix(x$SEs_cor_gamma[,,i]), 4), nsmall = 4, scientific = FALSE)
        pval    = format(round(2*(1 - pnorm(abs(as.matrix(x$cGamma[,,i])/as.matrix(x$SEs_cor_gamma[,,i])))), 4), nsmall = 4, scientific = FALSE)
        psignf  = matrix("   ", nrow(gammas), ncol(gammas))
        psignf[pval < 0.1 & pval >= 0.05]   = "*  "
        psignf[pval < 0.05 & pval >= 0.01]  = "** "
        psignf[pval < 0.01]                 = "***"
        
        for (j in 1:ncol(gammas)){
          
          G = paste0("G", i)
          C = paste0("C", 1 + j)
          logit_params = noquote(cbind(gammas[,j], SE[,j], Zscore[,j], matrix(paste0(pval[,j], psignf[,j]))))
          colnames(logit_params) = c("Gamma", "S.E.", "Z-score", "p-value")
          rownames(logit_params) = paste0(substr(rownames(as.matrix(x$cGamma[,,i])), 1, nchar(rownames(as.matrix(x$cGamma[,,i]))) - 1), j + 1, ",G", i, ")")
          
          cat("\nMODEL FOR", C, "(BASE C1) GIVEN", G, "\n\n")
          print(logit_params, right = TRUE)
          cat("\n")
          cat("\n", "*** p < 0.01, ** p < 0.05, * p < 0.1")
          cat("\n")
          
        }
        
      }
    } else{
      gammas  = format(round(as.matrix(x$mGamma_fixslope), 4), nsmall = 4, scientific = FALSE)
      SE      = format(round(as.matrix(x$SEs_cor_gamma), 4), nsmall = 4, scientific = FALSE)
      Zscore  = format(round(as.matrix(x$mGamma_fixslope)/as.matrix(x$SEs_cor_gamma), 4), nsmall = 4, scientific = FALSE)
      pval    = format(round(2*(1 - pnorm(abs(as.matrix(x$mGamma_fixslope)/as.matrix(x$SEs_cor_gamma)))), 4), nsmall = 4, scientific = FALSE)
      psignf  = matrix("   ", nrow(gammas), ncol(gammas))
      psignf[pval < 0.1 & pval >= 0.05]   = "*  "
      psignf[pval < 0.05 & pval >= 0.01]  = "** "
      psignf[pval < 0.01]                 = "***"
      
      iM = nrow(x$vOmega)
      
      for (i in 1:ncol(gammas)){
        
        C = paste0("C", 1 + i)
        logit_params = noquote(cbind(gammas[,i], SE[,i], Zscore[,i], matrix(paste0(pval[,i], psignf[,i]))))
        colnames(logit_params) = c("Gamma", "S.E.", "Z-score", "p-value")
        rownames(logit_params) = c(paste0(substr(rownames(as.matrix(x$mGamma_fixslope))[1:iM], 1, 17), 1 + i, substr(rownames(as.matrix(x$mGamma_fixslope))[1:iM], 18, nchar(rownames(as.matrix(x$mGamma_fixslope))[1:iM]))),
                                   paste0(substr(rownames(as.matrix(x$mGamma_fixslope))[-c(1:iM)], 1, nchar(rownames(as.matrix(x$mGamma_fixslope))[-c(1:iM)]) - 1), 1 + i, ")"))
        
        cat("\nMODEL FOR", C, "(BASE C1)\n\n")
        print(logit_params, right = TRUE)
        cat("\n")
        cat("\n", "*** p < 0.01, ** p < 0.05, * p < 0.1")
        cat("\n")
      }
    }
    
  } else if(x$spec == "Multilevel LC model with lower- and higher-level covariates"){
    
    ###########################################################################
    ### Multilevel structural model with lower- and higher-level covariates ###
    ###########################################################################
    
    cat("\nCALL:\n")
    print(x$call)
    cat("\nSPECIFICATION:\n")
    print(noquote(x$spec))
    cat("\nMISSING VALUES ON THE INDICATORS:\n")
    print(noquote(x$missing_values))
    cat("\nFINAL SAMPLE SIZE:\n")
    print(x$sample_size)
    
    iter            = t(matrix(c(x$iter, x$LLKSeries[1], tail(x$LLKSeries, 1))))
    colnames(iter)  = c("EMiter", "LLfirst", "LLlast")
    rownames(iter)  = ""
    
    cat("\nESTIMATION DETAILS:\n\n")
    print(iter)
    cat("\n---------------------------\n")
    cat("\nGROUP PROPORTIONS (SAMPLE MEAN):\n")
    print(round(x$vOmega_avg, 4))
    cat("\nCLASS PROPORTIONS (SAMPLE MEAN):\n\n")
    print(round(x$mPi_avg, 4))
    cat("\nRESPONSE PROBABILITIES:\n\n")
    print(round(x$mPhi, 4))
    cat("\n---------------------------\n")
    cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
    
    stat = as.matrix(c(as.character(round(x$R2entr_low, 4)), as.character(round(x$R2entr_high, 4)), as.character(round(x$BIClow, 4)), as.character(round(x$BIChigh, 4)), as.character(round(x$ICL_BIClow, 4)), as.character(round(x$ICL_BIChigh, 4)), as.character(round(x$AIC, 4))))
    rownames(stat) = c("R2entrlow", "R2entrhigh", "BIClow", "BIChigh", "ICLBIClow", "ICLBIChigh", "AIC")
    colnames(stat) = ""
    
    print(noquote(stat), right = TRUE)
    cat("\n")
    
    cat("\n---------------------------\n")
    cat("\nLOGISTIC MODEL FOR HIGHER-LEVEL CLASS MEMBERSHIP:\n\n")
    cat("Estimator:",noquote(x$estimator),"\n\n")
    
    alphas  = format(round(as.matrix(x$mAlpha), 4), nsmall = 4, scientific = FALSE)
    SE      = format(round(as.matrix(x$SEs_cor_alpha), 4), nsmall = 4, scientific = FALSE)
    Zscore  = format(round(as.matrix(x$mAlpha)/as.matrix(x$SEs_cor_alpha), 4), nsmall = 4, scientific = FALSE)
    pval    = format(round(2*(1 - pnorm(abs(as.matrix(x$mAlpha)/as.matrix(x$SEs_cor_alpha)))), 4), nsmall = 4, scientific = FALSE)
    psignf  = matrix("   ", nrow(alphas), ncol(alphas))
    psignf[pval < 0.1 & pval >= 0.05]   = "*  "
    psignf[pval < 0.05 & pval >= 0.01]  = "** "
    psignf[pval < 0.01]                 = "***"
    
    for (i in 1:ncol(alphas)){
      
      G = paste0("G", 1 + i)
      logit_params = noquote(cbind(alphas[,i], SE[,i], Zscore[,i], matrix(paste0(pval[,i], psignf[,i]))))
      colnames(logit_params) = c("Alpha", "S.E.", "Z-score", "p-value")
      rownames(logit_params) = paste0(substr(rownames(as.matrix(x$mAlpha)), 1, nchar(rownames(as.matrix(x$mAlpha))) - 1), 1 + i, ")")
      
      cat("\nMODEL FOR", G, "(BASE G1)\n\n")
      print(logit_params, right = TRUE)
      cat("\n")
      cat("\n", "*** p < 0.01, ** p < 0.05, * p < 0.1")
      cat("\n")
      
    }
    
    cat("\n---------------------------\n")
    cat("\nLOGISTIC MODEL FOR LOWER-LEVEL CLASS MEMBERSHIP:\n\n")
    cat("Estimator:",noquote(x$estimator),"\n\n")
    if(!is.null(x$mGamma_fixslope)){
      cat("Fixed slopes across higher-level classes\n\n")
    }
    
    if(is.null(x$mGamma_fixslope)){
      for (i in 1:dim(x$cGamma)[3]){
        
        gammas  = format(round(as.matrix(x$cGamma[,,i]), 4), nsmall = 4, scientific = FALSE)
        SE      = format(round(as.matrix(x$SEs_cor_gamma[,,i]), 4), nsmall = 4, scientific = FALSE)
        Zscore  = format(round(as.matrix(x$cGamma[,,i])/as.matrix(x$SEs_cor_gamma[,,i]), 4), nsmall = 4, scientific = FALSE)
        pval    = format(round(2*(1 - pnorm(abs(as.matrix(x$cGamma[,,i])/as.matrix(x$SEs_cor_gamma[,,i])))), 4), nsmall = 4, scientific = FALSE)
        psignf  = matrix("   ", nrow(gammas), ncol(gammas))
        psignf[pval < 0.1 & pval >= 0.05]   = "*  "
        psignf[pval < 0.05 & pval >= 0.01]  = "** "
        psignf[pval < 0.01]                 = "***"
        
        for (j in 1:ncol(gammas)){
          
          G = paste0("G", i)
          C = paste0("C", 1 + j)
          logit_params = noquote(cbind(gammas[,j], SE[,j], Zscore[,j], matrix(paste0(pval[,j], psignf[,j]))))
          colnames(logit_params) = c("Gamma", "S.E.", "Z-score", "p-value")
          rownames(logit_params) = paste0(substr(rownames(as.matrix(x$cGamma[,,i])), 1, nchar(rownames(as.matrix(x$cGamma[,,i]))) - 1), j + 1, ",G", i, ")")
          
          cat("\nMODEL FOR", C, "(BASE C1) GIVEN", G, "\n\n")
          print(logit_params, right = TRUE)
          cat("\n")
          cat("\n", "*** p < 0.01, ** p < 0.05, * p < 0.1")
          cat("\n")
          
        }
        
      }
    } else{
      gammas  = format(round(as.matrix(x$mGamma_fixslope), 4), nsmall = 4, scientific = FALSE)
      SE      = format(round(as.matrix(x$SEs_cor_gamma), 4), nsmall = 4, scientific = FALSE)
      Zscore  = format(round(as.matrix(x$mGamma_fixslope)/as.matrix(x$SEs_cor_gamma), 4), nsmall = 4, scientific = FALSE)
      pval    = format(round(2*(1 - pnorm(abs(as.matrix(x$mGamma_fixslope)/as.matrix(x$SEs_cor_gamma)))), 4), nsmall = 4, scientific = FALSE)
      psignf  = matrix("   ", nrow(gammas), ncol(gammas))
      psignf[pval < 0.1 & pval >= 0.05]   = "*  "
      psignf[pval < 0.05 & pval >= 0.01]  = "** "
      psignf[pval < 0.01]                 = "***"
      
      iM = nrow(x$vOmega)
      
      for (i in 1:ncol(gammas)){
        
        C = paste0("C", 1 + i)
        logit_params = noquote(cbind(gammas[,i], SE[,i], Zscore[,i], matrix(paste0(pval[,i], psignf[,i]))))
        colnames(logit_params) = c("Gamma", "S.E.", "Z-score", "p-value")
        rownames(logit_params) = c(paste0(substr(rownames(as.matrix(x$mGamma_fixslope))[1:iM], 1, 17), 1 + i, substr(rownames(as.matrix(x$mGamma_fixslope))[1:iM], 18, nchar(rownames(as.matrix(x$mGamma_fixslope))[1:iM]))),
                                   paste0(substr(rownames(as.matrix(x$mGamma_fixslope))[-c(1:iM)], 1, nchar(rownames(as.matrix(x$mGamma_fixslope))[-c(1:iM)]) - 1), 1 + i, ")"))
        
        cat("\nMODEL FOR", C, "(BASE C1)\n\n")
        print(logit_params, right = TRUE)
        cat("\n")
        cat("\n", "*** p < 0.01, ** p < 0.05, * p < 0.1")
        cat("\n")
      }
    }
    
  }
  
}
#
plot.multiLCA = function(x, horiz = FALSE, clab = NULL, ...){
  
  out = NULL
  
  iT        = ncol(x$mPhi)
  items     = 1:nrow(x$mPhi)
  itemnames = substr(rownames(x$mPhi), 3, nchar(rownames(x$mPhi)) - 3)
  
  if(horiz == TRUE){
    
    las = 1
    
  } else{
    
    las = 3
    
  }
  if(is.null(clab)){
    
    legend = paste0("Class ", 1:iT)
    
  } else{
    
    if(!is.vector(clab) | is.list(clab)){
      
      stop("Invalid clab.",call.=FALSE)
      
    } else if(length(clab) != iT | !is.character(clab)){
      
      stop("Invalid clab.",call.=FALSE)
      
    }
    
    legend = clab
    
  }
  
  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(5, 5, 3, 8), xpd = TRUE)
  plot(items, x$mPhi[,1], type = "b", xaxt = "n",yaxt = "n", pch = 1, 
       xlab = "", ylab = "Response probability", frame.plot = FALSE, ylim = c(0,1), lty = 1, ...)
  axis(1, at = items, labels = itemnames, las = las)
  axis(2, las = 1)
  for(h in 2:iT){
    lines(items, x$mPhi[,h], pch = h, type = "b", lty = h)
  }
  legend("topright", legend = legend, inset = c(-0.4,0),
         lty = (1:iT), pch = 1:iT, cex = 0.8)
  
}