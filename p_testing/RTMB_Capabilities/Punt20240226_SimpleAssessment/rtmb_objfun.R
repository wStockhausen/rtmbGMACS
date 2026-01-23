calcLogisticSel<-function(p1,p2){
  S <- 1.0/(1+exp(-log(19.0)*(Length-p1)/(p2-p1)));
  return(S);
}

calcFandZ<-function(logFullF,s){
  F<-matrix(0,nrow=Nyear,ncol=Nclass);#--NOTE: need to initialize values here
  Z<-matrix(0,nrow=Nyear,ncol=Nclass);#--NOTE: need to initialize values here
  for (Iyear in 1:Nyear){
    for (Ilen in 1:Nclass) {
       F[Iyear,Ilen] <- exp(logFullF[Iyear])*s[Ilen];
       Z[Iyear,Ilen] <- M + F[Iyear,Ilen];
    }
  }
  return(list(F=F,Z=Z));
}

calcN<-function(N0,F,Z,Rbar,Eps){
  N   <- matrix(0,nrow=Nyear+1,ncol=Nclass);#--NOTE: need to initialize values here
  CAL <- matrix(0,nrow=Nyear,  ncol=Nclass);#--NOTE: need to initialize values here
  CW <- rep(0,Nyear);
  for (Ilen in 1:Nclass) N[1,Ilen] = N0[Ilen];
  for (Iyear in 1:Nyear) {
    # Compute catch-at-length and -in-weight
    for (Ilen in 1:Nclass) {
      CAL[Iyear,Ilen] = F[Iyear,Ilen]/Z[Iyear,Ilen]*N[Iyear,Ilen]*(1.0-exp(-Z[Iyear,Ilen]));
      CW[Iyear] <- CW[Iyear] + CAL[Iyear,Ilen]*WL[Ilen];
     }#--Ilen
    CAL[Iyear,] <-CAL[Iyear,]/sum(CAL[Iyear,])

    # Update the dynamics
    for (Ilen in 1:Nclass){
      N[Iyear+1,Ilen] = 0;
      for (Jlen in 1:Nclass)
        N[Iyear+1,Ilen] <- N[Iyear+1,Ilen]+ N[Iyear,Jlen]*exp(-Z[Iyear,Jlen])*X[Jlen,Ilen];
     }#--Ilen
    # Add recruitment
    N[Iyear+1,1] <- N[Iyear+1,1] + Rbar*exp(Eps[Iyear]);
  }#--Iyear
  return(list(N=N,CW=CW,CAL=CAL));
}

objfun <- function(params){

  getAll(params)

  #--identify observations for simulation
  CWObs   <-OBS(CWObs);
  CALObs  <-OBS(CALObs);
  BioIndex<-OBS(BioIndex);

  # Extract parameters
  Rbar     <- exp(pars[1]);
  N0       <- exp(pars[2:(Nclass+1)]);
  LogFullF <- pars[(Nclass+2):(Nclass+1+Nyear)]
  Eps      <- pars[(Nclass+1+Nyear+1):(Nclass+1+2*Nyear)]
  ip = Nclass+1+2*Nyear+1;
  S50  <- exp(pars[ip]); ip = ip+1;
  S95  <- exp(pars[ip]); ip = ip+1;
  SS50 <- exp(pars[ip]); ip = ip+1;
  SS95 <- exp(pars[ip]); ip = ip+1;
  Sps  = rep(0,Nclass);
  for (i in 1: Nclass) {Sps[i] <- exp(pars[ip]); ip = ip+1;}
  SSps  = rep(0,Nclass);
  for (i in 1: Nclass) {SSps[i] <- exp(pars[ip]); ip = ip+1;}

  recs = Rbar*exp(Eps)
  ADREPORT(Rbar);
  ADREPORT(N0);
  ADREPORT(recs);
  ADREPORT(S50);
  ADREPORT(S95);
  ADREPORT(SS50);
  ADREPORT(SS95);
  ADREPORT(Sps);
  ADREPORT(SSps);


  # Fishery Selectivity
  if (useParFishSel){
    S = calcLogisticSel(S50,S95);#--values calculated internally need to be passed(?)
  } else {
    cat("Using non-parametric fishery selectivity\n")
    S = Sps;
  }
  ADREPORT(S);

  # F and Z
  lst = calcFandZ(LogFullF,S);#--values calculated internally need to be passed(?)
  F = lst$F;
  Z = lst$Z;
  ADREPORT(F);
  ADREPORT(Z);

  # Compute the N matrix
  lst = calcN(N0,F,Z,Rbar,Eps);#--values calculated internally need to be passed(?)
  N   = lst$N;
  CAL = lst$CAL;
  CW  = lst$CW;
  ADREPORT(N);
  ADREPORT(CW);
  ADREPORT(CAL);

  # Catch likelihood
  LikeCatch = -sum(dnorm(log(CWObs),log(CW),0.05,log=TRUE));
  # LikeCatch = -sum(dlnorm(CWObs,log(CW),0.05,log=TRUE));
  lnCWObs = log(CWObs);
  lnCW    = log(CW);
  lnCWObs %~% dnorm(lnCW,0.05); #--adds -log-likelihood to hidden variable `.nll`

  # Survey Selectivity
  if (useParSurvSel){
    SurveyS = calcLogisticSel(SS50,SS95);
  } else {
    cat("Using non-parametric fishery selectivity\n")
    SurveyS = SSps;
  }
  ADREPORT(SurveyS);

  # Biomass predictions
  BioPred <- rep(0,Nyear)
  for (Iyear in 1:Nyear)
   for (Ilen in 1:Nclass)
    BioPred[Iyear] <- BioPred[Iyear] + N[Iyear,Ilen]*SurveyS[Ilen]*WL[Ilen]
  ADREPORT(BioPred);

  #ML estimate of q
  Top = 0; Bot = 0;
  for (Iyear in 1:Nyear) {
    Top <- Top + log(BioIndex[Iyear]/BioPred[Iyear]);
    Bot <- Bot + 1.0;
  }
  q <- exp(Top/Bot);
  ADREPORT(q);

  # Likelihood
  LikeBio = -sum(dnorm(log(BioIndex),log(q*BioPred),BioSig,log=TRUE));
  lnBioIndex = log(BioIndex);
  lnBioIndex %~% dnorm(log(q*BioPred),BioSig);#--adds -log-likelihood to hidden variable `.nll`

  LikeCAL = 0;
  nEffCALObs = matrix(0,nrow=Nyear,ncol=Nclass);
  for (Iyear in 1:Nyear){
    LikeCAL = LikeCAL - dmultinom(Neff*CALObs[Iyear,],Neff,CAL[Iyear,],log=TRUE);
    nEffCALObs[Iyear,] = Neff*CALObs[Iyear,];
    nEffCALObs[Iyear,] %~% dmultinom(Neff,CAL[Iyear,]);#--adds -log-likelihood to hidden variable `.nll`
  }

  Penal = -sum(dnorm(Eps,0,0.6,log=TRUE));
  Eps %~% dnorm(0,0.6);#--adds -log-likelihood to hidden variable `.nll`
  ADREPORT(Eps);

  obj_fun <- LikeCatch + LikeCAL + LikeBio + Penal;
  REPORT(LikeCatch);
  REPORT(LikeCAL);
  REPORT(LikeBio);
  REPORT(Penal);
  REPORT(obj_fun);
  ADREPORT(LikeCatch);
  ADREPORT(LikeCAL);
  ADREPORT(LikeBio);
  ADREPORT(Penal);
  ADREPORT(obj_fun);

  ADREPORT(pars);

  #return(obj_fun);#--returns `.nll`
}
