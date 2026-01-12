#--test REs
require(ggplot2);
require(Matrix);
require(RTMB);
dirPrj = rstudioapi::getActiveProject();

#--helper functions----
source(file.path(dirPrj,"R/MiscFunctions_Matrices.R"));
logit<-function(x){log(x/(1-x))};
invlogit<-function(x){exp(x)/(1+exp(x))}

#--define fixed parameters
lstParsFEs = list(mu=10);

#--define covariance parameters related to random effects----
lstParsCMsREs = list(list(id=1,type="diag",n=10,pars=list(c(logsd=0,logit_rho=NA))),
               list(id=2,type="ar1", n=20,pars=list(c(logsd=0,logit_rho=logit(0.1)))));
dfrParsCMsREs = dplyr::bind_rows(lstParsCMsREs) |> dplyr::mutate(name=paste0("parCMREs_",dplyr::row_number()),.before=1);

#--create precision matrices for REs----
force=FALSE;
lstQs = list();
for (i in 1:nrow(dfrParsCMsREs)){
  #--testing: i = 2;
  #                      number of rows     parameter vector       process type
  lstQs[[i]] = mkPrecMat(dfrParsCMsREs$n[i],dfrParsCMsREs$pars[i][[1]],dfrParsCMsREs$type[i],force=force);
}
Qp = buildBlockDiagonalMat(lstQs,force=force);

#--calculate Cholesky factor----
Lp = Matrix::chol(Q); #--UPPER TRIANGULAR Cholesky factor used for simulation

#--create observations----
z   = rnorm(nrow(Q));                      #--random effects: iid normal random deviates
obs = lstParsFEs$mu + Matrix::solve(Lp,z); #--observed Normal process for mean mu,  precision matrix Q
lstObs = list();
for (i in 1:nrow(dfrParsCMsREs)) lstObs[[i]] = tibble::tibble(block=as.factor(rep(i,dfrParsCMsREs$n[i])));
dfrObs = dplyr::bind_cols(dplyr::bind_rows(lstObs),obs=obs) |> dplyr::mutate(i=dplyr::row_number(),.before=1);
ggplot(dfrObs,aes(x=i,y=obs,colour=block,group=block)) + geom_point() + geom_line();

#--set up parameters and input values for RTMB objective function----
#   obs[i] ~ mu + r[i]; mu: mean. r: random effects
#   r ~ N(0,Q^{-1});    Q: precision matrix
##--fixed effects parameters----
lstParamsFEs = list(mu=5);

##--random effects parameters
lstParamsREs = list(r=rep(0.0,length(z)));

##--covariance-related parameters related to random effects----
lstParamsCMsREs = NULL;
for (r in 1:nrow(dfrParsCMsREs)){
  lstParamsCMsREs[[dfrParsCMsREs$name[r]]] = dfrParsCMsREs$pars[r][[1]];
}
##--covariance-related parameters related to smooths----
lstParamsCMsSMs = NULL;
##--create combined parameter list
lstParams = c(lstParamsFEs,lstParamsREs,lstParamsCMsREs,lstParamsCMsSMs);

##--create parameter `map` to handle fixed (and mirrored) values
map = list();
for (p in names(lstParams)) {
  #--testing: p = 3;
  map[[p]] = 1:(length(lstParams[[p]]));
  nas = is.na(lstParams[[p]]);
  map[[p]][nas] = NA;
  map[[p]] = factor(map[[p]]);
}

##--define limits on parameters?----

#--define objective function----
objfun<-function(params){
  mu = params$mu; #--might want to make this a vector (needed as vector for dgmrf)
  ##--calculate precision matrix (REs + smooths)----
  lstQs = list();
  ###-calculate precision matrices for REs----
  if (exists("dfrParsCMsREs")){
    for (i in 1:nrow(dfrParsCMsREs)){
      #--testing: i = 2;
      nm = dfrParsCMsREs$name[i]
      lstQs[[nm]] = mkPrecMat(dfrParsCMsREs$n[i],params[[nm]],dfrParsCMsREs$type[i]);
    }
  }
  ###--calculate precision matrices for SMs----
  if (exists("dfrParsCMsSMs")){
    for (i in 1:nrow(dfrParsCMsSMs)){
      #--testing: i = 2;
      nm = dfrParsCMsSMs$name[i]
      lstQs[[nm]] = mkPrecMat(dfrParsCMsSMs$n[i],params[[name]],dfrParsCMsREs$type[i]);
    }
  }

  Q = buildBlockDiagonalMat(lstQs);

  #--predict observations
  prd = mu + params$r;

  #--calculate negative log-likelihood
  nll = RTMB::AD(0.0);
  nll = nll-RTMB::dgmrf(obs,rep(mu,length(obs)),Q,log=TRUE,scale=1);

  RTMB::ADREPORT(prd)
  RTMB::REPORT(Q)

  return(nll)
}

#--test objfun----
objfun(lstParams);


#--do estimation----
require(Matrix);
require(RTMB);
testing=FALSE;
objfun<-MakeADFun(objfun,
                  params,
                  random=names(lstParamsREs),
                  map=map,
                  hessian=TRUE,
                  LaplaceNonZeroGradient=TRUE,
                  ADreport=FALSE,
                  ridge.correct=FALSE, #--Note: apparently EXPERIMENTAL
                  silent=FALSE);
opt <- nlminb(objfun$par, objfun$fn, objfun$gr,objfun$he);
sdrep <- sdreport(objfun);
print(sdrep);
mPrd = summary(sdrep,"report");
