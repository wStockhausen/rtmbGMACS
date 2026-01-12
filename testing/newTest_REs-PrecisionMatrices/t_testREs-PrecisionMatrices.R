#--test REs
require(ggplot2);
require(Matrix);
require(RTMB);
dirPrj = rstudioapi::getActiveProject();

#--helper functions----
source(file.path(dirPrj,"R/MiscFunctions_Matrices.R"));
symlogit<-function(x){y = (x+1)/2.0; return(log(y/(1-y)))};
syminvlogit<-function(x){2*exp(x)/(1+exp(x))-1.0}

#--define fixed parameters
lstParsFEs = list(mu=10);

#--define covariance parameters related to random effects----
lstParsCMsREs = list(list(id=1,type="diag",n=10,pars=list(c(logsd=-1,logit_rho=0))),
                     list(id=2,type="ar1", n=20,pars=list(c(logsd=-1,logit_rho=symlogit(0))))
                    );
dfrParsCMsREs = dplyr::bind_rows(lstParsCMsREs) |> dplyr::mutate(name=paste0("parCMREs_",dplyr::row_number()),.before=1);

#--create precision matrices for REs----
force=FALSE; #--for testing RTMB AD objects (if TRUE)
lstQs = list();
for (i in 1:nrow(dfrParsCMsREs)){
  #--testing: i = 2;
  #                      number of rows     parameter vector       process type
  lstQs[[i]] = mkPrecMat(dfrParsCMsREs$n[i],dfrParsCMsREs$pars[i][[1]],dfrParsCMsREs$type[i],force=force);
}
Q = buildBlockDiagonalMat(lstQs,force=force);

#--calculate Cholesky factor----
Lp = chol(Q); #--UPPER TRIANGULAR Cholesky factor used for simulation

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
# for (nm in names(lstParams)){
#   lstParams[[nm]] = unname(lstParams[[nm]]);
# }

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
      lstQs[[nm]] = mkPrecMat(dfrParsCMsSMs$n[i],params[[nm]],dfrParsCMsSMs$type[i]);
    }
  }

  Q = buildBlockDiagonalMat(lstQs);

  d = determinant(as.matrix(Q));#--log determinant

  L = chol(as.matrix(Q));

  #--expected mean for REs
  zro = rep(RTMB::AD(0),length(obs));

  #--calculate negative log-likelihood
  nllFE = - dnorm(obs,mean=mu,sd=1,log=TRUE);     #  (1)
  nllRE = - dgmrf(params$r,0,Q,log=TRUE,scale=1); #  (2)
  nll = sum(nllFE) + sum(nllRE);

  #--predict observations
  prd = mu + params$r; #--NOTE: to $simulate() this, it must come after (1) and (2)

  RTMB::REPORT(params)
  RTMB::REPORT(obs)
  RTMB::REPORT(prd)
  RTMB::REPORT(lstQs)
  RTMB::REPORT(Q)
  RTMB::REPORT(d)
  RTMB::REPORT(L)
  RTMB::REPORT(nllFE)
  RTMB::REPORT(nllRE)
  RTMB::REPORT(nll)

  RTMB::ADREPORT(prd)

  return(nll);
}

#--test objfun----
objfun(lstParams);

##--create parameter `map` to handle fixed (and mirrored) values
map = list();
for (p in names(lstParams)) {
  #--testing: p = 3;
  map[[p]] = 1:(length(lstParams[[p]]));
  map[[p]] = factor(map[[p]]);
}
#--need to turn off estimation for rho in lstParams[[3]] because the associated Q is iid normal
if (FALSE){
  #--turn off estimation of rho in iid normal section
  tmp = map[[3]] |> as.numeric();
  tmp[2] = NA;
  map[[3]] = factor(tmp);
} else {
  #--turn off estimation of all variance components
  for (i in 3:4){
    tmp = map[[i]] |> as.numeric();
    tmp[] = NA;
    map[[i]] = factor(tmp);
  }
}
map = map[!(names(map) %in% c(names(lstParamsREs)))]; #--remove random effects variables

##--define limits on parameters?----
###???

#--do estimation----
require(Matrix);
require(RTMB);
testing=FALSE;
obj<-RTMB::MakeADFun(objfun,
                     parameters=lstParams,
                     random=names(lstParamsREs),
                     map=map,
                     silent=FALSE);
obj$fn();
obj$report();
opt <- nlminb(obj$par, obj$fn, obj$gr);
sdrep <- sdreport(obj);
print(sdrep);
mPrd = summary(sdrep,"report");
dfrCmp = dplyr::bind_cols(
           dfrObs,
           tibble::as_tibble(mPrd) |>
             dplyr::mutate(lwr=Estimate-1.966*`Std. Error`,
                           upr=Estimate+1.966*`Std. Error`)
         );
ggplot(dfrCmp,aes(x=i,y=obs,colour=block,group=block,fill=block)) + geom_point() + geom_line() +
  geom_line(aes(y=Estimate),linetype=2) +
  geom_ribbon(aes(y=Estimate,ymin=lwr,ymax=upr),color=NA,alpha=0.2)
