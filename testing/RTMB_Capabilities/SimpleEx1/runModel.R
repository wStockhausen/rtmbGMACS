#--script to test RTMB concepts
require(RTMB);
source("objfun.R",chdir=TRUE);

verbose=TRUE;
data = list(xobs=c(4,5));
params = list(p=c(10,4));
res = objfcn(params);

verbose = FALSE;
obj <- RTMB::MakeADFun(objfcn,
                       params,
                       silent=FALSE);
verbose=TRUE;#<-doesn't work on next line
obj$fn(obj$par)

opt = nlminb(obj$par,obj$fn,obj$gr);

