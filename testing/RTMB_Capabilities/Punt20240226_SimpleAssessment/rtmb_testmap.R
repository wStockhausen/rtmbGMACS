require(RTMB);
set.seed((1))
objfun<-function(params){
  getAll(params);
  y %~% dnorm(a[1]+a[2]*x+a[3]*(x-5)^2,1);
  ADREPORT(a);
}

a = c(1,2,0.3);
x = 0:10
y = a[1] + a[2]*x + a[3]*(x-5)^2 + rnorm(length(x),0,1);
plot(x,y);
params = list(a=a);
ap = 1:length(a); ap[2] = NA;
map = list(a=factor(ap));
# map$a[2] = NA;
# map = list(a=as.factor(c(1,NA)))
# params = list(a1=a[1],a2=a[2])
# map = 
  
obj = MakeADFun(objfun,parameters=params,map=map);
opt = nlminb(obj$par,obj$fn,obj$gr);
sdr = sdreport(obj);
summary(sdr)
dfr2 = getDFR_sdreport(obj,report=TRUE);

#rm(obj,opt,sdr);

