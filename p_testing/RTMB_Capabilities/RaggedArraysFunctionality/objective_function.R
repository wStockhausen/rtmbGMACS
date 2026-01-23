require(RTMB)
objective_function<-function(params){
  #--inputs
  #----data: list(y,i,j)
  #----opts: list(objfun)

  #--params
  #----beta : vector
  #----u    : vector
  #----sigma: scalar
  #----tau  : scalar

  getAll(inputs,params);

  y = data$y;
  i = data$i;
  j = data$j;

  nll = 0.0;
  nll = -sum(dnorm(u,0,tau,TRUE));

  z <- AD(0.0);
  n <- length(y);
  #t   <- AD(vector("numeric",n));
  t <- AD(y);
  eta <- AD(vector("numeric",n));
  if (testing){
    cat("class(t)   = ",class(t),"\n");
    cat("class(eta) = ",class(eta),"\n");
  }
  for (k in 1:n){
    eta[k] <- beta[i[k]] + u[j[k]];
    t[k]   <- AD(y[k]) - eta[k];
    if (opts$like_type=="normal"){
      nll <- nll - dnorm(y[k],eta[k],sigma,TRUE);
    } else {
      #nll <- nll - dnorm(t[k],z,sigma,TRUE);
      nll <- nll - dt(t[k],sigma,log=TRUE);
    }
  }
  if (testing){
    cat("class(t)   = ",class(t),"\n");
    cat("class(eta) = ",class(eta),"\n");
    cat("eta = \n"); print(eta);
    cat("t   = \n"); print(t);
  }
  REPORT(eta)
  REPORT(t)
  ADREPORT(t)
  return(nll);
}
