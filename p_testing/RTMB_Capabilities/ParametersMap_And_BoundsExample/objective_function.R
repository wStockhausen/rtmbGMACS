require(RTMB)
objective_function<-function(params){
  #DATA_VECTOR(y);
  #DATA_IVECTOR(i);
  #DATA_IVECTOR(j);

  #PARAMETER_VECTOR(beta);
  #PARAMETER_VECTOR(u);
  #PARAMETER(sigma);
  #PARAMETER(tau);

  getAll(inputs,params);

  nll = 0.0;

  nll = -sum(dnorm(u,0,tau,TRUE));
  for (k in 1:length(y)){
    eta = beta[i[k]] + u[j[k]];
    nll = nll - dnorm(y[k],eta,sigma,TRUE);
  }

  return(nll);
}
