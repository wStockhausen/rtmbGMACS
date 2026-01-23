#'
#' @title The objective function
#' @description The RTMB objective function for GMACS
#' @param params - list of input parameter values
#' @return the objective function value
#'
#' @details The function expects that `data` is a list with elements that describe
#' the model configuration and data
#' @export
#'
source("./testing.R");
objfcn<-function(params){
  getAll(data,params);
  if (verbose) cat("xobs =",xobs,"\np =",p,"\n");
  x = sqr(p);
  if (verbose) cat("x =",x,"\nxobs =",xobs,"\n");
  nll = sum(0.5*((xobs-x)^2+2*pi));
  if (verbose) cat("nll = ",nll,"\n");
  objfun = -sum(dnorm(xobs,x,1,log=TRUE));
  if (verbose) cat("objfun =",objfun,"\n");
  return(objfun);
}
