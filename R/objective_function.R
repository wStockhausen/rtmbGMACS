#'
#' @title Calculate the objective function
#' @description
#' Function to calculate the objective function
#' @param params - named list of parameter values (p=1 for now)
#'
#' @return The objective function (negative log-likelihood) value given the data and parameters
#'
#' @details Input data should be organized in the calling environment as a list with
#' the name "data"
#'
#' @import RTMB
#'
#' @export
#'
objfcn<-function(params){
  #--expand the list of parameters
  getAll(params);

  obj_fun = -dnorm(p,1,1,log=TRUE);
  return(obj_fun);
}
