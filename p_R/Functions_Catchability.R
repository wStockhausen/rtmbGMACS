#--(survey) catchability functions----

#'
#' @title Calculate catchability `Q` as a function of size
#' @description Function to calculate catchability `Q` as a function of size.
#' @param z - vector of sizes at which to calculate `Q`
#' @param pLnQ - vector of ln-scale values of `Q`
#' @param verbose - flag to print diagnostic info
#' @return arithmetic-scale vector of catchabilities the same size as `z`
#'
#' @details This function simply exponentiates `pLnQ` and returns an arithmetic-scale vector.
#' `pLnQ` can be a vector, but it must be conformable to `z`.
#'
lnQ<-function(z,pLnQ,verbose=FALSE){
  if (verbose) cat("in function lnQ.\n")
  return(RTMB::AD(array(exp(pLnQ),length(z))));
}
