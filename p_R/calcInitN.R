#'
#' @title Calculate initial population abundance
#' @description Calculate initial population abundance.
#'
#' @param dims - model dimensions list
#' @param options - model options list
#' @param params - model parameters list
#' @param verbose - flag to print diagnostic info
#'
#' @return (ad)vector with initial population abundance
#'
#' @details Options for calculating the initial N vector include:
#' \itemize{
#'  \item{zeroPop - all elements of the initial N vector are zero}
#'  \item{inputPop - all elements of the initial N vector are given}
#' }
#' @export
#'
calcInitN<-function(dims,options,params,verbose){
  if (verbose) cat("starting calcInitN\n");
  if (tolower(options$initN)==tolower("zeroPop")){
    ra = ragged_array(0.0,dims$dmsN);
    if (verbose) print(ra);
    initN = AD(as.vector(ra));
  } else   if (tolower(options$initN)==tolower("inputPop")){
    ra = ragged_array(params,dims$dmsN);
    initN = AD(as.vector(ra));
    if (verbose) print(ra);
    initN = AD(as.vector(ra));
  } else {
    stop(paste0("Initial N option '",options$initN,"' not implemented."))
  }
  if (verbose) cat("finished calcInitN\n");
  return(initN);
}
