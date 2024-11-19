#'
#' @title Calculate the objective function for a gmacs model
#' @description Function to calculate the objective function for a gmacs model.
#' @param params - named list of parameter values (p=1 for now)
#'
#' @return The objective function (negative log-likelihood) value given the data and parameters
#'
#' @details Input data and options should be organized in the calling environment as a list of lists
#' with the name "inputs". `inputs` should have the structure
#' \itemize{
#'   \item{dims - a list with model dimensions maps}
#'   \item{data - a named list with data-related entries}
#'   \item{options - a named list with selected options}
#'   \item{testing - a logical to indicate the function is run in a testing environment}
#' }
#'
#' Example
#' ```{r, eval=FALSE}
#' inputs = list(data=NULL,options=NULL,testing=TRUE);
#' params = list(dummy=0.0);
#' obj = MakeADFun(objfun,params,silent=TRUE);
#' opt = nlminb(obj$par, obj$fn, obj$gr);
#' sdreport(obj);
#' ````
#'
#' @import RTMB
#'
#' @export
#' @md
#'
obj_fun<-function(params){
  #--expand the list of parameters
  getAll(inputs,params);
  if (verbose) cat(objects(all.names=TRUE),"\n");

  #--inputs has elements
  #----dims, data, options, testing, verbose

  #--set up model processes----
  ##--initial N's----
  initN = calcInitN(dims,options,params,verbose);
  ##--Recruitment----
  ##--Natural mortality----
  ##--Growth----
  ##--Selectivity/retention----
  ##--Fishing mortality----
  ##--Survey catchability----

  #--run operating model----
  N = initN;
  for (y_ in dims$y_){
    for (s_ in dims$s_){
      if (verbose) cat("y, s =",y_,s_,"\n")
    }
  }

  #--run estimation model----

  #--calculate likelihoods

  if (testing){
    nll = -dnorm(1,dummy,1,log=TRUE);
  }
  if (verbose) cat("end objective function\n");
  return(nll);
}

#'
#' @title Calculate initial population abundance
#' @description Calculate initial population abundance.
#' @param dims - model dimensions list
#' @param options - model options list
#' @param params - model parameters list
#' @param verbose - flag to print diagnostic info
#' @return (ad)vector with initial population abundance
#' @details Options for calculating the initial N vector include:
#' \itemize{
#'  \item{noPop - all elements of the vector are zero}
#' }
#' @export
#'
calcInitN<-function(dims,options,params,verbose){
  if (verbose) cat("starting calcInitN\n");
  if (tolower(options$initN)==tolower("noPop")){
    ra = ragged_array(0.0,dims$dmsN);
    if (verbose) print(ra);
    initN = AD(as.vector(ra));
  } else {
    stop(paste0("Initial N option '",options$initN,"' not implemented."))
  }
  if (verbose) cat("finished calcInitN\n");
  return(initN);
}
