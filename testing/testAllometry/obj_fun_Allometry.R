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
#' inputs = list(dims=NULL,data=NULL,options=NULL,testing=TRUE);
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
  ##--Allometry----
  wAtZ = calcWatZ(dims,options,params$pWAtZ,verbose);
  ##--initial N's----
  # initN = calcInitN(dims,options,params,verbose);
  # REPORT(initN);
  # ADREPORT(initN);

  ##--Recruitment----
  # lstRec = calcRecruitment(dims,options,params,verbose);

  ##--Natural mortality----

  ##--Growth----
  ###--will calculate all growth transition matrices
  ###--includes molt probabilities, growth matrices,
  ###--maturity transition matrices, probability of terminal molt

  ##--Selectivity/Retention----
  ###--will calculate all selectivity curves
  #lstSelFcns = calcSelFcns(dims,options$sel_fcns,sel_params,verbose);

  ##--Fishery characteristics----
  ###--will calculate capture, handling mortality, and retention rates
  ###--for all fisheries
  #lstFRs = calcFshRates(dims,options,params,lstSelFcns,verbose);

  ##--Survey characteristics----
  ###--will calculate selectivity, q for all surveys

  ##--Movement----
  ###--to be implemented at some point in future

  #--run operating model----
  N = initN;
  for (y_ in dims$y_){ #--loop over years
    for (s_ in dims$s_){ #--loop over seasons win years
      if (verbose) cat("y, s =",y_,s_,"\n");
      #--instantaneous processes happen at beginning of season in order
      ##--growth
      ##--surveys
      ##--fishing mortality
      totFM = AD(vector("numeric",length(N)));
      for (f in 1:nFsh){
        totFCR = totFM + lstFMs[[f]];
      }
      N = totFM %*% N;
      ##--movement
      N = M %*% N;

      #--continuous processes happen throughout season
      ##--natural mortality rate
      diag_NMRs = getNaturalMortalityRates(y_,s_,lstNMRs);
      ##--fishery catch rates
      lst_FCRs = getFisheryCatchRates(y_,s_);#--list with diagonal matrices for fishery catch rates
      #--update N for continuous processes
      diag_MR = AD(vector("numeric",length(N)));
      N = exp(-diag_MR*dt[s_]) %*% N;
    } #--end seasons loop (s_)
  } #--end years loop (y_)

  #--run estimation model----

  #--calculate likelihoods

  if (testing){
    nll = -dnorm(1,dummy,1,log=TRUE);
  }
  if (verbose) cat("end objective function\n");
  return(nll);
}

#'
#' @title Calculate weight-at-size
#' @description Calculate weight-at-size.
#' @param dims - model dimensions dataframe
#' @param options - model options list
#' @param params - model parameters list
#' @param verbose - flag to print diagnostic info
#' @return (ad)vector with weight-at-size
#' @details Weight-at-size is parameterized as\n
#'   $$w(z) = a \cdot z^b$$
#' @md
#' @export
#'
calcWatZ<-function(dims,options,params,verbose){
  if (verbose) cat("starting calcWatZ\n");
  if (tolower(options$allometry$option)=="data"){
  } else {
    stop(paste0("Allometry option '",options$allometry$option,"' not implemented."))
  }
  if (verbose) cat("finished calcWatZ\n");
  return(initN);
}
