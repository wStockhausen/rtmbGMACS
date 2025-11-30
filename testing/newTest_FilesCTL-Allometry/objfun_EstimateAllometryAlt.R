  #--define allometric functions----
  pwrlaw1<-function(z,pA,pB){
    w = pA * (RTMB::AD(z)^pB);
    return(w);
  }

  pwrlaw2<-function(z,pLnA,pLnS,pB,pZ0){
    w = exp(pLnA+pLnS+pB*log(RTMB::AD(z)/pZ0));
    return(w);
  }

#'
#' @title Calculate the objective function for a gmacs model
#' @description Function to calculate the objective function for a gmacs model.
#' @param params - named list of parameter values
#'
#' @return The objective function (negative log-likelihood) value given the data and parameters
#'
#' @details TODO!!
#' @examplesIf FALSE
#' # example code: TODO
#'
#'
#' @import RTMB
#'
#' @export
#' @md
#'
objfun_EstimateAllometryAlt<-function(params){
  if (testing) {
    cat("starting objfun_EstimateAllometryAlt.\n");
    cat("model parameters are:\n", paste0("\t",names(params),"\n"));
  }

  if (testing) {
    cat("RTMB params: \n")
    for (nm in names(params)){
      cat("\t",nm,": ",paste(params[[nm]],colllapse=" "),"\n",sep="");
    }
  }

  #--predict individual weight-at-size based on model parameter values
  nll = RTMB::AD(0.0);
  for (r in 1:nrow(dfrZWppp)){
    #testing: r = 1;
    rw = dfrZWppp[r,];
    if (TYPE=="pwrLaw1"){
      pA = RTMB::AD(0.0); pB = RTMB::AD(0.0);
      if (rw$x=="male") {pA = params$pA[1]; pB = params$pB[1];}
      if (rw$x=="female") {
        if (rw$m=="imm") {pA = params$pA[2]; pB = params$pB[2];}
        if (rw$m=="mat") {pA = params$pA[3]; pB = params$pB[3];}
      }
      prd = pwrlaw1(RTMB::AD(as.numeric(rw$z)),pA=pA,pB=pB);
      dfrZWppp$pars[r] = list(list(pA=pA,pB=pB));#--no conversion necessary
      dfrZWppp$prd[r]  = prd;                    #--no conversion necessary
      if (rw$obs>0) {
        nllp = -dnorm(log(as.numeric(rw$obs)),log(prd),RTMB::AD(0.1),log=TRUE);
        dfrZWppp$nll[r] = nllp;
        nll = nll + nllp;
      }
      # if (testing){
      #   #--keep as numeric values
      #   dfrZWppp$pars[r] = list(list(pA=pA,pB=pB));#--no conversion necessary
      #   dfrZWppp$prd[r]  = prd;                    #--no conversion necessary
      # } else {
      #   #--convert to numeric values
      #   dfrZWppp$pars[r] = list(list(pA=RTMB:::getValues(pA),pB=RTMB:::getValues(pB)));
      #   dfrZWppp$prd[r]  = RTMB:::getValues(prd);
      # }
      # dfrZWppp$pars[r] = list(list(pA=pA,pB=pB));#--no conversion necessary
      # dfrZWppp$prd[r]  = prd;                    #--no conversion necessary
    }#--if: TYPE=="pwrLaw1"

    if (TYPE=="pwrLaw2"){
      pLnA = RTMB::AD(0.0); pB = RTMB::AD(0.0);
      if (rw$x=="male") {pLnA = params$pLnA[1]; pB = params$pB[1];}
      if (rw$x=="female") {
        if (rw$m=="imm") {pLnA = params$pLnA[2]; pB = params$pB[2];}
        if (rw$m=="mat") {pLnA = params$pLnA[3]; pB = params$pB[3];}
      }
      prd = pwrlaw2(RTMB::AD(as.numeric(rw$z)),pLnA,pLnS,pB,pZ0);
      dfrZWppp$pars[r] = list(list(pLnA=pLnA,pB=pB));#--no conversion necessary
      dfrZWppp$prd[r]  = prd;                    #--no conversion necessary
      if (rw$obs>0) {
        if (NLLTYPE=="dnorm")
          nllp = -dnorm(log(as.numeric(rw$obs)),log(prd),RTMB::AD(0.1),log=TRUE);
        if (NLLTYPE=="dlnorm")
          nllp = -dlnorm(as.numeric(rw$obs),log(prd),RTMB::AD(0.1),log=TRUE);
        dfrZWppp$nll[r] = nllp;
        nll = nll + nllp;
      }
    }#--if: TYPE=="pwrLaw2"
  }#--loop: r

  dfrZWs = dfrZWppp |> dplyr::select(obs_id,y,x,m,p,z,obs,prd,nll) |>
             dplyr::mutate(z=as.numeric(z),obs=as.numeric(obs));
  REPORT(dfrZWs)

  if (TYPE=="pwrLaw2"){
    pA = exp(params$pLnA+pLnS-params$pB*log(pZ0));
    ADREPORT(pA);
  }

  if (testing) {
    cat("end objective function\n");
    return(list(nll=nll,dfrZWs=dfrZWs))
  }
  return(nll);
}


