#--ALternative parameterizations for selectivity functions
#-----------------------------------------------------------------------------------
#'
#' @title Calculate a constant-valued selectivity curve
#' @description Function to calculate a constant-valued selectivity curve
#' @param z      - sizes at which to compute selectivity values
#' @param params - 1-element parameter vector
#' @param pRefZ   - reference size (dummy input)
#' @param verbose - flag (T/F) to print debugging messages
#'
#' @return named vector with selectivity values at the elements of z
#'
#' @details  The parameter vector has values
#' \itemize{
#'  \item{params[1]: the constant value}
#' }
#'
#' @export
#'
const_sel<-function(z,
                    pCnst,
                    pRefZ=0,
                    verbose=FALSE){
    if (verbose) message("Starting const_sel(...)");
    s = pCnst + 0*z;
    if (verbose) message("Finished const_sel(...)");
    return(s);
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and slope
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and slope.
#'
#'@param z     - vector of sizes at which to compute function
#'@param pZ50 - AD size at 50\% selectivity
#'@param pSlp - AD slope at `pZ50`
#'@param pRefZ   - reference size (AD, but a fixed value)
#'@param verbose - flag (T/F) to print debugging messages
#'
#' @return vector of selectivity values at the elements of `z`
#'
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'
#'If `pRefZ`>0, `pRefZ`=fully-selected size.
#'If `pRefZ`<0, function is normalized to max.
#'If `pRefZ`=0, no re-scaling is done.
#'
#'@export
#'
asclogistic<-function(z,
                      pZ50,
                      pSlp,
                      pRefZ,
                      verbose=FALSE){
  if (verbose) {
    cat("Starting asclogistic\n");
    cat(pZ50,pSlp,pRefZ,"\n")
  }
  z = as.numeric(z);
  res  <- AD(1.0)/(1.0+exp(-pSlp*(z-pZ50)));
  scl  <-AD(1);
  if (RTMB:::ad_context()){
    nrefZ = RTMB:::getValues(pRefZ);#--pRefZ is fixed, so this should not be a problem
  } else {
    nrefZ = pRefZ;
  }
  if (verbose) print(nrefZ);
  if (nrefZ>0){
      scl<-(AD(1.0)+exp(-pSlp*(pRefZ-pZ50)));
  } else if (nrefZ<0){
      scl<-AD(1.0)/res[length(res)];#--scale to max (last element)
  }
  res<-scl*res;
  if (verbose)  cat("Finished asclogistic\n");
  return(res)
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and width (1/slope)
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and width.
#'
#'@param z     - vector of sizes at which to compute selectivities
#'@param pZ50 - size at 50\% selectivity
#'@param pWdZ - width (1/slope) at `pZ50`
#'@param pRefZ   - reference size
#'@param verbose - flag (T/F) to print debugging messages
#'
#' @return named vector with selectivity values at the elements of `z`
#'
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'
#'If `pRefZ`>0, `pRefZ`=fully-selected size.
#'If `pRefZ`<0, function is normalized to max.
#'If `pRefZ`=0, no re-scaling is done.
#'
#'@export
#'
asclogistic1<-function(z,
                       pZ50,
                       pWdZ,
                       pRefZ=0,
                       verbose=FALSE){
  if (verbose) {
    cat("Starting asclogistic1\n");
    cat(pZ50,pWdZ,pRefZ,"\n")
  }
  z = as.numeric(z);
  res  <- AD(1.0)/(1.0+exp(-(z-pZ50)/pWdZ));
  scl  <-AD(1);
  if (RTMB:::ad_context()){
    nrefZ = RTMB:::getValues(pRefZ);#--pRefZ is fixed, so this should not be a problem
  } else {
    nrefZ = pRefZ;
  }
  if (verbose) print(nrefZ);
  if (nrefZ>0){
      scl<-(AD(1.0)+exp(-(pRefZ-pZ50)/pWdZ));
  } else if (nrefZ<0){
      scl<-AD(1.0)/res[length(res)];#--scale to max (last element)
  }
  res<-scl*res;
  if (verbose)  cat("Finished asclogistic1\n");
  return(res)
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and z95
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and z95.
#'
#'@param z   - vector of sizes at which to compute selectivities
#'@param pZ50 - ize at which selectivity = 0.50 (logit-scale mean)
#'@param pZ95 - size at which selectivity = 0.95
#'@param pRefZ - reference size
#'@param verbose - flag (T/F) to print debugging messages
#'
#' @return vector with selectivity values at the elements of z
#'
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'
#'If `pRefZ`>0, `pRefZ`=fully-selected size.
#'If `pRefZ`<0, function is normalized to max.
#'If `pRefZ`=0, no re-scaling is done.
#'
#'@export
#'
asclogistic5095<-function(z,
                          pZ50,
                          pZ95,
                          pRefZ=0,
                          verbose=FALSE){
  slope <- log(19.0)/(pZ95-pZ50);
  return(asclogistic(z,pZ50,slope,pRefZ,verbose));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and D95 (=z95-z50)
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and D95.
#'
#'@param z   - vector of sizes at which to compute selectivity values
#'@param pZ50 - size at which selectivity = 0.5 (logit-scale mean)
#'@param pZ9550 - z95-z50: difference between sizes at 95\% and 50\%-selected
#'@param pRefZ - reference size
#'@param verbose - flag (T/F) to print debugging messages
#'
#'@return vector with selectivity values at the elements of z
#'
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'
#'If `pRefZ`>0, `pRefZ`=fully-selected size.
#'If `pRefZ`<0, function is normalized to max.
#'If `pRefZ`=0, no re-scaling is done.
#'
#'@export
#'
asclogistic50D95<-function(z,
                           pZ50,
                           pZ9550,
                           pRefZ=0,verbose=FALSE){
    slope <- log(19.0)/pZ9550;
    return(asclogistic(z,pZ50,slope,pRefZ,verbose));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate a double logistic selectivity curve
#'@description Function to calculate a double logistic selectivity curve.
#'@param z - vector of sizes at which to compute selectivity curve
#'@param pAscZ50   - ascending limb size at which selectivity = 0.5 (logit-scale mean)}
#'@param pAscSlp - ascending limb slope at 50\%-selected}
#'@param pDscZ50   - descending limb size at which selectivity = 0.5 (logit-scale mean)}
#'@param pDscSlp - descending limb slope at 50\%-selected}
#'@param pRefZ - reference size
#'@param verbose - flag (T/F) to print debugging messages
#'
#'@return vector with selectivity values at the elements of z
#'
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'
#'If `pRefZ`>0, `pRefZ`=fully-selected size.
#'if `pRefZ`<0, function is normalized to max.
#'if `pRefZ`=0, no re-scaling is done.
#'
#'@export
#'
dbllogistic<-function(z,
                      pAscZ50,
                      pAscSlp,
                      pDscZ50,
                      pDscSlp,
                      pRefZ=0,
                      verbose=FALSE){
  if (verbose) cat('Starting dbllogistic(z,params)\n')

  res <- (1.0/(1.0+exp(-pAscSlp*(z-pAscZ50))))*(1.0/(1.0+exp(pDscSlp*(z-pDscZ50))));
  scl <- AD(1);
  if (RTMB:::ad_context()){
    nrefZ = RTMB:::getValues(pRefZ);#--pRefZ is fixed, so this should not be a problem
  } else {
    nrefZ = pRefZ;
  }
  if (nrefZ>0){
      scl<-(1.0+exp(-pAscSlp*(pRefZ-pAscZ50)))*(1.0+exp(pDscSlp*(pRefZ-pDscZ50)));
  } else if (nrefZ<0){
      scl<-1.0/max(res);
  }
  res<-scl*res;
  return(res);
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate a double logistic function parameterized by z50 and D95=z95-z50 for ascending/descending limbs
#'
#'@description Function to calculate a double logistic function parameterized by z50, z95 on ascending/descending limbs.
#'
#'@param z      - numeric vector of sizes at which to compute selectivities
#'@param pAscZ50 - ascending limb size at 50\% selected
#'@param pAscZ95 - ascending limb size at 95\% selected
#'@param pDscZ95 - descending limb size at 95\% selected
#'@param pDscZ50 - descending limb size at 50\% selected
#'@param pRefZ   - reference size
#'@param verbose - flag (T/F) to print debugging messages
#'
#'@return vector with selectivity values at the elements of z
#'
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'
#'If `pRefZ`>0, pRefZ is the fully-selected size.
#'if `pRefZ`<0, function is normalized to max.
#'if `pRefZ`=0, no re-scaling is done.
#'
#'@export
#'
dbllogistic5095<-function(z,
                          pAscZ50,
                          pAscZ95,
                          pDscZ95,
                          pDscZ50,
                          pRefZ=0,
                          verbose=FALSE){

  pAscSlp<-log(19.0)/(pAscZ95-pAscZ50);
  pDscSlp<-log(19.0)/(pDscZ50-pDscZ95);
  return(dbllogistic(z,pAscZ50,pAscSlp,pDscZ50,pDscSlp,pRefZ,verbose));
}

#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending normal selectivity curve
#'
#'@description Function to calculate an ascending normal selectivity curve.
#'
#'@param z      - vector of sizes at which to compute selectivity values
#'@param pAscZ1 - AD size at which ascending limb hits 1 (mean of a normal distribution)}
#'@param pAscWdZ - AD width of ascending limb (standard deviation of a normal distribution)}
#'@param dZ - bin size
#'@param verbose - flag (T/F) to print debugging messages
#'
#'@return vector with selectivity values at the elements of z
#'
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'This function uses a left square wave [squarewave_left()] that rises
#'from 0 to 1 "at" pAscZ1 to create a differentiable function.
#'
#'@export
#'
ascnormal1<-function(z,
                     pAscZ1,
                     pAscWdZ,
                     dZ,
                     verbose=FALSE){
  if (verbose) message(paste("Starting ascnormal1(...)"));
  ascN = exp(-0.5*((z-pAscZ1)/pAscWdZ)^2);
  sqwL = squarewave_left(pAscZ1,z,dZ);
  s    = sqwL*ascN+(1.0-sqwL);
  if (verbose) message(paste("Finished ascnormal1(...)"));
  return(s);
}
#-----------------------------------------------------------------------------------
#'
#' @title Calculate ascending normal selectivity curve
#' @description Function to calculate ascending normal selectivity curve.
#' @param z      - vector of sizes at which to compute selectivity values
#' @param pAscZ1 - size at which ascending limb reaches 1
#' @param pAscRefS - selectivity at size = pRefZ
#' @param pRefZ   - reference size (a constant) at which function reaches pAscRefS
#' @param dZ - bin size reference for join (square wave) function
#' @param verbose - flag (T/F) to print debugging messages
#'
#' @return vector with selectivity values at the elements of z
#'
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'This function uses a left square wave [squarewave_left()] that rises
#'from 0 to 1 "at" pAscZ1 to create a differentiable function.
#'
#' @export
#'
ascnormal2<-function(z,
                     pAscZ1,
                     pAscRefS,
                     pRefZ,
                     dZ,
                     verbose=FALSE){
  if (verbose) message("Starting SelFcns::ascnormal2(...)");
  ascN = exp(log(pAscRefS)*((z-pAscZ1)/(pRefZ-pAscZ1))^2);
  sqwL = squarewave_left(pAscZ1,z,dZ);
  s = sqwL*ascN+(1.0-sqwL);
  if (verbose) message("Finished ascnormal2(...)");
  return(s);
}

#-----------------------------------------------------------------------------------
#'
#' @title Calculate an ascending normal selectivity curve
#' @description Function to calculate an ascending normal selectivity curve
#' @param z      - vector of sizes at which to compute function values
#' @param pAscZ1  - size at which ascending limb reaches 1
#' @param pZatRefS - size at which selectivity = pRefS
#' @param pRefS    - reference selectivity value (constant)
#' @param dZ - bin size reference for join (square wave) function
#' @param verbose - flag (T/F) to print debugging messages
#'
#' @return vector with selectivity values at the elements of z
#'
#'@details The parameter pRefS is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'This function uses a left square wave [squarewave_left()] that rises
#'from 0 to 1 "at" pAscZ1 to create a differentiable function.
#'
#' @export
#'
ascnormal2a<-function(z,
                      pAscZ1,
                      pZatRefS,
                      pRefS,
                      dZ,verbose=FALSE){
  if (verbose) message("Starting SelFcns::ascnormal2a(...)");
  ascN = exp(log(pRefS)*((z-pAscZ1)/(pZatRefS-pAscZ1))^2);
  sqwL = squarewave_left(pAscZ1,z,dZ);
  s = sqwL*ascN + (1.0-sqwL);
  if (verbose) message("Finished SelFcns::ascnormal2a(...)");
  return(s);
}

#-----------------------------------------------------------------------------------
#' @title Calculate an ascending normal selectivity curve
#' @description Function to calculate an ascending normal selectivity curve.
#' @param pAscZ1 - size at which ascending limb reaches 1
#' @param pDZ2RefS - delta from size at 1 to size at which selectivity=pRefS
#' @param pRefS    - selectivity at pMnZ-pDZ2RefS
#'
#'@details The parameter pRefS is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'This function uses a left square wave [squarewave_left()] that rises
#'from 0 to 1 "at" pAscZ1 to create a differentiable function.
#'
#' @return - vector of selectivity values
#' @export
#'
ascnormal2b<-function(z,
                      pAscZ1,
                      pDZ2RefS,
                      pRefS,
                      dZ,
                      verbose=FALSE){
  if (verbose) message("Starting SelFcns::ascnormal2b(...)");
  ascN = exp(log(pRefS)*square((z-pAscZ1)/pDZ2RefS));
  sqwL = squarewave_left(pAscZ1,z,dZ);
  s = sqwL*ascN + (1.0-sqwL);
  if (verbose) message("Finished SelFcns::ascnormal2b(...)");
  return(s);
}

#-----------------------------------------------------------------------------------
#' @title Calculates a 2-parameter ascending normal selectivity curve
#' @description Function to calculate a 2-parameter ascending normal selectivity curve.
#' @param pDZ1 - max possible size at which the curve could reach 1
#' @param pSatZ2 - selectivity at size=pRefZ2
#' @param pMxZ1 - max possible size at which the curve could reach 1 (constant)
#' @param pRefZ2 - reference size at which curve reaches the value of pSatZ2 (constant)
#'
#'@details The parameters `pMxZ1` and `pRefZ2` are constanta and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'This function uses a left square wave [squarewave_left()] that rises
#'from 0 to 1 "at" `pMxZ1-pDZ1` to create a differentiable function.
#'
#' @return - vector of selectivity values
#' @export
#'
ascnormal3<-function(z,
                     pDZ1,
                     pSatZ2,
                     pMxZ1,
                     pRefZ2,
                     dZ,verbose=FALSE){
  if (verbose) message("Starting SelFcns::ascnormal3(...)");
  ascZ1   = pMxZ1-pDZ1;#--size at which ascending limb hits 1
  ascN = exp(log(pSatZ2)*square((z-ascZ1)/(pRefZ2-ascZ1)));
  sqwL = squarewave_left(ascZ1,z,dZ);
  s = sqwL*ascN + (1.0-sqwL);
  if (verbose) message("Finished SelFcns::ascnormal3(...)");
  return(s);
}

#-----------------------------------------------------------------------------------
#' @title Calculates a 4-parameter normal selectivity curve
#' @description Function to calculate a 4-parameter normal selectivity curve.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'
#' @param z      - dvector of sizes at which to compute function values
#' @param pAscZ1 - size at which ascending limb reaches 1
#' @param pAscWd - width of ascending limb
#' @param pDscDZ - offset to size at which descending limb departs from 1
#' @param pDscWd - width of descending limb
#' @param dZ - bin size reference for join (square wave) function
#' @param verbose - flag to print debugging info
#'
#' @return - vector of selectivity values
#' @export
#'
dblnormal4<-function(z,
                     pAscZ1,
                     pAscWd,
                     pDscDZ,
                     pDscWd,
                     dZ,
                     verbose=FALSE){
    if (verbose) message("Starting SelFcns::dblnormal4(...)");
    dscZ1 = pAscZ1+pDscDZ;
    ascN = exp(-0.5*square((z-pAscZ1)/pAscWd));
    dscN = exp(-0.5*square((z- dscZ1)/pDscWd));
    ascJ = squarewave_left(pAscZ1,z,dZ);
    dscJ = squarewave_right(dscZ1,z,dZ);
    s = elem_prod(elem_prod(ascJ,ascN)+(1.0-ascJ), elem_prod(dscJ,dscN)+(1.0-dscJ));
    if (verbose) message("Finished SelFcns::dblnormal4(...)");
    return(s);
}

#--dblnormal4a-----
#' @title Calculates a 4-parameter normal selectivity curve
#' @description Function to calculate a 4-parameter normal selectivity curve.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'
#' @param z      - vector of sizes at which to compute function values
#' @param pAscZ1 - size at which ascending limb reaches 1
#' @param pAscWd - width of ascending limb
#' @param pDscSclDZ - scaled increment to pAscZ1 at which descending limb departs from 1
#' @param pDscWd - width of descending limb
#' @param pRefZ   - max possible size (e.g., max(z))
#' @param verbose - flag to print debugging info
#'
#' @return - named vector of selectivity values
#' @export
#'
dblnormal4a<-function(z,
                      pAscZ1,
                      pAscWd,
                      pDscSclDZ,
                      pDscWd,
                      pRefZ,
                      dZ,verbose=FALSE){
    if (verbose) message("Starting SelFcns::dblnormal4a(...)");
    dscZ1 = (pRefZ-pAscZ1)*pDscSclDZ + pAscZ1;
    ascN = exp(-0.5*square((z-pAscZ1)/pAscWd));
    dscN = exp(-0.5*square((z- dscZ1)/pDscWd));
    ascJ = squarewave_left(pAscZ1,z,dZ);
    dscJ = squarewave_right(dscZ1,z,dZ);
    s = elem_prod(elem_prod(ascJ,ascN)+(1.0-ascJ), elem_prod(dscJ,dscN)+(1.0-dscJ));
    if (verbose) message("Finished SelFcns::dblnormal4a(...)");
    return(s);
}

#--dblnormal6-----
#' @title Calculates a 6-parameter normal selectivity curve
#' @description Function to calculate a 6-parameter normal selectivity curve.

#'@details Uses differentiable square wave functions to provide "joins"
#' between sills and limbs
#'
#' @param z      - dvector of sizes at which to compute function values
#' @param pAscZ1 - size at which ascending limb reaches 1
#' @param pAscWd - width of ascending limb
#' @param pDscDZ - offset to size at which descending limb departs from 1
#' @param pDscWd - width of descending limb
#' @param pAscFlr floor of ascending limb
#' @param pDscFlr floor of descending limb
#' @param dZ - bin size reference for join (square wave) functions
#' @param verbose - lag to print debugging info
#'
#' @return - vector of selectivity values
#' @export
#'
dblnormal6<-function(z,
                     pAscZ1,
                     pAscWd,
                     pDscZ1,
                     pDscWd,
                     pAscFlr,
                     pDscFlr,
                     dZ,
                     verbose=FALSE){
    if (verbose) message("Starting SelFcns::dblnormal6(...)");
    ascN = pAscFlr+(1.0-pAscFlr)*exp(-0.5*square((z-pAscZ1)/pAscWd));
    dscN = pDscFlr+(1.0-pDscFlr)*exp(-0.5*square((z-pDscZ1)/pDscWd));
    ascJ = squarewave_left( pAscZ1,z,dZ);
    dscJ = squarewave_right(pDscZ1,z,dZ);
    s = elem_prod(elem_prod(ascJ,ascN)+(1.0-ascJ), elem_prod(dscJ,dscN)+(1.0-dscJ));
    if (verbose) message("Finished SelFcns::dblnormal6(...)");
    return(s);
}

#--stackedLogistic1-----
#' @title Calculates a 5-parameter "stacked" logistic selectivity curve
#' @description Function to calculate a 5-parameter "stacked" logistic selectivity curve.
#'
#' @details Calculates 5-parameter "stacked" logistic selectivity curve.
#'
#' @param z      - vector of sizes at which to compute function values
#' @param pMnZ1 - size at inflection point for 1st logistic curve
#' @param pSdZ1 - sd for 1st logistic curve
#' @param pMnZ2 -  size at inflection point for 2nd logistic curve
#' @param pSdZ2 -sd for the 2nd logistic curve
#' @param pOmga  - weighting factor on the first curve
#' @param verbose   - flag to print debugging info
#'
#' @return - vector of selectivity values
#' @export
#'
stackedLogistic1<-function(z,
                           pMnZ1,
                           pSdZ1,
                           pMnZ2,
                           pSdZ2,
                           pOmga,
                           verbose=FALSE){
    if (verbose) message("Starting SelFcns::stackedLogistic(...)");
    s1<-asclogistic1(z,pMnZ1,pSdZ1,pRefZ=AD(0),FALSE);
    s2<-asclogistic1(z,pMnZ2,pSdZ2,pRefZ=AD(0),FALSE);
    s<-pOmga*s1 + (1-pOmga)*s2;
    if (verbose) message("Finished SelFcns::stackedLogistic1(...)");
    return(s);
}

#--stackedLogistic2 (not yet implemented in calcSelectivity)-----
#' @title Calculates a 6-parameter "stacked" logistic selectivity curve
#' @description Function to calculate a 6-parameter "stacked" logistic selectivity curve.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates 6-parameter "stacked" logistic selectivity curve parameterized by
#' \itemize{
#'      \item params[1]: asymptote for 1st logistic curve
#'      \item params[2]: size at inflection point for 1st logistic curve
#'      \item params[3]: sd for 1st logistic curve
#'      \item params[4]: asymptote for complete stacked curve
#'      \item params[5]: size at inflection point for 2nd logistic curve
#'      \item params[6]: sd for the 2nd logistic curve
#' }
#' @param z      - vector of sizes at which to compute function values
#' @param params - vector of function parameters
#' @param pRefZ   - ignored
#'
#' @return - vector of selectivity values
#' @export
#'
stackedLogistic2<-function(z,params,pRefZ=0,verbose=FALSE){
    if (verbose) message("Starting SelFcns::stackedLogistic(...)");
    asm1  = params[1];#--asymptote for 1st logistic curve
    mnZ1  = params[2];#--size at inflection point for 1st logistic curve
    sdZ1  = params[3];#--sd for 1st logistic curve
    asmT  = params[4];#--asymptote for complete stacked curve
    mnZ2  = params[5];#--size at inflection point for 2nd logistic curve
    sdZ2  = params[6];#--sd for 2nd logistic curve
    s1<-asclogistic(z,c(mnZ1,sdZ1),pRefZ=0,FALSE);
    s2<-asclogistic(z,c(mnZ2,sdZ2),pRefZ=0,FALSE);
    s<-asm1*s1 + (asmT-asm1)*s2;
    if (verbose) message("Finished SelFcns::stackedLogistic2(...)");
    return(s);
}

#--selSpline (not yet implemented in calcSelectivity)-----
#' @title Calculates a selectivity curve using a spline function
#' @description Function to calculate a selectivity curve using a spline function.
#' @details Calculates a selectivity curve using a spline function based on
#' [RTMB::splinefun-method]. The values in the `params` vector are the logit-scale
#' values at the knots. The first `n` values in the `knots` vector are the knots,
#' where `n` is th length of `params`.
#' @param z      - vector of sizes at which to compute function values
#' @param params - logit-scale vector of function parameters
#' @param knots - (AD constant) vector of knots
#'
#' @return - vector of selectivity values
#' @export
#'
selSpline<-function(z,params,knots,verbose=FALSE){
  if (verbose) message("Starting SelFcns::selSpline(...)");
  print(knots);
  if (RTMB:::ad_context()){
    xk = RTMB:::getValues(knots);
  } else {
    xk = knots;
  }
  yk = params;
  print(xk);
  print(yk);
  yk = exp(yk)/(1.0+exp(yk));
  print(yk);
  sf = RTMB::splinefun(xk,yk,"natural");#--returns a function
  str(sf)
  s = sf(z);#--values of sf evaluated at z
  if (verbose) message("Finished SelFcns::selSpline(...)");
  return(s);
}

#--selSplineClmpd (not yet implemented in calcSelectivity)-----
#' @title Calculates a selectivity curve using a clamped spline function
#' @description Function to calculate a selectivity curve using a clamped spline function.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates a selectivity curve using a clamped spline function based on
#' [RTMB::splinefun()]. The values in the `params` vector are the logit-scale
#' values at the knots. The first `n` values in the `knots` vector are the knots,
#' where `n` is th length of `params`. The spline is "clamped" at both ends
#' by adding `nr` extra knots at each end with values equal to the first (lower knots)
#' or last (upper knots) values in `param`. The separation of the extra knots is
#' taken as `1/(nr+1)` times the minimum separation in the `z` vector. Any values extrapolated
#' for `z` values less than the first extra knot (or greater than the last extra knot)
#' will be the same as the the inverse logit-transformed value of the first (last) value
#' of the `param` vector.
#' @param z      - vector of sizes at which to compute function values
#' @param params - logit-scale vector of function parameters (values at knots)
#' @param knots - (AD constant) vector of knots
#'
#' @return - vector of selectivity values at `z`
#' @importFrom RTMB splinefun
#' @export
#'
selSplineClmpd<-function(z,params,knots,verbose=FALSE){
  if (verbose) message("Starting SelFcns::selSplineClmpd(...)");
  nr = 10; #--number of values to repeat
  if (RTMB:::ad_context()){
    xkp = RTMB:::getValues(knots);
  } else {
    xkp = knots;
  }
  dk = min(diff(z))/(nr+1);
  nk = length(params);
  xk = c(xkp[1]-dk*rev(1:nr),xkp[1:nk],xk[nk]+dk*(1:nr));
  yk = c(rep(params[1],nr),params,rep(params[nk],nr));
  #--testing with stats:
  #  sf = stats::splinefun(xk,yk,"fmm");#--returns a function
  #  sf(z,2); sf(z,1);sf(z,0);
  sf = RTMB::splinefun(xk,yk,"natural");#--returns a function
  s = sf(z);#--logit-scale values of sf evaluated at z
  s = exp(s)/(1.0+exp(s));
  if (verbose) message("Finished SelFcns::selSplineClmpd(...)");
  return(s);
}

#--selSplineClmpdRight (not yet implemented in calcSelectivity)-----
#' @title Calculates a selectivity curve using a right-clamped spline function
#' @description Function to calculate a selectivity curve using a right-clamped spline function.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates a selectivity curve using a right-clamped spline function based on
#' [RTMB::splinefun()]. The values in the `params` vector are the logit-scale
#' values at the knots. The first `n` values in the `knots` vector are the knots,
#' where `n` is th length of `params`. The spline is "clamped" at the righthand end
#' by adding `nr` extra knots at the end with values equal to the
#' last (upper knots) values in `param`. The separation of the extra knots is
#' taken as `1/(nr+1)` the minimum separation in the `z` vector. Any values extrapolated
#' for `z` values greater than the last extra knot
#' will be the same as the the inverse logit-transformed value of the last value
#' of the `param` vector.
#' @param z      - vector of sizes at which to compute function values
#' @param params - logit-scale vector of function parameters (values at knots)
#' @param knots - (AD constant) vector of knots
#'
#' @return - vector of selectivity values at `z`
#' @importFrom RTMB splinefun
#' @export
#'
selSplineClmpdRight<-function(z,params,knots,verbose=FALSE){
  if (verbose) message("Starting SelFcns::selClmpSplineRight(...)");
  nr = 10; #--number of values to repeat
  if (RTMB:::ad_context()){
    xkp = RTMB:::getValues(knots);
  } else {
    xkp = knots;
  }
  dk = min(diff(z))/(nr+1);
  nk = length(params);
  xk = c(xkp[1:nk],xkp[nk]+dk*(1:nr));
  yk = c(params,rep(params[nk],nr));
  #--testing with stats:
  #  sf = stats::splinefun(xk,yk,"fmm");#--returns a function
  #  sf(z,2); sf(z,1);sf(z,0);
  sf = RTMB::splinefun(xk,yk,"natural");#--returns a function
  s = sf(z);#--logit-scale values of sf evaluated at z
  s = exp(s)/(1.0+exp(s));
  if (verbose) message("Finished SelFcns::selClmpSplineRight(...)");
  return(s);
}

#--selSplineClmpdLeft (not yet implemented in calcSelectivity)-----
#' @title Calculates a selectivity curve using a left-clamped spline function
#' @description Function to calculate a selectivity curve using a left-clamped spline function.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates a selectivity curve using a clamped spline function based on
#' [RTMB::splinefun()]. The values in the `params` vector are the logit-scale
#' values at the knots. The first `n` values in the `knots` vector are the knots,
#' where `n` is th length of `params`. The spline is "clamped" at the left end
#' by adding `nr` extra knots at each end with values equal to `param[1]`.
#' The separation of the extra knots is
#' taken as `1/(nr+1)` times the minimum separation in the `z` vector. Any values extrapolated
#' for `z` values less than the first extra knot
#' will be the same as the the inverse logit-transformed value of the first value
#' of the `param` vector.
#' @param z      - vector of sizes at which to compute function values
#' @param params - logit-scale vector of function parameters (values at knots)
#' @param knots - (AD constant) vector of knots
#'
#' @return - vector of selectivity values at `z`
#' @importFrom RTMB splinefun
#' @export
#'
selSplineClmpdLeft<-function(z,params,knots,verbose=FALSE){
  if (verbose) message("Starting SelFcns::selClmpSplineLeft(...)");
  nr = 10; #--number of values to repeat
  if (RTMB:::ad_context()){
    xkp = RTMB:::getValues(knots);
  } else {
    xkp = knots;
  }
  dk = min(diff(z))/(nr+1);
  nk = length(params);
  xk = c(xkp[1]-dk*rev(1:nr),xkp[1:nk]);
  yk = c(rep(params[1],nr),params);
  #--testing with stats:
  #  sf = stats::splinefun(xk,yk,"fmm");#--returns a function
  #  sf(z,2); sf(z,1);sf(z,0);
  sf = RTMB::splinefun(xk,yk,"natural");#--returns a function
  s = sf(z);#--logit-scale values of sf evaluated at z
  s = exp(s)/(1.0+exp(s));
  if (verbose) message("Finished SelFcns::selClmpSplineLeft(...)");
  return(s);
}

