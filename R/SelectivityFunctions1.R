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
const_sel<-function(z, pCnst,pRefZ=0,verbose=FALSE){
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
asclogistic<-function(z,pZ50,pSlp,pRefZ,verbose=FALSE){
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
asclogistic1<-function(z,pZ50,pWdZ,pRefZ=0,verbose=FALSE){
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
      scl<-(AD(1.0)+exp(-slp*(pRefZ-pZ50)/pWdZ));
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
asclogistic5095<-function(z,pZ50,pZ95,pRefZ=0,verbose=FALSE){
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
asclogistic50D95<-function(z,pZ50,pZ9550,pRefZ=0,verbose=FALSE){
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
  #cat(z,'\n')
  #cat('z50, lnD = ',z50,lnD,'\n')

  res <- 1.0/(1.0+exp(-pAscSlp*(z-pAscZ50)))*1.0/(1.0+exp(pDscSlp*(z-pDscZ50)));
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
  return(dbllogistic(z,pAscZ50,pAscSlp,pDscZ50,pDscSlope,pRefZ,verbose));
}

#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending normal selectivity curve
#'
#'@description Function to calculate an ascending normal selectivity curve.
#'
#'@param z      - vector of sizes at which to compute selectivity values
#'@param pAscWdZ - AD width of ascending limb (standard deviation of a normal distribution)}
#'@param pAscMnZ - AD size at which ascending limb hits 1 (mean of a normal distribution)}
#'@param dZ - bin size
#'@param verbose - flag (T/F) to print debugging messages
#'
#'@return named vector with selectivity values at the elements of z
#'
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#'This function uses a left square wave [squarewave_left()] that rises
#'from 0 to 1 "at" pAscMnZ to create a differentiable function.
#'
#'@export
#'
ascnormal1<-function(z,pAscWdZ,pAscMnZ,dZ,verbose=FALSE){
  if (verbose) message(paste("Starting ascnormal1(...)"));
  ascN = exp(-0.5*((z-pAscMnZ)/pAscWdZ)^2);
  sqwL = squarewave_left(pAscMnZ,z,dZ);
  s    = sqwL*ascN+(1.0-sqwL);
  if (verbose) message(paste("Finished ascnormal1(...)"));
  return(s);
}
#-----------------------------------------------------------------------------------
#'
#' @title Calculate ascending normal selectivity curve
#' @description Function to calculate ascending normal selectivity curve.
#' @param z      - vector of sizes at which to compute selectivity values
#' @param pAscSref - selectivity at size = pRefZ
#' @param pAscMnZ - size at which ascending limb reaches 1
#' @param pRefZ   - reference size (a constant) at which function reaches pAscSref
#' @param dZ - bin size reference for join (square wave) function
#' @param verbose - flag (T/F) to print debugging messages
#'
#' @return named vector with selectivity values at the elements of z
#'
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details  The parameter vector has values
#' \itemize{
#' }
#'
#' @export
#'
ascnormal2<-function(z, pAscSref,pAscMnZ, pRefZ,dZ,verbose=FALSE){
  if (verbose) message("Starting SelFcns::ascnormal2(...)");
  ascN = exp(log(pAscSref)*((z-ascMnZ)/(pRefZ-pAscMnZ))^2);
  sqwL = squarewave_left(pAscMnZ,z,dZ);
  s = sqwL*ascN+(1.0-sqwL);
  if (verbose) message("Finished ascnormal2(...)");
  return(s);
}

#-----------------------------------------------------------------------------------
#'
#' @title Calculate an ascending normal selectivity curve
#' @description Function to calculate an ascending normal selectivity curve
#' @param z      - vector of sizes at which to compute function values
#' @param params - vector of function parameters
#' @param refS    - reference selectivity value (default=0.5)
#' @param verbose - flag (T/F) to print debugging messages
#'
#' @return named vector with selectivity values at the elements of z
#'
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details  The parameter vector has values
#' \itemize{
#'  \item{params[1]: size at which selectivity = refS}
#'  \item{params[2]: size at which ascending limb reaches 1}
#' }
#'
#' @export
#'
ascnormal2a<-function(z, params, refS=0.5, verbose=FALSE){
    if (verbose) message("Starting SelFcns::ascnormal2a(...)");
    slp = 5.0;
    ascSref  = refS;     #--selectivity at ascZref
    ascZref  = params[1];#--size at which selectivity reaches ascSref
    ascZ1    = params[2];#--size at which ascending limb hits 1
    ascN = exp(log(ascSref)*((z-ascZ1)/(ascZref-ascZ1))^2);
    ascJ = 1.0/(1.0+exp(slp*(z-(ascZ1))));
    s = ascJ*ascN + (1.0-ascJ);
    if (verbose) message("Finished SelFcns::ascnormal2a(...)");
    names(s) = as.character(z);
    return(s);
}

#-----------------------------------------------------------------------------------
#' @title Calculate an ascending normal selectivity curve
#' @description Function to calculate an ascending normal selectivity curve.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Function is parameterized by
#' \itemize{
#'      \item params[1]: size at which ascending limb reaches 1
#'      \item params[2]: delta from size at 1 to size at which selectivity=refS
#' }
#' @param z      - dvector of sizes at which to compute function values
#' @param params - dvar_vector of function parameters
#' @param refS    - selectivity at params[1]-params[2]
#'
#' @return - named vector of selectivity values
#' @export
#'
ascnormal2b<-function(z,params,refS=0,verbose=FALSE){
    if (verbose) message("Starting SelFcns::ascnormal2b(...)");
    slp = 5.0;
    ascZ1    = params[1];#--size at which ascending limb hits 1
    ascSref  = refS;     #--selectivity at ascZref
    ascZref  = params[1]-params[2];#--size at which selectivity reaches ascSref
    ascN = exp(log(ascSref)*square((z-ascZ1)/(ascZref-ascZ1)));
    ascJ = 1.0/(1.0+exp(slp*(z-(ascZ1))));
    s = elem_prod(ascJ,ascN)+(1.0-ascJ);
    names(s) = as.character(z);
    if (verbose) message("Finished SelFcns::ascnormal2b(...)");
    return(s);
}

#-----------------------------------------------------------------------------------
#' @title Calculates a 2-parameter ascending normal selectivity curve
#' @description Function to calculate a 2-parameter ascending normal selectivity curve.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates 2-parameter normal function parameterized by
#' \itemize{
#'      \item params[1]: delta from max possible size (pRefZ[1]) at which ascending limb could reach 1
#'      \item params[2]: selectivity at size=pRefZ[2]
#' }:
#' `pRefZ` is a 2-element vector with elements
#' \itemize{
#'   \item pRefZ[1] - max possible size at which the curve could reach 1
#'   \item pRefZ[2] - reference size at which curve reaches the value of param[2]
#' }
#' @param z      - dvector of sizes at which to compute function values
#' @param params - dvar_vector of function parameters
#' @param pRefZ   - 2-element vector of reference sizes (see details)
#'
#' @return - named vector of selectivity values
#' @export
#'
ascnormal3<-function(z,params,pRefZ=0,verbose=FALSE){
    if (verbose) message("Starting SelFcns::ascnormal3(...)");
    slp = 5.0;
    ascZ1   = pRefZ[1]-params[1];#--size at which ascending limb hits 1
    ascSref = params[2];        #--selectivity at ascZref
    ascZref = pRefZ[2];          #--size at which selectivity reaches ascSref
    ascN = exp(log(ascSref)*square((z-ascZ1)/(ascZref-ascZ1)));
    ascJ = 1.0/(1.0+exp(slp*(z-(ascZ1))));
    # ggplot(tibble::tibble(z=z,ascN=ascN,ascJ=ascJ,mlt=elem_prod(ascJ,ascN)),aes(x=z)) +
    #   geom_line(aes(y=ascN)) + geom_line(aes(y=ascJ)) + geom_point(aes(y=mlt))
    s = elem_prod(ascJ,ascN)+(1.0-ascJ);
    names(s) = as.character(z);
    if (verbose) message("Finished SelFcns::ascnormal3(...)");
    return(s);
}

#-----------------------------------------------------------------------------------
#' @title Calculates a 4-parameter normal selectivity curve
#' @description Function to calculate a 4-parameter normal selectivity curve.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates 4-parameter normal function parameterized by
#' \itemize{
#'      \item params[1]: size at which ascending limb reaches 1
#'      \item params[2]: width of ascending limb
#'      \item params[3]: offset to size at which descending limb departs from 1
#'      \item params[4]: width of descending limb
#' }
#' @param z      - dvector of sizes at which to compute function values
#' @param params - dvar_vector of function parameters
#' @param pRefZ   - size at which function = 1 (i.e., fully-selected size) [double] NOTE: ignored!
#'
#' @return - named vector of selectivity values
#' @export
#'
dblnormal4<-function(z,params,pRefZ=0,verbose=FALSE){
    if (verbose) message("Starting SelFcns::dblnormal4(...)");
    slp = 5.0;
    ascMnZ = params[1];#--size at which ascending limb hits 1
    ascWdZ = params[2];#--width of ascending limb
    dscMnZ = params[1]+params[3];#--size at which descending limb departs from 1
    dscWdZ = params[4];#--width of descending limb
    ascN = exp(-0.5*square((z-ascMnZ)/ascWdZ));
    ascJ = 1.0/(1.0+exp(slp*(z-(ascMnZ))));
    dscN = exp(-0.5*square((z-dscMnZ)/dscWdZ));
    dscJ = 1.0/(1.0+exp(-slp*(z-(dscMnZ))));
    s = elem_prod(elem_prod(ascJ,ascN)+(1.0-ascJ), elem_prod(dscJ,dscN)+(1.0-dscJ));
    names(s) = as.character(z);
    if (verbose) message("Finished SelFcns::dblnormal4(...)");
    return(s);
}

#--dblnormal4a-----
#' @title Calculates a 4-parameter normal selectivity curve
#' @description Function to calculate a 4-parameter normal selectivity curve.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates 4-parameter normal function parameterized by
#' \itemize{
#'      \item params[1]: size at which ascending limb reaches 1
#'      \item params[2]: width of ascending limb
#'      \item params[3]: scaled increment to params[1] at which descending limb departs from 1
#'      \item params[4]: width of descending limb
#' }
#' @param z      - vector of sizes at which to compute function values
#' @param params - vector of function parameters
#' @param pRefZ   - max possible size (e.g., max(z))
#'
#' @return - named vector of selectivity values
#' @export
#'
dblnormal4a<-function(z,params,pRefZ,verbose=FALSE){
    if (verbose) message("Starting SelFcns::dblnormal4a(...)");
    slp = 5.0;
    ascMnZ = params[1];#--size at which ascending limb hits 1
    ascWdZ = params[2];#--width of ascending limb
    sclInc = params[3];#--scaled size at which descending limb departs from 1
    dscWdZ = params[4];#--width of descending limb
    dscMnZ = (pRefZ-ascMnZ)*sclInc + ascMnZ;
    ascN = exp(-0.5*square((z-ascMnZ)/ascWdZ));
    ascJ = 1.0/(1.0+exp(slp*(z-(ascMnZ))));
    dscN = exp(-0.5*square((z-dscMnZ)/dscWdZ));
    dscJ = 1.0/(1.0+exp(-slp*(z-(dscMnZ))));
    s = elem_prod(elem_prod(ascJ,ascN)+(1.0-ascJ), elem_prod(dscJ,dscN)+(1.0-dscJ));
    names(s) = as.character(z);
    if (verbose) message("Finished SelFcns::dblnormal4a(...)");
    return(s);
}

#--dblnormal6-----
#' @title Calculates a 6-parameter normal selectivity curve
#' @description Function to calculate a 6-parameter normal selectivity curve.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates 6-parameter normal function parameterized by
#' \itemize{
#'      \item params[1]: size at which ascending limb reaches 1
#'      \item params[2]: width of ascending limb
#'      \item params[3]: size at which descending limb departs from 1
#'      \item params[4]: width of descending limb
#'      \item params[5]: floor of ascending limb
#'      \item params[6]: floor of descending limb
#' }
#' @param z      - vector of sizes at which to compute function values
#' @param params - vector of function parameters
#' @param pRefZ    - ignored
#'
#' @return - named vector of selectivity values
#' @export
#'
dblnormal6<-function(z,params,pRefZ=0,verbose=FALSE){
    if (verbose) message("Starting SelFcns::dblnormal6(...)");
    slp = 5.0;
    ascMnZ = params[1];#--size at which ascending limb hits 1
    ascWdZ = params[2];#--width of ascending limb
    dscMnZ = params[3];#--size at which descending limb departs from 1
    dscWdZ = params[4];#--width of descending limb
    ascFlr = params[5];#--floor of ascending limb
    dscFlr = params[6];#--floor of descending limb
    ascN = ascFlr+(1.0-ascFlr)*exp(-0.5*square((z-ascMnZ)/ascWdZ));
    ascJ = 1.0/(1.0+exp(slp*(z-(ascMnZ))));
    dscN = dscFlr+(1.0-dscFlr)*exp(-0.5*square((z-dscMnZ)/dscWdZ));
    dscJ = 1.0/(1.0+exp(-slp*(z-(dscMnZ))));
    s = elem_prod(elem_prod(ascJ,ascN)+(1.0-ascJ), elem_prod(dscJ,dscN)+(1.0-dscJ));
    names(s) = as.character(z);
    if (verbose) message("Finished SelFcns::dblnormal6(...)");
    return(s);
}

#--stackedLogistic1-----
#' @title Calculates a 5-parameter "stacked" logistic selectivity curve
#' @description Function to calculate a 5-parameter "stacked" logistic selectivity curve.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates 5-parameter "stacked" logistic selectivity curve parameterized by
#' \itemize{
#'      \item params[1]: weighting factor on the first curve
#'      \item params[2]: size at inflection point for 1st logistic curve
#'      \item params[3]: sd for 1st logistic curve
#'      \item params[4]: size at inflection point for 2nd logistic curve
#'      \item params[5]: sd for the 2nd logistic curve
#' }
#' @param z      - vector of sizes at which to compute function values
#' @param params - vector of function parameters
#' @param pRefZ   - ignored
#'
#' @return - named vector of selectivity values
#' @export
#'
stackedLogistic1<-function(z,params,pRefZ=0,verbose=FALSE){
    if (verbose) message("Starting SelFcns::stackedLogistic(...)");
    omega = params[1];#--weighting factor on the first curve
    mnZ1  = params[2];#--size at inflection point for 1st logistic curve
    sdZ1  = params[3];#--sd for 1st logistic curve
    mnZ2  = params[4];#--size at inflection point for 2nd logistic curve
    sdZ2  = params[5];#--sd for 2nd logistic curve
    s1<-asclogistic1(z,c(mnZ1,sdZ1),pRefZ=0,FALSE);
    s2<-asclogistic1(z,c(mnZ2,sdZ2),pRefZ=0,FALSE);
    s<-omega*s1 + (1-omega)*s2;
    names(s) = as.character(z);
    if (verbose) message("Finished SelFcns::stackedLogistic1(...)");
    return(s);
}

#--stackedLogistic2-----
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
#' @return - named vector of selectivity values
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
    names(s) = as.character(z);
    if (verbose) message("Finished SelFcns::stackedLogistic2(...)");
    return(s);
}

#--selSpline-----
#' @title Calculates a selectivity curve using a spline function
#' @description Function to calculate a selectivity curve using a spline function.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates a selectivity curve using a spline function based on
#' [RTMB::splinefun-method]. The values in the `params` vector are the logit-scale
#' values at the knots. The first `n` values in the `consts` vector are the knots,
#' where `n` is th length of `params`.
#' @param z      - vector of sizes at which to compute function values
#' @param params - logit-scale vector of function parameters
#' @param consts - vector of knots
#'
#' @return - named vector of selectivity values
#' @export
#'
selSpline<-function(z,params,consts,verbose=FALSE){
  if (verbose) message("Starting SelFcns::selSpline(...)");
  nk = length(params);
  xk = consts[1:nk];
  yk = params;
  yk = exp(yk)/(1.0+exp(yk));
  sf = RTMB::splinefun(xk,yk,"natural");#--returns a function
  s = sf(z);#--values of sf evaluated at z
  names(s) = as.character(z);
  if (verbose) message("Finished SelFcns::selSpline(...)");
  return(s);
}

#--selSplineClmpd-----
#' @title Calculates a selectivity curve using a clamped spline function
#' @description Function to calculate a selectivity curve using a clamped spline function.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates a selectivity curve using a clamped spline function based on
#' [RTMB::splinefun()]. The values in the `params` vector are the logit-scale
#' values at the knots. The first `n` values in the `consts` vector are the knots,
#' where `n` is th length of `params`. The spline is "clamped" at both ends
#' by adding `nr` extra knots at each end with values equal to the first (lower knots)
#' or last (upper knots) values in `param`. The separation of the extra knots is
#' taken as `1/(nr+1)` times the minimum separation in the `z` vector. Any values extrapolated
#' for `z` values less than the first extra knot (or greater than the last extra knot)
#' will be the same as the the inverse logit-transformed value of the first (last) value
#' of the `param` vector.
#' @param z      - vector of sizes at which to compute function values
#' @param params - logit-scale vector of function parameters (values at knots)
#' @param consts - vector of knots
#'
#' @return - named vector of selectivity values at `z`
#' @importFrom RTMB splinefun
#' @export
#'
selSplineClmpd<-function(z,params,consts,verbose=FALSE){
  if (verbose) message("Starting SelFcns::selClmpSpline(...)");
  nr = 10; #--number of values to repeat
  dk = min(diff(z))/(nr+1);
  nk = length(params);
  xk = c(consts[1]-dk*rev(1:nr),consts[1:nk],consts[nk]+dk*(1:nr));
  yk = c(rep(params[1],nr),params,rep(params[nk],nr));
  #--testing with stats:
  #  sf = stats::splinefun(xk,yk,"fmm");#--returns a function
  #  sf(z,2); sf(z,1);sf(z,0);
  sf = RTMB::splinefun(xk,yk,"natural");#--returns a function
  s = sf(z);#--logit-scale values of sf evaluated at z
  s = exp(s)/(1.0+exp(s));
  names(s) = as.character(z);
  if (verbose) message("Finished SelFcns::selClmpSpline(...)");
  return(s);
}

#--selSplineClmpdRight-----
#' @title Calculates a selectivity curve using a right-clamped spline function
#' @description Function to calculate a selectivity curve using a right-clamped spline function.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates a selectivity curve using a right-clamped spline function based on
#' [RTMB::splinefun()]. The values in the `params` vector are the logit-scale
#' values at the knots. The first `n` values in the `consts` vector are the knots,
#' where `n` is th length of `params`. The spline is "clamped" at the righthand end
#' by adding `nr` extra knots at the end with values equal to the
#' last (upper knots) values in `param`. The separation of the extra knots is
#' taken as `1/(nr+1)` the minimum separation in the `z` vector. Any values extrapolated
#' for `z` values greater than the last extra knot
#' will be the same as the the inverse logit-transformed value of the last value
#' of the `param` vector.
#' @param z      - vector of sizes at which to compute function values
#' @param params - logit-scale vector of function parameters (values at knots)
#' @param consts - vector of knots
#'
#' @return - named vector of selectivity values at `z`
#' @importFrom RTMB splinefun
#' @export
#'
selSplineClmpdRight<-function(z,params,consts,verbose=FALSE){
  if (verbose) message("Starting SelFcns::selClmpSplineRight(...)");
  nr = 10; #--number of values to repeat
  dk = min(diff(z))/(nr+1);
  nk = length(params);
  xk = c(consts[1:nk],consts[nk]+dk*(1:nr));
  yk = c(params,rep(params[nk],nr));
  #--testing with stats:
  #  sf = stats::splinefun(xk,yk,"fmm");#--returns a function
  #  sf(z,2); sf(z,1);sf(z,0);
  sf = RTMB::splinefun(xk,yk,"natural");#--returns a function
  s = sf(z);#--logit-scale values of sf evaluated at z
  s = exp(s)/(1.0+exp(s));
  names(s) = as.character(z);
  if (verbose) message("Finished SelFcns::selClmpSplineRight(...)");
  return(s);
}

#--selSplineClmpdLeft-----
#' @title Calculates a selectivity curve using a left-clamped spline function
#' @description Function to calculate a selectivity curve using a left-clamped spline function.
#'@details The parameter pRefZ is a constant and must be identified as such in the
#'`map` list when MakeADFun'ing an objective function that uses this function.
#' @details Calculates a selectivity curve using a clamped spline function based on
#' [RTMB::splinefun()]. The values in the `params` vector are the logit-scale
#' values at the knots. The first `n` values in the `consts` vector are the knots,
#' where `n` is th length of `params`. The spline is "clamped" at the left end
#' by adding `nr` extra knots at each end with values equal to `param[1]`.
#' The separation of the extra knots is
#' taken as `1/(nr+1)` times the minimum separation in the `z` vector. Any values extrapolated
#' for `z` values less than the first extra knot
#' will be the same as the the inverse logit-transformed value of the first value
#' of the `param` vector.
#' @param z      - vector of sizes at which to compute function values
#' @param params - logit-scale vector of function parameters (values at knots)
#' @param consts - vector of knots
#'
#' @return - named vector of selectivity values at `z`
#' @importFrom RTMB splinefun
#' @export
#'
selSplineClmpdLeft<-function(z,params,consts,verbose=FALSE){
  if (verbose) message("Starting SelFcns::selClmpSplineLeft(...)");
  nr = 10; #--number of values to repeat
  dk = min(diff(z))/(nr+1);
  nk = length(params);
  xk = c(consts[1]-dk*rev(1:nr),consts[1:nk]);
  yk = c(rep(params[1],nr),params);
  #--testing with stats:
  #  sf = stats::splinefun(xk,yk,"fmm");#--returns a function
  #  sf(z,2); sf(z,1);sf(z,0);
  sf = RTMB::splinefun(xk,yk,"natural");#--returns a function
  s = sf(z);#--logit-scale values of sf evaluated at z
  s = exp(s)/(1.0+exp(s));
  names(s) = as.character(z);
  if (verbose) message("Finished SelFcns::selClmpSplineLeft(...)");
  return(s);
}

