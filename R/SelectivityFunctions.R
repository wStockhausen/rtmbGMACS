#'
#'@title Calculate values for a selectivity curve
#'
#'@description Calculate values for a selectivity curve.
#'
#'@param type - the type of selectivity function to calculate
#'@param z - vector of values at which to calculate the function
#'@param params - the selectivity function parameters, as a vector
#'
#'@return vector matching size of z, with names given by elements of z
#'
#'@details "type" may be one of
#'\itemize{
#' \item{"sel_const"}
#' \item{"asclogistic"}
#' \item{'asclogistic5095'}
#' \item{'asclogistic50D95'}
#' \item{'dbllogistic'}
#' \item{'dbllogistic5095'}
#' \item{"ascnormal"}
#' \item{"ascnormal2"}
#' \item{"ascnormal2a"}
#' \item{"ascnormal3"}
#' \item{"dblnormal4"}
#' \item{"dblnormal4a"}
#' \item{"dblnormal6"}
#' }
#'
#'@export
#'
calcSelFcn<-function(type,z,params,ref=0,debug=FALSE){
    if (debug) message('sel function =',type);
    if (type=='const_sel'){
        res<-const_sel(z,params,ref,debug);
    } else if (type=='asclogistic'){
        res<-asclogistic(z,params,ref,debug);
    } else if (type=='asclogistic5095'){
        res<-asclogistic5095(z,params,ref,debug);
    } else if (type=='asclogistic50D95'){
        res<-asclogistic50D95(params,ref,debug);
    # } else if (type=='asclogistic50LnD95'){
    #     res<-asclogistic50LnD95(params,ref,debug);
    # } else if (type=='asclogisticLn50LnD95'){
    #     res<-asclogisticLn50LnD95(params,ref,debug);
    } else if (type=='dbllogistic'){
        res<-dbllogistic(z,params,ref,debug);
    } else if (type=='dbllogistic5095'){
        res<-dbllogistic5095(z,params,ref,debug);
    } else if (type=='ascnormal'){
        res<-ascnormal(z,params,ref,debug);
    } else if (type=='ascnormal2'){
        res<-ascnormal2(z,params,ref,debug);
    } else if (type=='ascnormal2a'){
        res<-ascnormal2a(z,params,ref,debug);
    } else {
        stop(paste0('Selectivity/retention function type "',type,'" not recognnized.\n',
                    'Aborting...\n'));
    }
    return(res);
}
#-----------------------------------------------------------------------------------
#'
#' @title Calculate a constant-valued selectivity curve
#' @description Function to calculate a constant-valued selectivity curve
#' @param z      - sizes at which to compute selectivity values
#' @param params - 1-element parameter vector
#' @param refZ   - reference size (dummy input)
#' @param debug - flag (T/F) to print debugging messages
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
const_sel<-function(z, params,ref=0,debug=FALSE){
    if (debug) message("Starting const_sel(...)");
    s = params[1] + 0*z;
    names(s) = as.character(z);
    if (debug) message("Finished const_sel(...)");
    return(s);
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and slope
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and slope.
#'
#'@param z     - vector of sizes at which to compute selectivities
#'@param params - 2-element vector with selectivity function parameters
#'@param refZ   - reference size
#'@param debug - flag (T/F) to print debugging messages
#'
#' @return named vector with selectivity values at the elements of z
#'
#'@details The parameter values are
#'\itemize{
#' \item{params[1] - pZ50 - size at which selectivity = 0.5 (logit-scale mean)}
#' \item{params[2] - pWdZ - width (1/slope) at z50}
#'}
#'
#'If `refZ`>0, `refZ`=fully-selected size. if `refZ`<0, function is normalized to max.
#'if `refZ`=0, no re-scaling is done.
#'
#'@export
#'
asclogistic<-function(z,params,fsz=0,debug=FALSE){
    z50   = params[1];
    slope = 1.0/params[2];
    res <- 1.0/(1.0+exp(-slope*(z-z50)));
    scl <-1;
    if (fsz>0){
        scl<-(1.0+exp(-slope*(fsz-z50)));
    } else if (fsz<0){
        scl<-1.0/max(res);
    }
    res<-scl*res;
    names(res)<-as.character(z);#--TODO: does this work with ADs?
    return(res)
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and z95
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and z95.
#'
#'@param z   - vector of sizes at which to compute selectivities
#'@param params - 2-element vector with selectivity function parameters
#'@param refZ - reference size
#'@param debug - flag (T/F) to print debugging messages
#'
#' @return named vector with selectivity values at the elements of z
#'
#'@details The parameter values are
#'\itemize{
#' \item{params[1] - z50 - size at which selectivity = 0.50 (logit-scale mean)}
#' \item{params[2] - z95 - size at which selectivity = 0.95}
#'}
#'
#'If `refZ`>0, `refZ`=fully-selected size. if `refZ`<0, function is normalized to max.
#'if `refZ`=0, no re-scaling is done.
#'
#'@export
#'
asclogistic5095<-function(z,params,refZ=0,debug=FALSE){
  z50 = params[1];
  slope = log(19.0)/(params[2]-params[1]);
  return(asclogistic(z,c(z50,slope),refZ,debug));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and D95 (=z95-z50)
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and D95.
#'
#'@param z   - vector of sizes at which to compute selectivity values
#'@param params - 2-element vector with selectivity function parameters
#'@param refZ - reference size
#'@param debug - flag (T/F) to print debugging messages
#'
#'@return named vector with selectivity values at the elements of z
#'
#'@details The parameter values are
#'\itemize{
#' \item{params[1] - z50   - size at which selectivity = 0.5 (logit-scale mean)}
#' \item{params[2] - z95-z50 - difference between sizes at 95\% and 50\%-selected}
#'}
#'
#'If `refZ`>0, `refZ`=fully-selected size. if `refZ`<0, function is normalized to max.
#'if `refZ`=0, no re-scaling is done.
#'
#'@export
#'
asclogistic50D95<-function(z,params,refZ=0,debug=FALSE){
    slope = log(19.0)/params[2];
    return(asclogistic(z,c(params[1],slope),refZ,debug));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate a logistic selectivity curve
#'@description Function to calculate a logistic selectivity curve.
#'@param z - vector of sizes at which to compute selectivity curve
#'@param params - 4-element vector of selectivity function parameters
#'@param refZ - reference size
#'@param debug - flag (T/F) to print debugging messages
#'
#'@return named vector with selectivity values at the elements of z
#'
#'@details The parameter values are
#'\itemize{
#' \item{params[1] - ascZ50   - ascending limb size at which selectivity = 0.5 (logit-scale mean)}
#' \item{params[2] - ascSlope - ascending limb slope at 50\%-selected}
#' \item{params[3] - dscZ50   - descending limb size at which selectivity = 0.5 (logit-scale mean)}
#' \item{params[4] - dscSlope - descending limb slope at 50\%-selected}
#'}
#'
#'If `refZ`>0, `refZ`=fully-selected size. if `refZ`<0, function is normalized to max.
#'if `refZ`=0, no re-scaling is done.
#'
#'@export
#'
dbllogistic<-function(z,params,fsz=0,debug=FALSE){
  #cat(z,'\n')
  #cat('z50, lnD = ',z50,lnD,'\n')
  ascZ50   = params[1];
  ascSlope = params[2];
  dscZ50   = params[3];
  dscSlope = params[4];
  res <- 1.0/(1.0+exp(-ascSlope*(z-ascZ50)))*1.0/(1.0+exp(dscSlope*(z-dscZ50)));
  scl <-1;
  if (fsz>0){
      scl<-(1.0+exp(-ascSlope*(fsz-ascZ50)))*(1.0+exp(dscSlope*(fsz-dscZ50)));
  } else if (fsz<0){
      scl<-1.0/max(res);
  }
  res<-scl*res;
  names(res)<-as.character(z);
  return(res);
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate a logistic function parameterized by z50 and D95=z95-z50 for ascending/descending limbs
#'
#'@description Function to calculate a logistic function parameterized by z50, z95 on ascending/descending limbs.
#'
#'@param z      - vector of sizes at which to compute selectivities
#'@param params - 4-element parameter vector
#'@param refZ   - reference size
#'@param debug - flag (T/F) to print debugging messages
#'
#'@return named vector with selectivity values at the elements of z
#'
#'@details The parameters are
#' \itemize{
#'  \item{params[1] - ascending limb size at 50\% selected}
#'  \item{params[2] - ascending limb size at 95\% selected}
#'  \item{params[3] - descending limb size at 95\% selected}
#'  \item{params[4] - descending limb size at 50\% selected}
#' }
#'
#'If `refZ`>0, fsz=fully-selected size. if `refZ`<0, function is normalized to max.
#'if `refZ`=0, no re-scaling is done.
#'
#'@export
#'
dbllogistic5095<-function(z,params,refZ=0,debug=FALSE){
  ascZ50 = params[1];
  ascZ95 = params[2];
  dscZ95 = params[3];
  dscZ50 = params[4];
  ascSlope<-log(19.0)/(ascZ95-ascZ50);
  dscSlope<-log(19.0)/(dscZ50-dscZ95);
  return(dbllogistic(z,c(ascZ50,ascSlope,dscZ50,dscSlope),refZ,debug));
}

#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending normal selectivity curve
#'
#'@description Function to calculate an ascending normal selectivity curve.
#'
#'@param z      - vector of sizes at which to compute selectivity values
#'@param params - 2-element vector of selectivity function parameters
#'@param refZ   - reference size (not used)
#'@param debug - flag (T/F) to print debugging messages
#'
#'@return named vector with selectivity values at the elements of z
#'
#'@details The parameters are
#' \itemize{
#' \item{params[1] - width of ascending limb (standard deviation of a normal distribution)}
#' \item{params[2] - size at which ascending limb hits 1 (mean of a normal distribution)}
#' }
#'
#'@export
#'
ascnormal<-function(z,params,refZ=0,debug=FALSE){
  if (debug) message(paste("Starting ascnormal(...)"));
  slp = 5.0;
  ascWdZ = params[1];#--width of ascending limb
  ascMnZ = params[2];#--size at which ascending limb hits 1
  ascN   = exp(-0.5*((z-ascMnZ)/ascWdZ)^2);
  ascJ   = 1.0/(1.0+exp(slp*(z-(ascMnZ))));
  s    = ascJ*ascN+(1.0-ascJ);
  names(s) = as.character(z);
  if (debug) message(paste("Finished ascnormal(...)"));
  return(s);
}
#-----------------------------------------------------------------------------------
#'
#' @title Calculate ascending normal selectivity curve
#' @description Function to calculate ascending normal selectivity curve.
#' @param z      - vector of sizes at which to compute selectivity values
#' @param params - 2-element vector of selectivity function parameters
#' @param refZ   - size at which function reaches params[2]
#' @param debug - flag (T/F) to print debugging messages
#'
#' @return named vector with selectivity values at the elements of z
#'
#' @details  The parameter vector has values
#' \itemize{
#'  \item{params[1]: selectivity at size = refZ}
#'  \item{params[2]: size at which ascending limb reaches 1}
#' }
#'
#' @export
#'
ascnormal2<-function(z, params, refZ=0,debug=FALSE){
  if (debug) message("Starting SelFcns::ascnormal2(...)");
  slp = 5.0;
  ascSref = params[1];#--selectivity at ascZref
  ascZref = refZ;      #--size at which selectivity reaches ascSref
  ascZ1   = params[2];#--size at which ascending limb hits 1
  ascN = exp(log(ascSref)*((z-ascZ1)/(ascZref-ascZ1))^2);
  ascJ = 1.0/(1.0+exp(slp*(z-(ascZ1))));
  s = ascJ*ascN+(1.0-ascJ);
  names(s) = as.character(z);
  if (debug) message("Finished ascnormal2(...)");
  return(s);
}

#-----------------------------------------------------------------------------------
#'
#' @title Calculate an ascending normal selectivity curve
#' @description Function to calculate an ascending normal selectivity curve
#' @param z      - vector of sizes at which to compute function values
#' @param params - vector of function parameters
#' @param refS    - reference selectivity value (default=0.5)
#' @param debug - flag (T/F) to print debugging messages
#'
#' @return named vector with selectivity values at the elements of z
#'
#' @details  The parameter vector has values
#' \itemize{
#'  \item{params[1]: size at which selectivity = refS}
#'  \item{params[2]: size at which ascending limb reaches 1}
#' }
#'
#' @export
#'
ascnormal2a<-function(z, params, refS=0.5, debug=FALSE){
    if (debug) message("Starting SelFcns::ascnormal2a(...)");
    slp = 5.0;
    ascSref  = refS;     #--selectivity at ascZref
    ascZref  = params[1];#--size at which selectivity reaches ascSref
    ascZ1    = params[2];#--size at which ascending limb hits 1
    ascN = exp(log(ascSref)*((z-ascZ1)/(ascZref-ascZ1))^2);
    ascJ = 1.0/(1.0+exp(slp*(z-(ascZ1))));
    s = ascJ*ascN + (1.0-ascJ);
    if (debug) message("Finished SelFcns::ascnormal2a(...)");
    names(s) = as.character(z);
    return(s);
}

#-----------------------------------------------------------------------------------
#' @title Calculate an ascending normal selectivity curve
#' @description Function to calculate an ascending normal selectivity curve.
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
ascnormal2b<-function(z,params,refS=0,debug=FALSE){
    if (debug) message("Starting SelFcns::ascnormal2b(...)");
    slp = 5.0;
    ascZ1    = params[1];#--size at which ascending limb hits 1
    ascSref  = refS;     #--selectivity at ascZref
    ascZref  = params[1]-params[2];#--size at which selectivity reaches ascSref
    ascN = exp(log(ascSref)*square((z-ascZ1)/(ascZref-ascZ1)));
    ascJ = 1.0/(1.0+exp(slp*(z-(ascZ1))));
    s = elem_prod(ascJ,ascN)+(1.0-ascJ);
    names(s) = as.character(z);
    if (debug) message("Finished SelFcns::ascnormal2b(...)");
    return(s);
}

#-----------------------------------------------------------------------------------
#' @title Calculates a 2-parameter ascending normal selectivity curve
#' @description Function to calculate a 2-parameter ascending normal selectivity curve.
#' @details Calculates 2-parameter normal function parameterized by
#' \itemize{
#'      \item params[1]: delta from max possible size (refZ[1]) at which ascending limb could reach 1
#'      \item params[2]: selectivity at size=refZ[2]
#' }:
#' `refZ` is a 2-element vector with elements
#' \itemize{
#'   \item refZ[1] - max possible size at which the curve could reach 1
#'   \item refZ[2] - reference size at which curve reaches the value of param[2]
#' }
#' @param z      - dvector of sizes at which to compute function values
#' @param params - dvar_vector of function parameters
#' @param refZ   - 2-element vector of reference sizes (see details)
#'
#' @return - named vector of selectivity values
#' @export
#'
ascnormal3<-function(z,params,refZ=0,debug=FALSE){
    if (debug) message("Starting SelFcns::ascnormal3(...)");
    slp = 5.0;
    ascZ1   = refZ[1]-params[1];#--size at which ascending limb hits 1
    ascSref = params[2];        #--selectivity at ascZref
    ascZref = refZ[2];          #--size at which selectivity reaches ascSref
    ascN = exp(log(ascSref)*square((z-ascZ1)/(ascZref-ascZ1)));
    ascJ = 1.0/(1.0+exp(slp*(z-(ascZ1))));
    # ggplot(tibble::tibble(z=z,ascN=ascN,ascJ=ascJ,mlt=elem_prod(ascJ,ascN)),aes(x=z)) +
    #   geom_line(aes(y=ascN)) + geom_line(aes(y=ascJ)) + geom_point(aes(y=mlt))
    s = elem_prod(ascJ,ascN)+(1.0-ascJ);
    names(s) = as.character(z);
    if (debug) message("Finished SelFcns::ascnormal3(...)");
    return(s);
}

#-----------------------------------------------------------------------------------
#' @title Calculates a 4-parameter normal selectivity curve
#' @description Function to calculate a 4-parameter normal selectivity curve.
#' @details Calculates 4-parameter normal function parameterized by
#' \itemize{
#'      \item params[1]: size at which ascending limb reaches 1
#'      \item params[2]: width of ascending limb
#'      \item params[3]: offset to size at which descending limb departs from 1
#'      \item params[4]: width of descending limb
#' }
#' @param z      - dvector of sizes at which to compute function values
#' @param params - dvar_vector of function parameters
#' @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double] NOTE: ignored!
#'
#' @return - named vector of selectivity values
#' @export
#'
dblnormal4<-function(z,params,refZ=0,debug=FALSE){
    if (debug) message("Starting SelFcns::dblnormal4(...)");
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
    if (debug) message("Finished SelFcns::dblnormal4(...)");
    return(s);
}

#--dblnormal4a-----
#' @title Calculates a 4-parameter normal selectivity curve
#' @description Function to calculate a 4-parameter normal selectivity curve.
#' @details Calculates 4-parameter normal function parameterized by
#' \itemize{
#'      \item params[1]: size at which ascending limb reaches 1
#'      \item params[2]: width of ascending limb
#'      \item params[3]: scaled increment to params[1] at which descending limb departs from 1
#'      \item params[4]: width of descending limb
#' }
#' @param z      - vector of sizes at which to compute function values
#' @param params - vector of function parameters
#' @param refZ   - max possible size (e.g., max(z))
#'
#' @return - named vector of selectivity values
#' @export
#'
dblnormal4a<-function(z,params,refZ,debug=FALSE){
    if (debug) message("Starting SelFcns::dblnormal4a(...)");
    slp = 5.0;
    ascMnZ = params[1];#--size at which ascending limb hits 1
    ascWdZ = params[2];#--width of ascending limb
    sclInc = params[3];#--scaled size at which descending limb departs from 1
    dscWdZ = params[4];#--width of descending limb
    dscMnZ = (refZ-ascMnZ)*sclInc + ascMnZ;
    ascN = exp(-0.5*square((z-ascMnZ)/ascWdZ));
    ascJ = 1.0/(1.0+exp(slp*(z-(ascMnZ))));
    dscN = exp(-0.5*square((z-dscMnZ)/dscWdZ));
    dscJ = 1.0/(1.0+exp(-slp*(z-(dscMnZ))));
    s = elem_prod(elem_prod(ascJ,ascN)+(1.0-ascJ), elem_prod(dscJ,dscN)+(1.0-dscJ));
    names(s) = as.character(z);
    if (debug) message("Finished SelFcns::dblnormal4a(...)");
    return(s);
}

#--dblnormal6-----
#' @title Calculates a 6-parameter normal selectivity curve
#' @description Function to calculate a 6-parameter normal selectivity curve.
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
#' @param refZ    - ignored
#'
#' @return - named vector of selectivity values
#' @export
#'
dblnormal6<-function(z,params,refZ=0,debug=FALSE){
    if (debug) message("Starting SelFcns::dblnormal6(...)");
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
    if (debug) message("Finished SelFcns::dblnormal6(...)");
    return(s);
}

