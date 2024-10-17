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
#' \item{"asclogistic"}
#' \item{'asclogistic5095'}
#' \item{'asclogistic50D95'}
#' \item{'asclogistic50LnD95'}
#' \item{'asclogisticLn50LnD95'}
#' \item{'dbllogistic'}
#' \item{'dbllogistic5095'}
#' }
#'
#'@export
#'
calcSelectivity<-function(type,z,params){
    fsz <- 0; #fully-selected size for re-scaling selectivity function
    if (type=='asclogistic'){
        cat('sel function=asclogistic\n')
        if (length(params)>2) fsz<-params[3];
        res<-asclogistic(z,params[1],params[2],fsz);
    } else if (type=='asclogistic5095'){
        cat('sel function =',type,'\n')
        if (length(params)>2) fsz<-params[3];
        res<-asclogistic5095(z,params[1],params[2],fsz);
    } else if (type=='asclogistic50D95'){
        cat('sel function =',type,'\n')
        if (length(params)>2) fsz<-params[3];
        res<-asclogistic50D95(z,params[1],params[2],fsz);
    } else if (type=='asclogistic50LnD95'){
        cat('sel function =',type,'\n')
        if (length(params)>2) fsz<-params[3];
        res<-asclogistic50LnD95(z,params[1],params[2],fsz);
    } else if (type=='asclogisticLn50LnD95'){
        cat('sel function =',type,'\n')
        if (length(params)>2) fsz<-params[3];
        res<-asclogisticLn50LnD95(z,params[1],params[2],fsz);
    } else if (type=='dbllogistic'){
        cat('sel function =',type,'\n')
        if (length(params)>4) fsz<-params[5];
        res<-dbllogistic(z,params[1],params[2],params[3],params[4],fsz);
    } else if (type=='dbllogistic5095'){
        cat('sel function =',type,'\n')
        if (length(params)>4) fsz<-params[5];
        res<-dbllogistic5095(z,params[1],params[2],params[3],params[4],fsz);
    } else {
        cat('Selectivity/retention function type "',type,'" not recognnized.\n',sep='');
        cat('Aborting...\n');
        stop();
    }
    return(res)
}
#-----------------------------------------------------------------------------------
#'
#' @title Calculate a constant-valued selectivity curve
#' @description Function to calculate a constant-valued selectivity curve
#' @param z      - sizes at which to compute selectivity values
#' @param params - parameter vector
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
const_sel<-function(z, params,debug=FALSE){
    if (debug) message("Starting const_sel(...)");
    s = 0*z + params[1];
    if (debug) message("Finished const_sel(...)");
    return(s);
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate a logistic selectivity curve
#'
#'@description Function to calculate a logistic selectivity curve.
#'
#'@param z   - vector of size class midpoints at which to compute the selectivity curve
#'@param z50 - size at which selectivity = 0.5 (logit-scale mean)
#'@param sd  - standard deviation in selectivity (logit-scale standard deviation)
#'@param fsz - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'@param debug - flag (T/F) to print debugging messages
#'
#'@return named vector with selectivity values at the elements of z
#'
#'@details None.
#'
#'@export
#'
plogis<-function(z,z50,sd,fsz=0){
    if (debug) message(paste('starting plogis: z50, sd =',z50,sd));
    res<-1.0/(1.0+exp(-(z-z50)/sd));
    scl <-1;
    if (fsz>0){
        scl<-(1.0+exp(-(fsz-z50)/sd));
    } else if (fsz<0){
        scl<-1/max(res);
    }
    res<-scl*res;
    names(res)<-as.character(z);
    if (debug) message(paste('finished plogis. res =',res));
    return(res);
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and slope
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and slope.
#'
#'@param z     - vector of sizes at which to compute selectivities
#'@param z50   - size at which selectivity  = 0.5 (logit-scale mean)
#'@param slope - slope at z50
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'@param debug - flag (T/F) to print debugging messages
#'
#' @return named vector with selectivity values at the elements of z
#'
#'@details None.
#'
#'@export
#'
asclogistic<-function(z,z50,slope,fsz=0,debug=FALSE){
    #cat(z,'\n')
    #cat('z50, lnD = ',z50,lnD,'\n')
    res <- 1.0/(1.0+exp(-slope*(z-z50)));
    scl <-1;
    if (fsz>0){
        scl<-(1.0+exp(-slope*(fsz-z50)));
    } else if (fsz<0){
        scl<-1.0/max(res);
    }
    res<-scl*res;
    #print(res);
    names(res)<-as.character(z);
    #print(res)
    return(res)
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and z95
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and z95.
#'
#'@param z   - vector of sizes at which to compute selectivities
#'@param z50 - size at which selectivity  = 0.5 (logit-scale mean)
#'@param z95 - size at which selectivity  = 0.95
#'@param fsz - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'@param debug - flag (T/F) to print debugging messages
#'
#' @return named vector with selectivity values at the elements of z
#'
#'@details None.
#'
#'@export
#'
asclogistic5095<-function(z,z50,z95,fsz=0){
    slope<-log(19.0)/(z95-z50);
    return(asclogistic(z,z50,slope,fsz=fsz));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and D95 (=z95-z50)
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and D95.
#'
#'@param z   - vector of sizes at which to compute selectivities
#'@param z50 - size at which selectivity  = 0.5 (logit-scale mean)
#'@param D95 - z95-z50
#'@param refZ - reference size
#'@param debug - flag (T/F) to print debugging messages
#'
#'@return named vector with selectivity values at the elements of z
#'
#'@details None.
#'
#'@export
#'
asclogistic50D95<-function(z,z50,D95,refZ=0,debug=FALSE){
    slope<-log(19.0)/params[2];
    return(asclogistic(z,c(params[1],slope),refZ,debug));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate a double logistic selectivity curve
#'@description Function to calculate a double logistic selectivity curve.
#'@param z - vector of sizes at which to compute selectivity curve
#'@param params - 4-element vector of selectivity function parameters
#'@param refZ - reference size
#'@param debug - flag (T/F) to print debugging messages
#'
#'@return named vector with selectivity values at the elements of z
#'
#'@details The parameter values are
#'\itemize{
#' \item{params[1] - ascZ50   - ascending logistic size at which selectivity  = 0.5 (logit-scale mean)}
#' \item{params[2] - ascSlope - ascending logistic slope at z50}
#' \item{params[3] - dscZ50   - descending logistic size at which selectivity  = 0.5 (logit-scale mean)}
#' \item{params[4] - dscSlope - descending logistic slope at z50}
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
#'@title Calculate a double logistic function parameterized by z50 and D95=z95-z50 for ascending/descending limbs
#'
#'@description Function to calculate a double logistic function parameterized by z50, z95 on ascending/descending limbs.
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
  sel    = ascJ*ascN+(1.0-ascJ);
  if (debug) message(paste("Finished ascnormal(...)"));
  return(sel);
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
