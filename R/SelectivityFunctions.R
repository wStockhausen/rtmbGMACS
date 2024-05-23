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
#'@title Calculate the logistic function
#'
#'@description Function to calculate the logistic function.
#'
#'@param z   - vector of size class midpoints at which to compute selectivity
#'@param z50 - size at which selectivity = 0.5 (logit-scale mean)
#'@param sd  - standard deviation in selectivity (logit-scale standard deviation)
#'@param fsz - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
#'@details None.
#'
#'@export
#'
plogis<-function(z,z50,sd,fsz=0){
    #cat(z,'\n')
    #cat('z50, sd = ',z50,sd,'\n')
    res<-1.0/(1.0+exp(-(z-z50)/sd));
    scl <-1;
    if (fsz>0){
        scl<-(1.0+exp(-(fsz-z50)/sd));
    } else if (fsz<0){
        scl<-1/max(res);
    }
    res<-scl*res;
    #print(res);
    names(res)<-as.character(z);
    #print(res)
    return(res)
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
#'
#'@return vector with selectivity values at the elements of z
#'
#'@details None.
#'
#'@export
#'
asclogistic<-function(z,z50,slope,fsz=0){
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
#'
#'@return vector with selectivity values at the elements of z
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
#'@param fsz - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
#'@details None.
#'
#'@export
#'
asclogistic50D95<-function(z,z50,D95,fsz=0){
    slope<-log(19.0)/D95;
    return(asclogistic(z,z50,slope,fsz=fsz));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate a double logistic function parameterized by z50 and slope for ascending/descending limbs
#'
#'@description Function to calculate a double logistic function parameterized by z50 and slope for ascending/descending limbs.
#'
#'@param z     - vector of sizes at which to compute selectivities
#'@param ascZ50   - ascending logistic size at which selectivity  = 0.5 (logit-scale mean)
#'@param ascSlope - ascending logistic slope at z50
#'@param dscZ50   - descending logistic size at which selectivity  = 0.5 (logit-scale mean)
#'@param dscSlope - descending logistic slope at z50
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
#'@details None.
#'
#'@export
#'
dbllogistic<-function(z,ascZ50,ascSlope,dscZ50,dscSlope,fsz=0){
    #cat(z,'\n')
    #cat('z50, lnD = ',z50,lnD,'\n')
    res <- 1.0/(1.0+exp(-ascSlope*(z-ascZ50)))*1.0/(1.0+exp(dscSlope*(z-dscZ50)));
    scl <-1;
    if (fsz>0){
        scl<-(1.0+exp(-ascSlope*(fsz-ascZ50)))*(1.0+exp(dscSlope*(fsz-dscZ50)));
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
#'@title Calculate a double logistic function parameterized by z50 and D95=z95-z50 for ascending/descending limbs
#'
#'@description Function to calculate a double logistic function parameterized by z50, z95 on ascending/descending limbs.
#'
#'@param z      - vector of sizes at which to compute selectivities
#'@param ascZ50 - ascending z50
#'@param ascZ95 - z95 on ascending limb
#'@param dscZ95 - z95 on descending limb
#'@param dscZ50 - descending z50
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
#'@details None.
#'
#'@export
#'
dbllogistic5095<-function(z,ascZ50,ascZ95,dscZ95,dscZ50,fsz=0){
    ascSlope<-log(19.0)/(ascZ95-ascZ50);
    dscSlope<-log(19.0)/(dscZ50-dscZ95);
    return(dbllogistic(z,ascZ50,ascSlope,dscZ50,dscSlope,fsz));
}

#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending normal selectivity function parameterized by size-at-1 and width
#'
#'@description Function to calculate an ascending normal selectivity function parameterized by size-at-1 and width.
#'
#'@param z      - vector of sizes at which to compute selectivities
#'@param params - selectivity parameters
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
#'@details The parameters are
#' \itemize{
#' \item{p[1] - size at which ascending limb hits 1 (mean of normal distribution)}
#' \item{p[2] - width of ascending limb (standard deviation ofnormal distribution)}
#' }
#'
#'@export
#'
ascnormal<-function(z,params,fsz=0){
  slp = 5.0;
  ascMnZ = params[1];#--size at which ascending limb hits 1
  ascWdZ = params[2];#--width of ascending limb
  ascN   = exp(-0.5*square((z-ascMnZ)/ascWdZ));
  ascJ   = 1.0/(1.0+mfexp(slp*(z-(ascMnZ))));
  sel    = ascJ*ascN+(1.0-ascJ);
  return(sel);
}
