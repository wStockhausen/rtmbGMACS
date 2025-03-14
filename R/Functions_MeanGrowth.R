#--mean growth functions
#'
#'@title Power law model 1 for mean growth
#'@description Function to calculate power law model 1 for mean growth
#'@param z - vector of pre-molt sizes (bin mid-points)
#'@param params - vector of parameters (potentially estimated)
#'@param consts - vector of constants (not estimated)
#'@param debug - debugging flag (T/F)
#'@return named vector of mean molt increments or post-molt size
#'@details Mean growth (molt increment *or* post-molt size) `g` is described as
#'a power law function of pre-molt size `z` by
#' g(z|p, c) = exp(p[1]) * z^p[2]
#' where `p` represents the parameter vector.
#' exp(`p[1]`) represents the mean growth at pre-molt size z=1.
#'
#'The interpretation of `g` as representing a molt increment or post-molt size
#'is up to the user.
#'
#'@export
#'
calcMnGrowthPwrLaw1<-function(z,params,consts=NULL,debug=FALSE){
  grA = params[1];
  grB = params[2];
  mnZs = mfexp(grA+grB*log(z));
  names(mnZs) = as.character(z);
  return(mnZs);
}

#'
#'@title Power law model 2 for mean growth
#'@description Function to calculate power law model 2 for mean growth
#'@param z - vector of pre-molt sizes (bin mid-points)
#'@param params - vector of parameters (potentially estimated)
#'@param consts - vector of constants (not estimated)
#'@param debug - debugging flag (T/F)
#'@return named vector of mean molt increments or post-molt size
#'@details Mean growth (molt increment *or* post-molt size) `g` is described as
#'a power law function of pre-molt size `z` by
#' g(z|p, c) = p[1] * (z/c[1])^(log(p[2]/p[1])/log(c[2]/c[1]))
#' where `p` represents the parameter vector and `c` represents the consts vector.
#' `p[1]` is the mean growth at pre-molt size `c[1]`.
#' `p[2]` is the mean growth at pre-molt size `c[2]`.
#'
#'The interpretation of `g` as representing a molt increment or post-molt size
#'is up to the user.
#'
#'@export
#'
calcMnGrowthPwrLaw3<-function(z,params,consts,debug=FALSE){
  grA = params[1];
  grB = params[2];
  zGrA = consts[1]; #--pre-molt size corresponding to grA as mean post-molt size
  zGrB = consts[2]; #--pre-molt size corresponding to grB as mean post-molt size
  mnZs = grA*mfexp(log(grB/grA)/log(zGrB/zGrA)*log(z/zGrA));
  names(mnZs) = as.character(z);
  return(mnZs);
}

#'
#'@title Power law model 3 for mean growth
#'@description Function to calculate power law model 3 for mean growth
#'@param z - vector of pre-molt sizes (bin mid-points)
#'@param params - vector of parameters (potentially estimated)
#'@param consts - vector of constants (not estimated)
#'@param debug - debugging flag (T/F)
#'@return named vector of mean molt increments or post-molt size
#'@details Mean growth (molt increment *or* post-molt size) `g` is described as
#'a power law function of pre-molt size `z` by
#' g(z|p, c) = p[1] * (z/c[1])^p[2]
#' where `p` represents the parameter vector and `c` represents the consts vector.
#' `p[1]` is the mean growth at pre-molt size `c[1]`.
#'
#'The interpretation of `g` as representing a molt increment or post-molt size
#'is up to the user.
#'
#'@export
#'
calcMnGrowthPwrLaw3<-function(z,params,consts,debug=FALSE){
  grA = params[1];
  grB = params[2];
  zGrA = consts[1]; #--pre-molt size corresponding to grA as mean post-molt size
  mnZs = grA*mfexp(grB*log(z/zGrA));
  names(mnZs) = as.character(z);
  return(mnZs);
}
