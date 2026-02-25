#--transforms and inverse transforms

#'@title list of valid inverse transforms
#'@details Given the function name `fcn` as a string, its inverse transform function is obtained by
#' `.valid_invtforms[fcn]`.
.valid_invtforms = list(
  exp=log,
  mfexp=log,
  log=exp,
  g2kg=kg2g,
  kg2g=g2kg,
  kg2kt=kt2kg,
  kt2kg=kg2kt,
  ident=ident,
  identity=identity,
  none=none,
  logistic=logit,
  logit=logistic,
  symlogistic=symlogit,
  symlogit=symlogistic
)

#' @title Compute transforms on a vector of values
#' @description Function to compute transforms on a vector of values.
#' @param fcn - vector of function names to apply
#' @param val - vector of values to apply transform(s) to
#' @return numeric (or adsparse) vector of transformed values
#' @details If the length of `fcn` is not equal to the length of `val`, it
#' is repeated until this condition is reached.
#' @examplesIf FALSE
#' # example code
#' tfcns = c('symlogistic','symlogit',   'ident','log');
#' ifcns = c('symlogit',   'symlogistic','ident','exp');
#' vals = c(0.5,0.5,10,10);
#' rest = tform(tfcns,vals);
#' tform(ifcns,rest);
#' tform_inv(fcns,rest);
#' @importFrom rlang env_get
#' @importFrom RTMB AD
#' @export
#'
tform<-function(fcn,val,force=FALSE){
  nv = length(val);
  res = RTMB::AD(numeric(nv),force=force);
  if (length(fcn)<nv) rep_len(fcn,nv);
  for (i in 1:nv){
    res[i] = rlang::env_get(nm=fcn[i],inherit=TRUE)(val[i]);
  }
  return(res);
}


#' @title Compute inverse transforms on a vector of values
#' @description Function to compute inverse transforms on a vector of values.
#' @param fcn - vector of function names used to apply "forward" transform
#' @param val - vector of values to apply inverse transform(s) to
#' @return numeric (or adsparse) vector of transformed values
#' @details If the length of `fcn` is not equal to the length of `val`, it
#' is repeated until this condition is reached.
#' @examplesIf FALSE
#' # example code
#' tfcns = c('symlogistic','symlogit',   'ident','log');
#' ifcns = c('symlogit',   'symlogistic','ident','exp');
#' vals = c(0.5,0.5,10,10);
#' rest = tform(tfcns,vals);
#' tform(ifcns,rest);
#' tform_inv(fcns,rest);
#' rest = tform(tfcns,vals,force=TRUE);
#' tform(ifcns,rest,force=TRUE);
#' tform_inv(fcns,rest,force=TRUE);
#' @importFrom RTMB AD
#' @export
#'
tform_inv<-function(fcn,val,force=FALSE){
  nv = length(val);
  res = RTMB::AD(numeric(nv),force=force);
  if (length(fcn)<nv) fcn = rep_len(fcn,nv);
  for (i in 1:nv){
    res[i] = .valid_invtforms[[fcn[i]]](val[i]);
  }
  return(res);
}
