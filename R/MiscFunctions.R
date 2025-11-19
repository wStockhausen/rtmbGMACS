#--miscellaneuos functions, including ADMB replacements

#'
#' @title Return the input with no changes
#' @description A pass-through transformation function (doesn't change values of input).
#' @param x - vector of values to "convert"
#' @return vector of "transformed" values
#' @details This constitutes a function that simply returns the input.
#' @export
#'
getLinkFcn<-function(txt){
  if (tolower(txt) %in% c("ident","identity","add","none"))
    return(list(link=identity,link_inv=identity));
  if (tolower(txt) %in% c("log"))
    return(list(link=log,link_inv=exp));
  if (tolower(txt) %in% c("logit"))
    return(list(link=logit,link_inv=logistic));
}

#' @title Zero function
#' @description Function that simply returns the input object multiplied by zero.
#' @param x - a numeric object to "convert"
#' @return a numeric object of the same class and shape as the input filled with zeros
#' @details This function simply returns a numeric object of the same class and shape as the input filled with zeros.
#' @export
#'
zero<-function(x){return(0.0*x)}

#' @title Identity function
#' @description Function that simply returns the input object
#' @param x - object to "convert"
#' @return the input object
#' @details This function that simply returns the input.
#' @export
#'
identity<-function(x){return(x);}

#'
#' @title Convert values to the logit scale
#' @description Function to convert values to the logit scale.
#' @param x - numeric object
#' @return log(x/(1-x)), with same names and class as `x`
#' @details inverse of logistic function. Values of `x` outside the
#' (0,1) range generate NaNs.
#' @export
#'
logit<-function(x){
  y = log(x/(1-x));
  names(y) = names(x);
  return(y);
}

#'
#' @title Square a variable
#' @description ADMB-equivalent function to square a variable.
#' @param x - variable to square
#' @return vector x*x with same names as x
#' @details Substitution for ADMB `square` function.
#' @export
#'
square<-function(x){
  y = x*x;
  names(y) = names(x);
  return(y);
}

#'
#' @title Exponentiate a variable
#' @description ADMB-equivalent  to exponentiate a variable.
#' @param x - variable to exponentiate
#' @return vector x*x with same names as x
#' @details Substitution for ADMB `mfexp` function.
#' @export
#'
mfexp<-function(x){
  y = exp(x);
  names(y) = names(x);
  return(y);
}

#'
#' @title Calculate the element-wise product of two variables
#' @description ADMB-equivalent to calculate the element-wise product of two variables.
#' @param x - variable (vector, matrix, etc.)
#' @param y - variable conformable to x
#' @return (unnamed) object x*y
#' @details Substitution for ADMB `elem_prod` function.
#' @export
#'
elem_prod<-function(x,y){
  #--TBD: check dimensions are conformable?
  z = x * y;
  return(z);
}

#'
#' @title Calculate the element-wise division of two variables
#' @description ADMB-equivalent to calculate the element-wise division of two variables.
#' @param x - numerator variable (vector, matrix, etc.)
#' @param y - divisor variable conformable to x
#' @return (unnamed) object x/y
#' @details Substitution for ADMB `elem_div` function.
#' @export
#'
elem_div<-function(x,y){
  #--TBD: check dimensions are conformable?
  z = x / y;
  return(z);
}

#'
#' @title Create a "square wave" opened to the right
#' @description Function to create a "square wave" opened to the right.
#' @param z0 - smallest size at which function is non-zero
#' @param - zs - object with sizes at which to calculate the function
#' @param dz - width to adjust zs by (typically, binwidth; default=1)
#' @param shftfac - factor to scale dz by when adjusting location by so w(z==z0) = 1 (default=20)
#' @param sclfac - factor to scale dz by to obtain `plogis` scale
#' @return object same size as zs
#' @details Provides a differentiable function such that `w(z>=z0) = 1` and
#' `w(z<z0) = 0`. Uses [plogis()] to calculate the logistic curve, with
#' `location = z0` and `scale = dz/sclfac`. The `zs` are shifted to the left
#' by `dz/shftfac` so that w(z==z0) = 1 and w(z>z0) = 1.
#' @export
#'
squarewave_right<-function(z0,zs,dz=1,shftfac=2,sclfac=20){
  # cat("in square_wave_right\n");
  # cat("z0",is.numeric(z0),"\n");
  # cat("zs",is.numeric(zs),"\n");
  # cat("dz",is.numeric(dz),"\n");
  # cat("shftfac",is.numeric(shftfac),"\n");
  # cat("sclfac",is.numeric(sclfac),"\n");
  q = (zs+dz/shftfac - z0)/(dz/sclfac);#--rescale for RTMB plogis
  # cat("q",is.numeric(q),"\n");
  # cat("zs",zs,"\n");
  # cat("q:",q,"\n");
  w = plogis(q);
  # cat("sqw_R:",w,"\n");
  return(w);
}

#'
#' @title Create a "square wave" opened to the left
#' @description Function to create a "square wave" opened to the left.
#' @param z0 - largest size at which function is non-zero
#' @param - zs - object with sizes at which to calculate the function
#' @param dz - width to adjust zs by (typically, binwidth; default=1)
#' @param shftfac - factor to scale dz by when adjusting location by so w(z==z0) = 1 (default=20)
#' @param sclfac - factor to scale dz by to obtain `plogis` scale
#' @return object same size as zs
#' @details Provides a differentiable function such that `w(z<=z0) = 1` and
#' `w(z>z0) = 0`. Uses [plogis()] to calculate the logistic curve, with
#' `location = z0` and `scale = dz/sclfac`. The `zs` are shifted to the right
#' by `dz/shftfac` so that w(z==z0) = 1 and w(z<z0) = 1.
#' @export
#'
squarewave_left<-function(z0,zs,dz=1,shftfac=2,sclfac=20){
  # cat("in square_wave_left\n");
  # cat("z0",is.numeric(z0),"\n");
  # cat("zs",is.numeric(zs),"\n");
  # cat("dz",is.numeric(dz),"\n");
  # cat("shftfac",is.numeric(shftfac),"\n");
  # cat("sclfac",is.numeric(sclfac),"\n");
  q = (zs-dz/shftfac - z0)/(dz/sclfac);#--rescale for RTMB plogis
  # cat("zs",zs,"\n");
  # cat("q:",q,"\n");
  w = 1-plogis(q);
  # cat("sqw_L:",w,"\n")
  return(w);
}

#'
#' @title Create a "square wave"
#' @description Function to create a "square wave".
#' @param zL - smallest size at which function is non-zero
#' @param zR - largest size at which function is non-zero
#' @param -zs - object with sizes at which to calculate the function
#' @param dz - width to adjust zs by (typically, binwidth; default=1)
#' @param shftfac - factor to scale dz by when adjusting location by so w(z==z0) = 1 (default=20)
#' @param sclfac - factor to scale dz by to obtain `plogis` scale
#' @return object same size as zs
#' @details Provides a differentiable function such that `w(zL<=z<=zR) = 1` and
#' 0 otherwise.
#' @export
#'
squarewave<-function(zL,zR,zs,dz=1,shftfac=2,sclfac=20){
  w = squarewave_right(zL,zs,dz=dz,shftfac=shftfac,sclfac=sclfac)*
        squarewave_left(zR,zs,dz=dz,shftfac=shftfac,sclfac=sclfac);
  return(w);
}

addList<-function(lst1,lst2){
  if (is.null(lst2)) return(lst1); #--do nothing to lst1
  if (is.null(lst1)) {
    if (class(lst2)!="list") stop(paste0("lst2 is class '",class(lst2),"'. Class must be 'list'."));
    return(lst2); #--return lst2
  } else {
    if (class(lst1)!="list") stop(paste0("lst1 is class '",class(lst1),"'. Class must be 'list'."));
    if (class(lst2)!="list") stop(paste0("lst2 is class '",class(lst2),"'. Class must be 'list'."));
    for (nm in names(lst2)){
      lst1[[nm]] = lst2[[nm]];
    }
  }
  return(lst1);
}
