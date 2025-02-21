#--miscellaneuos functions, including ADMB replacements

#'
#' @title Calculate variable on logit scale
#' @description Function to calculate variable on logit scale.
#' @param x - variable
#' @return vector log(x/(1-x)) with same names as x
#' @details inverse of logisitc function
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
#' @param dz - width to use for rise of function from 0 to 1
#' @return object same size as zs
#' @details Provides a differentiable function such that w(z<z0) = 0 and
#' w(z>z0) = 1.
#' @export
#'
squarewave_right<-function(z0,zs,dz=0.5){
  w = exp((zs-z0)/dz)/(1+exp((zs-z0)/dz));
  return(w);
}
#'
#' @title Create a "square wave" opened to the left
#' @description Function to create a "square wave" opened to the left.
#' @param z0 - largestest size at which function is non-zero
#' @param - zs - object with sizes at which to calculate the function
#' @param dz - width to use for the decline of the function from 1 to 0
#' @return object same size as zs
#' @details Provides a differentiable function such that w(z<z0) = 1 and
#' w(z>z0) = 0.
#' @export
#'
squarewave_left<-function(z0,zs,dz=0.5){
  w = exp(-(zs-z0)/dz)/(1+exp(-(zs-z0)/dz));
  return(w);
}

#'
#' @title Create a "square wave"
#' @description Function to create a "square wave".
#' @param zL - smallest size at which function is non-zero
#' @param zR - largest size at which function is non-zero
#' @param - zs - object with sizes at which to calculate the function
#' @param dz - width to use for rise of function from 0 to 1
#' @return object same size as zs
#' @details Provides a differentiable function such that w(zL<z<zR) = 1 and
#' 0 otherwise.
#' @export
#'
squarewave<-function(zL,zR,zs,dz=0.05){
  w = squarewave_right(zL,zs,dz)-squarewave_right(zR,zs,dz);
  return(w);
}
