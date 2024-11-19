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
