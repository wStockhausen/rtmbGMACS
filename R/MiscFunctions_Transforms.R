#'
#' @title Return the input with no changes
#' @description A pass-through transformation function (doesn't change values of input).
#' @param x - vector of values to "convert"
#' @return vector of "transformed" values
#' @details This constitutes a function that simply returns the input.
#' @export
#'
tf_none<-function(x){
  return(x);
}
#'
#' @title Convert values from grams to kilograms
#' @description Function to convert values from grams to kilograms.
#' @param x - vector of values to convert
#' @return vector of transformed values
#' @export
#'
tf_g2kg<-function(x){
  return(x/1000.0);
}
