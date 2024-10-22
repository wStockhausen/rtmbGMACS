
#'
#' @title Generic S3 method to expand an object by another
#' @description Generic S3 method to expand an object by another.
#' @param x - object to expand
#' @param y - object to expand by
#' @return object of `class(x)` expanded by object y
#' @export
#'
expand<-function(x,y){UseMethod("expand")}
