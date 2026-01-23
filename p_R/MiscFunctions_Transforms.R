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
#'
#' @title Transform values to log-space
#' @description Function to convert values to logspace.
#' @param x - vector of values to convert
#' @return vector of transformed values
#' @export
#'
tf_log<-function(x){return(log(x));}


#'
#' @title Convert values using a transform function
#' @description Function to convert values sing a transform function.
#' @param tform - name of transform function
#' @param x - vector of values to convert
#' @return vector of transformed values
#' @export
#'
tf_apply<-function(tform,x){
  tf = eval(parse(text=paste0(ifelse(stringr::str_starts(tform,"tf_"),tform,paste0("tf_",tform)),
                              "(",x,")")));
  return(tf);
}
