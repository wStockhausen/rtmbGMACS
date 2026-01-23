#'
#' @title Calculate midpoints of a vector
#' @description Function to calculate the midpoints of a vector.
#' @param v - a numeric vector
#' @return a numeric vector with 1 less element than v with corresponding midpoints of intervals
#'
#' @export
#'
calcMidPoints<-function(v){
  n = length(v);
  m = 0.5*(v[2:n]+v[1:(n-1)]);
  return(m);
}
