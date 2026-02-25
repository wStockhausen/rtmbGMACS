#' @title Create a function info object
#' @description Function to create a function info object.
#' @param fcn_id_ - label used to identify particular instance of the function
#' @param fcn_ - the name of a function
#' @param model_frame -  the model frame which the function applies to
#' @param vars - function variable (as character string)
#' @param par_ids - character vector with par_ids
#' @return a dataframe (at present--this will probably change)
#' @details The name of the function and the function id are appended to the model-frame. The
#' same function name, `fcn`, can be associated with multiple `fcn_id`'s, hence the necessity
#' for supplying a unique function "label".
#' @importFrom dplyr cross_join mutate
#' @importFrom tibble tibble
#' @export
createFunctionInfo<-function(fcn_id_,fcn_,model_frame,vars_,par_ids_){
  #--TODO: lots to add in here: function-level REs, priors, etc.
  mf = model_frame |> dplyr::cross_join(tibble::tibble(fcn_id=fcn_id_,fcn=fcn_,vars=vars_,par_ids=list(par_ids_)));
  return(mf);
}
