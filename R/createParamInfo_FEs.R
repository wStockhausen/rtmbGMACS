#' @title Get parameter information for fixed effects
#' @description Function to get parameter information for fixed effects.
#' @param formula - a formula from which to extract parameter info for fixed effects
#' @param model_frame - the model frame
#' @param contrasts - a list of contrasts
#' @param drop.unused.levels - flag (T/F) to drop unused levels
#' @param offset - the name of a column in `fr` to be used as a fixed offset
#' @param sparse - flag (T/F) to return sparse model matrix?
#' @param debug - flag (T/F) to print debugging info
#' @return a list composed of
#' \itemize{
#'  \item{formula}{input formula}
#'  \item{model_frame - the model frame}
#'  \item{fe_form}{the fixed effects formula}
#'  \item{contrasts}{the list of specified contrasts (if any)}
#'  \item{offset}{offset vector, or vector of zeros if offset not specified}
#'  \item{vars}{fixed effects variables names}
#'  \item{fe_pvs_lbls}{labels for fixed effects parameter values}
#'  \item{X}{design matrix for fixed effects}
#'  \item{model_matrix}{model matrix associated with `X`}
#' }
#'
#' @importFrom dplyr select
#' @importFrom reformulas extractForm inForm noSpecials
#' @importFrom RTMB AD
#' @importFrom stats model.matrix contrasts
#' @importFrom tidyselect all_of
#'
#' @export
#'
createParamInfo_FEs<-function(formula,
                                model_frame,
                                contrasts = NULL,
                                drop.unused.levels=FALSE,
                                offset = NULL,
                                sparse=TRUE,
                                debug=FALSE) {
  #debug=TRUE;
  #--get number of "observation" (rows) in model frame
  nobs = nrow(model_frame);
  #--extract only FE parts of a formula
  fe_form <- getFEs(formula);
  if (is.null(fe_form)){
    #--no fixed effects, return NULL (TODO: better idea??)
    return(NULL);
  }

  #--keep only columns in the model_frame corresponding to parameters
  ##--in the fixed effects formula
  vars = identifyVars(fe_form);#--model_frame columns included in FE formula
  if (debug) cat("FE formula:",deparse1(fe_form),"FE vars:",vars,"\n");
  if (length(vars)==0){
    model_matrix = model_frame; #--hack: problems joining model matrices later (better solution?)
  } else {
    model_matrix = model_frame |> dplyr::select(tidyselect::all_of(vars));
  }
  model_matrix = model_matrix |> dplyr::distinct();

  ##--determine FE design matrix----
  if (!sparse) {
      X <- model.matrix(reformulas::noSpecials(fe_form, specials = "offset"),
                                               model_matrix, contrasts);
  } else {
      X <- Matrix::sparse.model.matrix(reformulas::noSpecials(fe_form, specials = "offset"),
                                                              model_matrix, contrasts);
  }
  if (drop.unused.levels){
    idx = which(abs(colSums(as.matrix(X)))>0);
    X = X[,idx,drop=FALSE];
  }
  if (debug) str(X);
  if (debug) cat("colnames(X) = ",paste(colnames(X),collapse=" "),"\n");

  fe_pvs_lbls = createParamValLabels(fe_form,colnames(X),debug=FALSE)$revd; #--labels for parameter values
  colnames(X) = fe_pvs_lbls;

  ## will be 0-column matrix if fixed formula is empty
  offset <- rep(0,nobs);
  if (reformulas::inForm(fe_form,quote(offset))) {
      ## hate to match offset terms with model frame names
      ##  via deparse, but since that was what was presumably done
      ##  internally to get the model frame names in the first place ...
      for (o in reformulas::extractForm(fe_form,quote(offset))) {
          offset_nm <- deparse1(o);
          ## don't think this will happen, but ...
          if (length(offset_nm)>1) {
              stop("trouble reconstructing offset name")
          }
          offset <- offset + model_frame[[offset_nm]]; #--offset column might not be in model_matrix
      }
  }

  return(namedList(formula,
                   model_frame,
                   contrasts,
                   offset,
                   fe_form,
                   vars=vars,
                   fe_pvs_lbls, #--labels for parameter values
                   X,
                   model_matrix));
} #--CreateParamInfo_FEs

