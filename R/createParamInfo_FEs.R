#' @title Get parameter information for fixed effects
#' @description Function to get parameter information for fixed effects.
#' @param formula - a formula from which to extract parameter info for fixed effects
#' @param model_frame - the model frame
#' @param contrasts - a list of contrasts
#' @param drop.unused.levels - flag (T/F) to drop unused levels
#' @param offset - the name of a column in `fr` to be used as a fixed offset
#' @param sparse - flag (T/F) to return sparse model matrix?
#' @return a list composed of
#' \itemize{
#'  \item{p - a suitable parameter vector for the fixed effects}
#'  \item{X}{design matrix for fixed effects}
#'  \item{fe_pvs_lbls}{labels for fixed effects parameter values}
#'  \item{fe_form}{the fixed effects formula}
#'  \item{contrasts}{the list of specified contrasts (if any)}
#'  \item{offset}{offset vector, or vector of zeros if offset not specified}
#'  \item{formula}{input formula}
#'  \item{model_frame - the model frame}
#'  \item{model_matrix - combined model frame and X}
#' }
#'
#' @importFrom stats model.matrix contrasts
#' @importFrom reformulas extractForm inForm noSpecials
#' @importFrom RTMB AD
#'
#' @export
#'
createParamInfo_FEs<-function(formula,
                                model_frame,
                                contrasts = NULL,
                                drop.unused.levels=FALSE,
                                offset = NULL,
                                sparse=TRUE) {

  fe_form <- getFEs(formula);
  nobs <- nrow(model_frame)
  if (is.null(fe_form)){
    #--no fixed effects, return NULL (TODO: better idea??)
    return(NULL);
  }

  ##--determine FE model matrix----
  if (!sparse) {
      X <- model.matrix(reformulas::noSpecials(fe_form, specials = "offset"),
                                               model_frame, contrasts);
  } else {
      X <- Matrix::sparse.model.matrix(reformulas::noSpecials(fe_form, specials = "offset"),
                                                              model_frame, contrasts);
  }
  if (drop.unused.levels){
    idx = which(abs(colSums(as.matrix(X)))>0);
    X = X[,idx];
  }

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
          offset <- offset + model_frame[[offset_nm]];
      }
  }

  p           = RTMB::AD(numeric(ncol(X))); #--placeholder for FE parameters
  model_matrix = dplyr::bind_cols(model_frame,as.matrix(X));
  return(namedList(p,
                   X,
                   fe_pvs_lbls, #--labels for parameter values
                   fe_form,
                   contrasts,
                   offset,
                   formula,
                   model_frame,
                   model_matrix));
} #--CreateParamInfo_FEs

