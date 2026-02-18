#' @title Test if formula includes fixed effects
#' @description Function to test if formula includes fixed effects.
#' @param formula - formula to test
#' @return logical (T/F)
#' @examplesIf FALSE
#' # example code
#' hasFEs(~1+a:b+(r|x)+s(x));
#' hasFEs(~a:b);
#' hasFEs(~(r|x)+s(x));   #--has implied intercept
#' hasFEs(~0+(r|x)+s(x)); #--explicitly has no intercept
#'
#' @export
#'
hasFEs<-function(formula){
  fixedform <- getFEs(formula);
  if (is.language(fixedform[[2]])) return(TRUE);
  return(fixedform[[2]]!=0);
}

#' @title Get terms in a fixed effects formula
#' @description Function to get terms in a fixed effects formula.
#' @param formula - formula to extract terms from
#' @param spc - spacer for debug printing
#' @param debug - flag to print debugging info
#' @return list of terms from fixed effects formula in REVERSE order
#' @examplesIf FALSE
#' # example code
#' getTermsFEs(~a*b+1);
#' getTermsFEs(~1);
#' getTermsFEs(~0+a:b+x);
#' getTermsFEs(~a:b + x + 0);
#'
#' @export
#'
getTermsFEs<-function(trm,spc="",debug=TRUE){
    nt = length(trm);
    if (debug) {cat(spc,"length(term) =",nt,"\n");}
    if (is.symbol(trm)) {
      if (debug) cat(spc,"term: '",deparse(trm),"' is a symbol\n",sep="");
      return(trm);
    }
    if (debug) cat(spc,"term 1: '",deparse(trm[[1]]),"'\n",sep="");
    if (deparse(trm[[1]]) %in% c("*",":","/")){
      if (debug) cat(spc,"term 1: '",deparse(trm[[1]]),"'\n",sep="");
      return(trm);
    } else
    if (deparse(trm[[1]]) %in% c(0,1)){
      if (debug) cat(spc,"term 1: '",deparse(trm[[1]]),"'\n",sep="");
      return(trm);
    } else
    if (deparse(trm[[1]])=="~") {
      return(getTermsFEs(trm[[2]],spc=paste0(spc," "),debug=debug));
    } else
    if (deparse(trm[[1]])=="+") {
      lst = list();
      if (debug) {
        cat(spc,"term ",2,": '",deparse(trm[[2]]),"'\n",sep="");
        cat(spc,"keeping term ",3,": '",deparse(trm[[3]]),"'\n",sep="");
      }
      lst = c(lst,trm[[3]]);
      f2 = getTermsFEs(trm[[2]],spc=paste0(spc," "));
      lst = c(lst,f2);
      return(lst);
    }
    return(NULL);
}

#' @title Get fixed effects portion of formula
#' @description Function to get fixed effects portion of formula.
#' @param formula - formula to extract fixed effects from
#' @return fixed effects formula
#' @examplesIf FALSE
#' # example code
#' getFEs(~1+a:b+(r|x)+s(x));
#' getFEs(~a:b);
#' getFEs(~(r|x)+s(x));   #--has implied intercept
#' getFEs(~0+(r|x)+s(x)); #--explicitly has no intercept
#' getFEs(~1+a:b+(r|x)+s(x),returnTerms=TRUE);
#'
#' @importFrom reformulas nobars noSpecials RHSForm
#'
#' @export
#'
getFEs<-function(formula,returnTerms=FALSE,debug=FALSE){
  fixedform <- formula;
  ###--*remove* RE terms
  reformulas::RHSForm(fixedform) <- reformulas::nobars(reformulas::RHSForm(fixedform));
  ###--*remove* s() terms from fixed formula-
  fixedform <- reformulas::noSpecials(fixedform, specials = .valid_smooths);
  if (returnTerms && !is.null(fixedform)) return(rev(getTermsFEs(fixedform,debug=debug)));
  return(fixedform);
}


#' @title Get parameter information for fixed effects
#' @description Function to get parameter information for fixed effects.
#' @param formula - a formula from which to extract parameter info for fixed effects
#' @param model_frame - the model frame
#' @param contrasts - a list of contrasts (see ?glmmTMB)
#' @param offset - the name of a column in `fr` to be used as a fixed offset
#' @param sparse - (logical) return sparse model matrix?
#' @return a list composed of
#' \itemize{
#'  \item{p - a suitable parameter vector for the fixed effects}
#'  \item{X}{design matrix for fixed effects}
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
createParamsInfo_FEs <- function(formula, model_frame, contrasts = NULL, offset = NULL, sparse=TRUE) {

  fixedform <- getFEs(formula);
  nobs <- nrow(model_frame)
  if (is.null(fixedform)){
    #--no fixed effects, return NULL (TODO: better idea??)
    return(NULL);
  }

  ##--determine FE model matrix----
  if (!sparse) {
      X <- model.matrix(reformulas::noSpecials(fixedform, specials = "offset"),
                                               model_frame, contrasts);
  } else {
      X <- Matrix::sparse.model.matrix(reformulas::noSpecials(fixedform, specials = "offset"),
                                                              model_frame, contrasts);
  }

  ## will be 0-column matrix if fixed formula is empty
  offset <- rep(0,nobs);
  if (reformulas::inForm(fixedform,quote(offset))) {
      ## hate to match offset terms with model frame names
      ##  via deparse, but since that was what was presumably done
      ##  internally to get the model frame names in the first place ...
      for (o in reformulas::extractForm(fixedform,quote(offset))) {
          offset_nm <- deparse1(o);
          ## don't think this will happen, but ...
          if (length(offset_nm)>1) {
              stop("trouble reconstructing offset name")
          }
          offset <- offset + model_frame[[offset_nm]];
      }
  }

  p           = RTMB::AD(numeric(ncol(X))); #--placeholder for FE parameters
  fe_form     = fixedform;
  model_matrix = dplyr::bind_cols(model_frame,as.matrix(X));
  return(namedList(p,
                   X,
                   fe_form,
                   contrasts,
                   offset,
                   formula,
                   model_frame,
                   model_matrix));
} #--CreateParamsInfo_FEs

