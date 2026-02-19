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


