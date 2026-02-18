#' @title Get valid "smooth" class names
#' @return character vector with names of "smooths" in .valid_smooths
#' @export
findSmTrmClasses <- function() {
    .valid_smooths;
}

#' @title Test if formula includes `mgcv`-type smooths
#' @description Function to test if formula includes `mgcv`-type smooths.
#' @param formula - formula to test
#' @return logical (T/F)
#'
#' @examplesIf FALSE
#' # example code
#' hasSMs(~1+a:b+(r|x)+s(x));
#' hasSMs(~a:b+s:x);      #--FEs only
#' hasSMs(~(r|x)+s(x));   #--REs plus smooth
#' hasSms(~s(x));         #--smooth only
#' hasSMs(~ar1(r|x));     #--"special" REs
#'
#'@importFrom reformulas anySpecial
#'
#' @export
#'
hasSMs<-function(formula){
  return(reformulas::anySpecial(formula, specials=.valid_smooths));
}

#' @title Extract smooth terms from a formula
#' @param formula - formula (or term) from which to extract smooth terms
#' @param targets - character vector of smooth function names (default=.valid_smooths)
#' @param debug - flag to print debugging info
#' @return formula with only smooth terms, or list of smooth terms, or NULL if no smooths
#' @details Valid smooths are found in rtmbGMACS:::.valid_smooths.
#'
#' @examplesIf FALSE
#' # example code
#' f = ~ s(z) + s:x;
#' sms = getTermsSMs(f,debug=TRUE);
#'
#' f = ~ s:x;            #--no smooths
#' sms = getTermsSMs(f);
#'
#' f = ~ 1 + s:x + s(x, bs = "tp") + s(y, bs = "gp") + ti(z) + (g|p)
#' sms = getTermsSMs(f);
#' sms = getTermsSMs(f,targets="ti");
#'
#' f = ~ ti(x,y, bs = "gp") + s(z) + (g|p)
#' sms = getTermsSMs(f);
#'
#' @export
#'
getTermsSMs<-function(trm,targets=.valid_smooths,spc="",debug=FALSE){
  nt = length(trm);
  if (debug) {cat(spc,"length(term) =",nt,"\n");}
  if (nt>1){
    if (debug) cat(spc,"term 1: '",deparse(trm[[1]]),"'\n",sep="");
    spc=paste0(spc," ");
    if (deparse(trm[[1]]) %in% targets) {
      if (debug) cat(spc,"!!found match in term 1: '",deparse(trm[[1]]),"'\n",sep="");
      return(trm);
    }

    lst = list();
    if (debug) cat(spc,"no match for term 1\n",sep="");
    for (i in 2:nt){
      if (debug) cat(spc,"checking term ",i,": '",deparse(trm[[i]]),"'\n",sep="");
      f2 = getTermsSMs(trm[[i]],targets=targets,spc=paste0(spc," "),debug=debug);
      lst = c(lst,f2);
    }
    if (length(lst)==0) lst=NULL;
    return(lst);
  }
  return(NULL);
}

#' @title Extract a formula for smooth terms from a formula
#' @param formula - formula (or term) from which to extract smooth terms
#' @param targets - character vector of smooth function names (default=.valid_smooths)
#' @param returnTerms - flag to return a list of individual terms rather than a formula
#' @param debug - flag to print debugging info
#' @return formula with only smooth terms, or list of smooth terms, or NULL if no smooths
#' @details Valid smooths are found in rtmbGMACS:::.valid_smooths.
#'
#' @examplesIf FALSE
#' # example code
#' f = ~ s(z) + s:x;
#' sms = getSMs(f,targets="s",debug=TRUE);
#'
#' f = ~ s:x;
#' sms = getSMs(f);
#'
#' f = ~ 1 + s:x + s(x, bs = "tp") + ti(y, bs = "gp") + s(z) + (g|p)
#' sms = getSMs(f,returnTerms=TRUE);
#'
#' f = ~ s(x,y, bs = "gp") + ti(z) + (g|p)
#' sms = getSMs(f,"s");
#'
#' @importFrom reformulas sumTerms
#' @export
#'
getSMs <- function(formula,
                    targets = .valid_smooths,
                    returnTerms=FALSE,
                    debug=FALSE) {
  #--get individual smooth terms as list
  lst = getTermsSMs(formula,targets=targets,spc="",debug=debug);

  if (!is.null(lst)) {
    if (returnTerms) return(lst);
    frmla = reformulas::sumTerms(c(0,lst));#--add in zero intercept
    return(makeRHSFormula(frmla))
  }
  return(NULL);
}

##' Parse a formula into smooth terms,
##' treating 'special' terms (of the form foo(x|g\[,m\])) appropriately
##'
##' Hacked from [reformulas::splitForm()] to consider only RE terms.
##' @title Split formula containing special random effect terms
##' @param formula a formula containing special random effect terms
##' @param debug flag (TRUE/FALSE) to print debugging info
##' @param sm_fcns - vector of names of valid "smooth" functions
##' @return a list containing elements
##' \itemize{
##'  \item{\code{smTrmFormulas} list of \code{x | g} formulas for each term}
##'   \item{\code{smTrmAddArgs} list of function+additional arguments, i.e. \code{list()} (non-special), \code{foo()} (no additional arguments), \code{foo(addArgs)} (additional arguments)}
##'   \item{\code{smTrmClasses} (vector of special functions/classes, as character)}
##' }
##' @examples
##' splitForm_SM(~x+y)                     ## fixed effects only
##' splitForm_SM(~x+y+(f|g))               ## fixed and random effects
##' splitForm_SM(~x+y+diag(f|g))           ## one special
##' splitForm_SM(~1+diag(f|g)+s(x, bs = "tp"))
##' splitForm_SM(~1+diag(f|g)+s(x, k=7, bs = "tp"))
##'
##' @author Steve Walker
##' @export
splitForm_SM <- function(formula,
                          debug=FALSE,
                          sm_fcns = findSmTrmClasses()) {

    ## split formula into separate "smooth" terms

    fbxx <- findbars_x(formula, debug, sm_fcns)
    formSplits <- expandAllGrpVar(fbxx)

    if (length(formSplits)>0) {
        formSplitID <- sapply(lapply(formSplits, "[[", 1), as.character)
                                        # warn about terms without a
                                        # setReTrm method

        parenTerm <- formSplitID == "("
                                        # capture additional arguments
        smTrmAddArgs <- lapply(formSplits, "[", -2)[!parenTerm]
                                        # remove these additional
                                        # arguments
        formSplits <- lapply(formSplits, "[", 1:2)
                                        # standard RE terms
        formSplitStan <- formSplits[parenTerm]
                                        # structured RE terms
        formSplitSpec <- formSplits[!parenTerm]

        smTrmFormulas <- c(lapply(formSplitStan, "[[", 2),
                           lapply(formSplitSpec, "[[", 2))
        smTrmFormulas <- unlist(smTrmFormulas) # Fix me:: added for rr structure when it has n = 2, gives a list of list... quick fix
        # smTrmClasses <- c(rep(defaultTerm, length(formSplitStan)),
        #                   sapply(lapply(formSplitSpec, "[[", 1), as.character))
        smTrmClasses <- c(sapply(lapply(formSplitSpec, "[[", 1), as.character))
    } else {
        smTrmFormulas <- smTrmAddArgs <- smTrmClasses <- NULL
    }

    list(smTrmFormulas = smTrmFormulas,
         smTrmAddArgs  = smTrmAddArgs,
         smTrmClasses  = smTrmClasses)
}

