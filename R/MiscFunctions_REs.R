#' @title Convert a vector to a factor
#' @description Function to convert a vector to a factor
#' @param v - the vector to convert
#' @param char.only - flag to convert `v` only if it is a character vector
#' @details Copied from [reformulas:::makeFac()].
#' @export
#'
makeFac <- function(v,char.only=FALSE) {
    if (!is.factor(v) && (!char.only || is.character(v))) factor(v) else v
}

#' @title infix interaction operator to combine levels from two factors
#' @description Operator to combine levels from two factors.
#' @param f1 - factor
#' @param f2 - factor
#' @param order - flag to re-order levels
#' @return factor with merged levels from `f1` and `f2`
#' @details Copied from [reformulas:::`%i%`].
#' @export
`%i%` <- function(f1, f2, fix.order = TRUE) {
    if (!is.factor(f1) || !is.factor(f2)) stop("both inputs must be factors")
    f12 <- paste(f1, f2, sep = ":")
    ## explicitly specifying levels is faster in any case ...
    u <- which(!duplicated(f12))
    if (!fix.order) return(factor(f12, levels = f12[u]))
    ## deal with order of factor levels
    levs_rank <- length(levels(f2))*as.numeric(f1[u])+as.numeric(f2[u])
    return(factor(f12, levels = (f12[u])[order(levs_rank)]))
}

#' @title Get valid RE covariance structure names
#' @return character vector with names of RE covariance structures in .valid_covstruct
#' @export
findValidCovStructs <- function() {
    c(names(.valid_covstruct)); #--originally included "s"
}

#' @title Test if formula includes random effects
#' @description Function to test if formula includes random effects.
#' @param formula - formula to test
#' @return logical (T/F)
#'
#' @examplesIf FALSE
#' # example code
#' hasREs(~1+a:b+(r|x)+s(x));
#' hasREs(~a:b);          #--FEs only
#' hasREs(~(r|x)+s(x));   #--REs plus smooth
#' hasREs(~s(x));         #--smooth only
#'
#'@importFrom reformulas findbars_x
#'
#' @export
#'
hasREs<-function(formula){
  return(!is.null(reformulas::findbars_x(formula)));
}

#' @title Extract RE terms from a formula
#' @description Function to extract RE terms from a formula
#' @param term - formula or term to extract from
#' @param pass - counter for debugging
#' @param termno - counter for debugging
#' @param debug - flag to print debugginh info
#' @return a formula, list, or NULL
#'
#' @examplesIf FALSE
#' # example code
#' res = getTermsREs(~1+a:b+(r|x)+s(x));#--FEs + "ordinary" REs + a smooth
#' res = getTermsREs(~a:b);             #--FEs only
#' res = getTermsREs(~ar1(r|x)+s(x));   #--"special" REs + a smooth
#' res = getTermsREs(~s(x));            #--smooth only
#'
#' @importFrom reformulas RHSForm sumTerms
#'
#' @export
#'
getTermsREs <- function(term,
                        specials=.valid_covstruct,
                        pass=0,termno=0,debug=FALSE) {
  bar = '|';
  dbl = '||';
  if (debug) {
    cat("pass:",pass," termno:",termno,"\n");
    cat(deparse(term),"\n");
  }
    if (is.name(term) || !is.language(term)) {
      if (debug) cat("first check: returning NULL\n\n");
      return(NULL);
    }
    if (list(term[[1]]) %in% lapply(specials,as.name)) {
        if (debug) cat("special: ",deparse(term),"\n");
        return(term);
    }
    if (term[[1]] == as.name(bar)) {  ## found x | g
        if (debug) cat(sprintf("%s term: %s\n", "bar", deparse(term)));
        return(term);
    }
    if (term[[1]] == as.name(dbl)) {
        if (debug) cat(sprintf("%s term: %s\n", "double-bar", deparse(term)));
        return(term);
    }
    if (term[[1]] == as.name("(")) {  ## found (...)
        if (debug) cat("paren term:",deparse(term),"\n");
        #return(getTermsREs(term[[2]],specials=specials,pass=pass+1,termno=2,debug=debug));  #--recursion here
        return(term);
    }
    stopifnot(is.call(term))
    if (length(term) == 2) {
        ## unary operator, decompose argument
        if (debug) cat("unary operator:",deparse(term[[1]]),",",deparse(term[[2]]),"\n")
        return(getTermsREs(term[[2]],specials=specials,pass=pass+1,termno=2,debug=debug));  #--recursion here
    }
    ## binary operator, decompose both arguments
    f2 <- getTermsREs(term[[2]],specials=specials,pass=pass+1,termno=2,debug=debug); #--recursion
    f3 <- getTermsREs(term[[3]],specials=specials,pass=pass+1,termno=3,debug=debug); #--recursion

    if (debug) { cat("binary operator:",
                     deparse(term[[1]]),",",
                     deparse(term[[2]]),",",
                     deparse(term[[3]]),"\n")
                     cat("term 2: ", deparse(f2), "\n")
                     cat("term 3: ", deparse(f3), "\n")
                     cat("\n");
    }
    c(f2, f3);
}#--getTermsREs

#' @title Extract RE terms from a formula
#' @description Function to extract RE terms from a formula
#' @param term - formula or term to extract from
#' @param specials - character vector of text strings to consider "special"
#' @param returnTerms - flag to return a list of terms rather than a formula
#' @param debug - flag to print debugginh info
#' @return a formula, list, or NULL
#'
#' @examplesIf FALSE
#' # example code
#' res = getREs(~1+a:b+(r|x)+s(x));#--FEs+"ordinary" REs+a smooth
#' res = getREs(~a:b);             #--FEs only
#' res = getREs(~ar1(r|x)+s(x));   #--"special" REs plus smooth
#' res = getREs(~s(x));            #--smooth only
#'
#' @importFrom reformulas RHSForm sumTerms
#'
#' @export
#'
getREs <- function(term,
                   specials=names(.valid_covstruct),
                   returnTerms=FALSE,
                   debug=FALSE) {

  ## drop LHS from two-sided formula
  if (length(term) == 3 && identical(term[[1]], quote(`~`))) {
      term <- reformulas::RHSForm(term, as.form = TRUE)
  }

  fbx_terms <- getTermsREs(term,specials=specials,debug=debug);
  if (debug) cat("fbx(term): ", deparse(fbx_terms),"\n\n");
  if (!is.null(fbx_terms)) {
    if (returnTerms) return(fbx_terms);
    fbx_terms = reformulas::sumTerms(c(0,fbx_terms));#--add in zero intercept
    return(makeRHSFormula(fbx_terms))
  }
  return(NULL);
}#--getREs


##' Parse a formula into random effect terms, identifying covariance structures
##'
##' Hacked from [reformulas::splitForm()] to consider only RE terms.
##' reformulas::splitForm taken from Steve Walker's lme4ord,
##' ultimately from the flexLambda branch of [lme4]
##' <https://github.com/stevencarlislewalker/lme4ord/blob/master/R/formulaParsing.R>.
##' @title Split formula containing special random effect terms
##' @param reform - a formula containing random effect terms, or a single RE term
##' @param defaultTerm - cov_struct name for non-special RE terms (default="diag")
##' @param cov_structs - vector of names of valid covariance structures
##' @param debug - flag (TRUE/FALSE) to print debugging info
##' @return a list containing elements
##' \itemize{
##'  \item{\code{reTrms} list with each RE term as an element}
##'  \item{\code{reTrmFormulas} list of \code{x | g} formulas for each term}
##'   \item{\code{reTrmAddArgs} list of function+additional arguments, i.e. \code{list()} (non-special), \code{foo()} (no additional arguments), \code{foo(addArgs)} (additional arguments)}
##'   \item{\code{reTrmClasses} (vector of special functions/classes, as character)}
##' }
##' @examples
##' splitForm_RE(~x+y)                     ## no REs
##' splitForm_RE(~x+y+(f|g))               ## no covariance structure specified
##' splitForm_RE(~x+y+diag(f|g))           ## covariance structure specified for RE term
##' splitForm_RE(~x+y+(diag(f|g)))         ## 'hidden' term with covariance structure specified
##' splitForm_RE(~x+y+(f|g)+cs(1|g))       ## combination of terms with covariance structures un/specified
##' splitForm_RE(~x+y+(1|f/g))                      ## 'slash' term
##' splitForm_RE(~x+y+(1|f/g/h))                    ## 'slash' term
##' splitForm_RE(~x+y+(1|(f/g)/h))                  ## 'slash' term
##' splitForm_RE(~x+y+(f|g)+cs(1|g)+cs(a|b,stuff))  ## complex covariance structures specified
##' splitForm_RE(~(((x+y))))                        ## lots of parentheses
##' splitForm_RE(~1+rr(f|g,n=2))                    ## covariance structure specified with arguments
##' splitForm_RE(~1+s(x, bs = "tp"))                ## smooth special (ignored)
##'
##' @author Steve Walker (hacked by William Stockhausen)
##' @export
splitForm_RE <- function(reform,
                         defaultTerm="diag",
                         cov_structs = findValidCovStructs(),
                         debug=FALSE) {
  #--split formula into separate random effects terms
  ## with covariance structure information ("specials") retained
  bars = reformulas::findbars_x(reform,specials=findValidCovStructs(),default.special="diag",debug=debug);
  formSplits = reformulas::expandAllGrpVar(bars);

  ##--old: fbxx <- reformulas::findbars_x(formula, debug, cov_structs,default.special=defaultTerm)
  ##--old: formSplits <- reformulas::expandAllGrpVar(fbxx);

  if (length(formSplits)>0) {
      formSplitID <- sapply(lapply(formSplits, "[[", 1), as.character)
                                      # warn about terms without a
                                      # setReTrm method

      parenTerm <- formSplitID == "("
                                      # capture additional arguments
      reTrmAddArgs <- lapply(formSplits, "[", -2)[!parenTerm]
                                      # remove these additional
                                      # arguments
      formSplits <- lapply(formSplits, "[", 1:2)
                                      # standard RE terms
      formSplitStan <- formSplits[parenTerm]
                                      # structured RE terms
      formSplitSpec <- formSplits[!parenTerm]

      reTrmFormulas <- c(lapply(formSplitStan, "[[", 2),
                         lapply(formSplitSpec, "[[", 2))
      reTrmFormulas <- unlist(reTrmFormulas) # Fix me:: added for rr structure when it has n = 2, gives a list of list... quick fix
      reTrmClasses <- c(rep(defaultTerm, length(formSplitStan)),
                        sapply(lapply(formSplitSpec, "[[", 1), as.character))
  } else {
      reTrmFormulas <- reTrmAddArgs <- reTrmClasses <- NULL
  }

  names(bars) = vapply(bars,deparse,"");

  return(list(reTrms = bars,
              reTrmFormulas = reTrmFormulas,
              reTrmAddArgs  = reTrmAddArgs,
              reTrmClasses  = reTrmClasses));
}



#' @title Replace part (or all) of a formula term
#' @description Function to replace part (or all) of a formula term.
#' @param term -  the term in which to make the replacement
#' @param target -  the part to replace
#' @param repl - the replacement
#' @result term with part replaced
#' @details Copied from [reformulas:::replaceTerm()].
#' @importFrom reformulas inForm
#' @export
#'
replaceTerm <- function(term,target,repl) {
    if (identical(term,target)) return(repl)
    if (!reformulas::inForm(term,target)) return(term)
    if (length(term) == 2) {
        return(substitute(OP(re_term),list(OP=replaceTerm(term[[1]],target,repl),
                                      re_term=replaceTerm(term[[2]],target,repl))))
    }
    return(substitute(OP(re_term,y),list(OP=replaceTerm(term[[1]],target,repl),
                                         re_term=replaceTerm(term[[2]],target,repl),
                                         y=replaceTerm(term[[3]],target,repl))))
}


splitTerm_RE<-function(term,specials=findValidCovStructs()){
  dp1 = deparse1(term[[1]]);
  if (dp1 %in% specials) term = term[[2]];
  list(covstr=dp1,re=term[[2]],group=term[[3]]);
}

#' @param re_term a language object of the form  effect | groupvar
#' @param model_frame model frame
#' @param drop.unused.levels - flag to drop unused levels from model matrix
#' @param reorder.vars - flag to reorder variables
#' @param sparse (logical) set up sparse model matrices?
#' @return list, see details
#' @details Hacked from [reformulas::mkBlist()]
#'
#' Return list has elements:
#' \itemize{
#'  \item{ff - facotrs}
#'  \item{termZt - transposed Z for term}
#'  \item{ngroups - number of factor groups in term
#'  \item{npars number of RE parameters in term
#'  \item{colnms column names in termZt
#'  \item{re_term - expression for RE term
#'  \item{covstr = covariance function string
#'  \item{factorized_model_frame - model frame with non-numeric columns as factors
#'  \item{Qtp - template for precision matrix
#'  \item{fillQ - function to fill precision matrix
#' }
#'
#' @importFrom Matrix KhatriRao fac2sparse sparse.model.matrix
#' @importFrom stats model.matrix
#'
#' @export
#'
mkBlist <- function(re_term,
                    model_frame,
                    drop.unused.levels=TRUE,
                    reorder.vars=FALSE,
                    sparse = NULL) {
  #--factorize non-numeric elements of model frame
  frloc <- reformulas:::factorize(re_term,model_frame);

  ##--try to evaluate grouping factor within model frame ...
  sp_trm = splitTerm_RE(re_term);
  ff0 <- replaceTerm(sp_trm$group, quote(`:`), quote(`%i%`));
  ff <- try(eval(substitute(makeFac(fac),
                            list(fac = ff0)),
                 frloc), silent = TRUE)
  if (inherits(ff, "try-error")) {
      stop("couldn't evaluate grouping factor ",
           deparse1(sp_trm$group)," within model frame:",
           "error =",
           c(ff),
           " Try adding grouping factor to data ",
           "frame explicitly if possible",call.=FALSE)
  }
  if (all(is.na(ff)))
      stop("Invalid grouping factor specification, ",
           deparse1(sp_trm$group),call.=FALSE)

  if (drop.unused.levels) ff <- factor(ff, exclude=NA)
  ngrps <- length(levels(ff)); #--number of levels in grouping factor

  ##--this section implements eq. 6 of the JSS lmer paper----
  ###--model matrix based on LHS of random effect term (X_i)----
  ###--    sp_trm$re is the LHS (terms) of the a|b formula----
  has.sparse.contrasts <- function(re_term) {
    cc <- attr(re_term, "contrasts")
    !is.null(cc) && is(cc, "sparseMatrix")
  }
  any.sparse.contrasts <- any(vapply(frloc, has.sparse.contrasts, FUN.VALUE = logical(1)))
  mMatrix <- if (!isTRUE(sparse) && !any.sparse.contrasts) model.matrix else Matrix::sparse.model.matrix
  mm <- mMatrix(eval(substitute( ~ foo, list(foo = sp_trm$re))), frloc)
  if (reorder.vars) {
      mm <- mm[colSort(colnames(mm)),]
  }

  ##--this is J^T (see p. 9 of JSS lmer paper)----
  ###--construct indicator matrix for groups by observations
  ###--use fac2sparse() rather than as() to allow *not* dropping
  ###--unused levels where desired
  sm <- Matrix::fac2sparse(ff, to = "d",
                           drop.unused.levels = drop.unused.levels)
  termZt <- Matrix::KhatriRao(sm, t(mm));
  dimnames(termZt) <- list(rep(levels(ff),each=ncol(mm)),
                           rownames(mm));

  ##--determine precision/covariance matrix structure----
  npars = ncol(mm); #--number of random effects parameters/group
  q = ngrps*npars;
  if (sp_trm$covstr %in% c("diag","iid")){
    #--covariance parameters: ln-scale standard deviation
    Qtp = buildTemplateQ.diag(npars,ngrps); #--template for Q
    idx = Qtp@x;
    fillQ = function(covpars){fillQ.ar1(covpars,ngrps,idx,sequential=TRUE);}
  } else
  if (sp_trm$covstr %in% c("ar1")){
    #--covariance parameters: ln-scale standard deviation, correlation coefficient on symlogit scale (-inf,inf)->(-1,1)
    Qtp = buildTemplateQ.ar1(npars,ngrps,sequential=TRUE);
    idx = Qtp@x;
    fillQ = function(covpars){fillQ.ar1(covpars,ngrps,idx,sequential=TRUE);}
  } else {
    stop(paste0("unrecognized covstr '",sp_trm$covstr,"' when creating covariance/precision matrices\n"))
  }

  return(list(ff = ff,
              termZt = termZt,
              ngroups = ngrps,
              npars = npars,
              colnms = colnames(mm),
              re_term = re_term,
              covstr = sp_trm$covstr,
              factorized_model_frame = frloc,
              Qtp = Qtp,
              fillQ = fillQ));
} #--mkBlist
