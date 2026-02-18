#' @title Convert a vector to a factor
#' @description Function to convert a vector to a factor
#' @param v - the vector to convert
#' @param char.only - flag to convert `v` only if it is a character vector
#' @details Copied from `reformulas:::makeFac()` (char.only option changed to TRUE).
#' @export
#'
makeFac <- function(v,char.only=TRUE) {
    if (!is.factor(v) && ((!char.only && is.numeric(v)) || is.character(v))) factor(v) else v
}

#' @title Convert columns in a model frame to factors if included in a formula
#' @description Function to convert columns in a dataframe to factors if included in a formula.
#' @param frmla - formula or term
#' @param model_frame - the dataframe to convert
#' @param char.only - flag to convert only character columns to factors
#' @details Copied from `reformulas:::factorize()` (char.only option changed to TRUE).
#' @export
#'
factorize<-function(frmla, model_frame, char.only = TRUE) {
    for (i in all.vars(reformulas::RHSForm(frmla))) {
        if (!is.null(curf <- model_frame[[i]]))
            model_frame[[i]] <- makeFac(curf, char.only)
    }
    return(model_frame);
}

#' @title infix interaction operator to combine levels from two factors
#' @description Operator to combine levels from two factors.
#' @param f1 - factor
#' @param f2 - factor
#' @param order - flag to re-order levels
#' @return factor with merged levels from `f1` and `f2`
#' @details Copied from reformulas:::`%i%`.
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
##' ultimately from the flexLambda branch of the `lme4` package
##' <https://github.com/stevencarlislewalker/lme4ord/blob/master/R/formulaParsing.R>.
##' @title Split formula containing random effect terms
##' @param re_form - a formula containing random effect terms, or a single RE term
##' @param defaultCovStruct - covariance structure name for non-"special" RE terms (default="diag")
##' @param debug - flag (TRUE/FALSE) to print debugging info
##' @return NULL, or a list containing elements
##' \itemize{
##'  \item{\code{reTrms} named list with each RE term as an element}
##'  \item{\code{reTrmFormulas} named list of \code{x | g} formulas for each term}
##'   \item{\code{reTrmAddArgs} named list (see details)}
##'   \item{\code{reTrmCovTypes} named character vector of covariance structure types}
##' }
##'
##' @details The names associated with all elements in the returned list are the same names as used for the `reTrms` list.
##'
##' Each element in the `reTrmAddArgs` list is a language object, which itself has a list-like structure (use `[[` to access elements).
##' For each language object, the first element will be the covariance structure type, as a `symbol`.
##' Any additional elements represent additional arguments determining the covariance structure.
##' Additional arguments of the form "x=value"can be recovered as `value = el[["x"]]`, where `el` is list element.
##'
##' @examples
##' sf = splitForm_RE(re_form=~x+y)                     ## no REs
##' sf = splitForm_RE(re_form=~x+y+(f|g))               ## no covariance structure specified
##' sf = splitForm_RE(re_form=~x+y+ar1(f|g))            ## covariance structure specified for RE term
##' sf = splitForm_RE(re_form=~x+y+(diag(f|g)))         ## 'hidden' term with covariance structure specified
##' sf = splitForm_RE(re_form=~x+y+((ar1(f|g))))        ## doubly-'hidden' term with covariance structure specified
##' sf = splitForm_RE(re_form=~x+y+(f|g)+cs(1|g))       ## combination of terms with covariance structures un/specified
##' sf = splitForm_RE(re_form=~x+y+(1|f/g))                      ## 'slash' term: equivalent to (1|g:f) + (1|f)
##' sf = splitForm_RE(re_form=~x+y+(1|f/g/h))                    ## 'slash' term: equivalent to (1|h:g:f) + (1|g:f) + (1|f)
##' sf = splitForm_RE(re_form=~x+y+(1|(f/g)/h))                  ## 'slash' term: equivalent to (1|h:g:f) + (1|g:f) + (1|f)
##' sf = splitForm_RE(re_form=~x+y+(f|g)+cs(1|g)+cs(a|b,stuff))  ## complex covariance structures specified
##' sf = splitForm_RE(re_form=~(((x+y))))                        ## lots of parentheses
##' sf = splitForm_RE(re_form=~1+cs(f|g,n=2))                    ## covariance structure with one additional argument
##' sf = splitForm_RE(re_form=~1+cs(f|g,n=2,stuff))              ## covariance structure with multiple additional arguments (TODO: returns NULL!)
##' sf = splitForm_RE(re_form=~1+cs(f|g,n=2,m=3))                ## covariance structure with multiple additional arguments (TODO: returns NULL!)
##' sf = splitForm_RE(re_form=~1+cs(f|g,list(n=2,m=3)))          ## covariance structure with multiple additional arguments given as list (WORKS!)
##' sf = splitForm_RE(re_form=~1+s(x, bs = "tp"))                ## smooth special (ignored)
##'
##' @author Steve Walker (hacked by William Stockhausen)
##' @importFrom reformulas expandAllGrpVar findbars_x
##' @export
splitForm_RE <- function(re_form,
                         defaultCovStruct="diag",
                         debug=FALSE) {

  #--split formula into separate random effects terms
  ##--with covariance structure information ("specials") retained
  ##--`formSplits0` is a list of language elements) representing the RE terms with the
  ##--default covariance structure specified for any non-"special" term.
  ##
  ##--Each language element (accessed as a list) in `formSplits0` has two elements:
  ###--[[1]]  a symbol indicating the covariance structure
  ###--[[2]] a language element representing the "raw" RE term (i.e. (expr|group,additional arguments))
  formSplits0 = reformulas::findbars_x(re_form,specials=findValidCovStructs(),default.special=defaultCovStruct,debug=debug);
  if (is.null(formSplits0)) return(NULL);

  #--`formSplits1` seems to be identical to `formSplits0`. TODO: remove and rename `bars` to `formSplits`?
  formSplits1 = reformulas::expandAllGrpVar(formSplits0);

  if (!is.null(formSplits0)) names(formSplits0) = vapply(formSplits0,deparse,"");
  if (!is.null(formSplits1)) names(formSplits1) = vapply(formSplits1,deparse,"");

  #--identify terms without covariance structure information
  ##--TODO: don't think this can happen because default cov. structure info is applied when
  ##--creating `formSplits0` so no element of `formSplitID` is "(")
  formSplitID <- sapply(lapply(formSplits1, "[[", 1), as.character)
  parenTerm <- formSplitID == "("; #--vector elements should be all FALSE

  # capture additional arguments
  reTrmAddArgs <- lapply(formSplits1, "[", -2)[!parenTerm]

  # remove these additional arguments
  formSplits2 <- lapply(formSplits1, "[", 1:2); # standard RE terms, no additional arguments
  formSplitStan <- formSplits2[parenTerm];      # (standard) RE terms without covariance structure info (not sure this is possible at this point)
  formSplitSpec <- formSplits2[!parenTerm]      # (special) RE terms with covariance structure info

  #--extract "raw" formula for each RE term (no covariance structure info, no additional arguments)
  reTrmFormulas <- c(lapply(formSplitStan, "[[", 2),
                     lapply(formSplitSpec, "[[", 2));
  reTrmFormulas <- unlist(reTrmFormulas); # Fix me:: added for rr structure when it has n = 2, gives a list of list... quick fix

  #--extract character vector of covariance structure types
  reTrmCovTypes <- c(rep(defaultCovStruct, length(formSplitStan)),
                    sapply(lapply(formSplitSpec, "[[", 1), as.character));
  names(reTrmCovTypes) <- c(names(formSplits2)[parenTerm],names(formSplits2)[!parenTerm]);

  return(list(reTrms = formSplits0,
              reTrmFormulas = reTrmFormulas,
              reTrmAddArgs  = reTrmAddArgs,
              reTrmCovTypes = reTrmCovTypes));
}



#' @title Replace part (or all) of a formula term
#' @description Function to replace part (or all) of a formula term.
#' @param term -  the term in which to make the replacement
#' @param target -  the part to replace
#' @param repl - the replacement
#' @return term with part replaced
#' @details Copied from `reformulas:::replaceTerm()`.
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
  list(covstr=dp1,re=term[[2]],group=term[[3]]); #--TODO: want to take account of (or drop) "additional arguments" here!
}

#' @title Make an info list for a RE term
#' @description Function to make an info list for a RE term.
#' @param re_term a language object of the form  (effect | groupvar)
#' @param model_frame model frame
#' @param drop.unused.levels - flag to drop unused levels from model matrix
#' @param sparse (logical) set up sparse model matrices?
#' @return list, see details
#' @details Hacked from `reformulas:::mkBlist()`.
#'
#' Return list has elements:
#' \itemize{
#'  \item{re_term - expression for RE term}
#'  \item{re_term_lbl - character string for RE term}
#'  \item{factorized_model_frame - model frame with non-numeric columns as factors}
#'  \item{ff - factors}
#'  \item{termZt - transposed Z for term}
#'  \item{ngrps - number of factor groups in term}
#'  \item{nprs - number of RE parameters (i.e., per group) in term}
#'  \item{npvs - number of RE parameter values in term}
#'  \item{pr_lbls - labels for parameters}
#'  \item{pv_lbls - labels for parameter values}
#'  \item{covstr - covariance structure string}
#'  \item{Qtp - template for precision matrix}
#'  \item{fillQ - function to fill precision matrix}
#' }
#'
#' @importFrom Matrix KhatriRao fac2sparse sparse.model.matrix
#' @importFrom stats model.matrix
#'
#' @export
#'
mkReTermInfoList <- function(re_term,
                            model_frame,
                            drop.unused.levels=TRUE,
                            sparse = NULL) {
  #--get label for re_term
  re_term_lbl = deparse(re_term);

  #--factorize non-numeric elements of model frame
  fr_mf <- factorize(re_term,model_frame);

  #--try to evaluate grouping factor within model frame ...
  ##--split term to identify grouping factor
  sp_trm = splitTerm_RE(re_term);
  ##--convert grouping factor expression from using ":" to using "%i%
  gf0 <- replaceTerm(sp_trm$group, quote(`:`), quote(`%i%`));
  ##--evaluate grouping factor
  gf <- try(eval(substitute(makeFac(fac),list(fac = gf0)),fr_mf),
            silent = TRUE
           );
  if (inherits(gf, "try-error")) {
      stop("couldn't evaluate grouping factor ",
           deparse1(sp_trm$group)," within model frame:",
           "error =",
           c(gf),
           " Try adding grouping factor to data ",
           "frame explicitly if possible",call.=FALSE)
  }
  if (all(is.na(gf))){
      stop("Invalid grouping factor specification, ",
           deparse1(sp_trm$group),call.=FALSE)
  }

  if (drop.unused.levels) gf <- factor(gf, exclude=NA)
  ngrps <- length(levels(gf)); #--number of levels in grouping factor

  ##--this section implements eq. 6 of the JSS lmer paper----
  ###--model matrix based on LHS of random effect term (X_i)----
  ###--    sp_trm$re is the LHS (terms) of the a|b formula----
  has.sparse.contrasts <- function(re_term) {
    cc <- attr(re_term, "contrasts")
    !is.null(cc) && is(cc, "sparseMatrix")
  }
  any.sparse.contrasts <- any(vapply(fr_mf, has.sparse.contrasts, FUN.VALUE = logical(1)))
  mMatrix <- if (!isTRUE(sparse) && !any.sparse.contrasts) model.matrix else Matrix::sparse.model.matrix
  mmX_i <- mMatrix(eval(substitute( ~ foo, list(foo = sp_trm$re))), fr_mf); #--"raw" RE model matrix X_i in Table 3, JSS lmer paper
  nprs    = ncol(mmX_i);    #--number of RE parameters/group
  pr_lbls = colnames(mmX_i);#--parameter labels
  #cat("nprs:",nprs,"\npr_lbls:",pr_lbls,"\n");

  ##--construct transpose of indicator matrix J_i in Table 3, JSS lmer paper
  ###--with dimensions grouping factor levels (ngrps) by observations (nobs).
  ####--use fac2sparse() rather than as() to allow *not* dropping unused levels where desired.
  imJt_i <- Matrix::fac2sparse(gf, to = "d",
                               drop.unused.levels = drop.unused.levels);
  ##--construct transpose of term-wise random effects model matrix Z_i in Table 3, JSS lmer paper
  termZt <- Matrix::KhatriRao(imJt_i, t(mmX_i));
  dimnames(termZt) <- list(rep(levels(gf),each=nprs),
                           rownames(mmX_i));
  pv_lbls = paste(re_term_lbl,
                  paste(pr_lbls,rep(levels(gf),each=nprs),sep="@"),
                  sep="@"); #--cparameter value labels: re_term @ parameter label @-grouping factor value
  #cat("rownames(termZt): ","\n",paste0(rownames(termZt),collase="\n"),"\n");
  #cat("pv_lbls:","\n",paste0(pv_lbls,collase="\n"),"\n");
  rownames(termZt) = pv_lbls;
  npvs = length(pv_lbls);

  ##--determine precision/covariance matrix structure----
  q = ngrps*nprs;  #--total number of RE parameters
  if (sp_trm$covstr %in% c("diag","iid")){
    #--covariance parameters: ln-scale standard deviation
    Qtp = buildTemplateQ.diag(nprs,ngrps); #--template for Q
    idx = Qtp@x;
    fillQ = function(covpars){fillQ.ar1(covpars,ngrps,idx,sequential=TRUE);}
  } else
  if (sp_trm$covstr %in% c("ar1")){
    #--covariance parameters: ln-scale standard deviation, correlation coefficient on symlogit scale (-inf,inf)->(-1,1)
    Qtp = buildTemplateQ.ar1(nprs,ngrps,sequential=TRUE);
    idx = Qtp@x;
    fillQ = function(covpars){fillQ.ar1(covpars,ngrps,idx,sequential=TRUE);}
  } else
  if (sp_trm$covstr %in% c("us")){
    #--covariance parameters: ln-scale standard deviation, correlation coefficient on symlogit scale (-inf,inf)->(-1,1)
    Qtp = buildTemplateQ.us(nprs,ngrps,sequential=TRUE);
    idx = Qtp@x;
    fillQ = function(covpars){fillQ.ar1(covpars,ngrps,idx,sequential=TRUE);}
  } else {
    stop(paste0("unrecognized covstr '",sp_trm$covstr,"' when creating covariance/precision matrices\n"))
  }

  return(list(re_term = re_term,
              factorized_model_frame = fr_mf,
              grouping_factor = gf,
              termZt = termZt,
              re_term_lbl = re_term_lbl,
              ngrps = ngrps,          #--number of groups
              grp_lbls = levels(gf),  #--group labels
              nprs = nprs,            #--number opf parameters
              pr_lbls = pr_lbls,      #--parameter labels
              npvs = npvs,            #--number of parameter values
              pv_lbls = pv_lbls,      #--parameter value labels
              covstr = sp_trm$covstr,
              Qtp = Qtp,
              fillQ = fillQ));
} #--mkBlist
