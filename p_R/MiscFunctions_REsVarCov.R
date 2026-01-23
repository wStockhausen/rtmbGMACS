#'
#' @title Create the model frame from a parameter formula and base dataframe
#' @description Function to create the model frame from a parameter formula and base dataframe.
#' @param formula - RHS formula describing model for which to create the model frame
#' @param data - base dataframe from which to extract model frame
#' @return A data.frame
#' @details The formula is processed to substitute '+' for RE grouping symbols ('|' and '||')
#' and drop "specials" (see [.valid_covstruct] and [.valid_smooths]) from the formula.
#' Column names for derived columns (e.g., `cos(x)`) are a back-quoted version of the deriving formula.
#' @importFrom reformulas inForm noSpecials sub_specials
#' @export
#'
createModelFrame <- function(formula,data){

    mf_call <- match.call();

    if (reformulas::inForm(formula, quote(`$`))) {
        warning("use of the ", sQuote("$"), " operator in formulas is not recommended");
    }

    environment(formula) <- parent.frame();

    ## now work on evaluating model frame
    mf_call$drop.unused.levels <- TRUE;
    mf_call[[1]] <- as.name("model.frame");
    mf_call$data <- data ## propagate ..offset modification?

    ## want the model frame to contain the union of all variables used in any of the terms
    ## combine all formulas
    f <- reformulas::noSpecials(reformulas::sub_specials(formula),
                                delete=FALSE,
                                specials = c(.valid_covstruct, .valid_smooths))
    environment(f) <- environment(formula);

    mf_call$formula <- f;
    mf <- eval(mf_call,envir=environment(f),enclos=parent.frame());

    return(mf);
}

#' @title Create a RHS formula from a call or a (list of) language object(s)
#' @param terms - a call or a (list of) language object(s)
#' @return a RHS formula
#' @examplesIf FALSE
#' # example code
#' f = ~ 1 + x + s(x, bs = "tp") + s(y, bs = "gp") + s(z)
#' bars = reformulas::findbars_x(f,target="s",default.special=NULL);
#' makeRHSFormula(bars);
#'
makeRHSFormula<-function(terms){
  if (is.list(terms))
    terms = reformulas::sumTerms(terms);
  return(as.formula(reformulas::makeOp(terms,as.symbol("~"))))
}

hasREs<-function(formula){
  return(!is.null(reformulas::findbars_x(formula)));
}

#' @title Extract RE terms from a formula
#' @description Function to extract RE terms from a formula
#' @param term - formula or term to extract from
#' @param specials - character vector of text strings to consider "special"
#' @param debug - flag to print debugginh info
#' @return a call (TODO: want RHS formula with only RE terms??, not a "call")
#'
getREs <- function(term,
                   specials=names(.valid_covstruct),
                   debug=TRUE) {

    ## drop LHS from two-sided formula
    if (length(term) == 3 && identical(term[[1]], quote(`~`))) {
        term <- reformulas::RHSForm(term, as.form = TRUE)
    }

    ## base function
    ## defining internally in this way makes debugging slightly
    ## harder, but (1) allows easy propagation of the top-level
    ## arguments down the recursive chain; (2) allows the top-level
    ## expandAllGrpVar() operation (which also handles cases where
    ## a naked term rather than a list is returned)

    fbx <- function(term,pass=0,termno=0) {
      bar = '|';
      dbl = '||';
      if (debug) {
        cat("pass:",pass," termno:",termno,"\n");
        cat(deparse(term),"\n");
      }
        if (is.name(term) || !is.language(term)) {
          if (debug) cat("first check: returning NULL\n\n")
          return(NULL)
        }
        if (list(term[[1]]) %in% lapply(specials,as.name)) {
            if (debug) cat("special: ",deparse(term),"\n")
            return(term)
        }
        if (term[[1]] == as.name(bar)) {  ## found x | g
            if (debug) cat(sprintf("%s term: %s\n", "bar", deparse(term)))
            return(term)
        }
        if (term[[1]] == as.name(dbl)) {
            if (debug) cat(sprintf("%s term: %s\n", "double-bar", deparse(term)))
            return(term)
        }
        if (term[[1]] == as.name("(")) {  ## found (...)
            if (debug) cat("paren term:",deparse(term),"\n")
            return(fbx(term[[2]],pass=pass+1,termno=2))  #--recursion here
        }
        stopifnot(is.call(term))
        if (length(term) == 2) {
            ## unary operator, decompose argument
            if (debug) cat("unary operator:",deparse(term[[1]]),",",deparse(term[[2]]),"\n")
            return(fbx(term[[2]],pass=pass+1,termno=2))  #--recursion here
        }
        ## binary operator, decompose both arguments
        f2 <- fbx(term[[2]],pass=pass+1,termno=2) #--recursion
        f3 <- fbx(term[[3]],pass=pass+1,termno=3) #--recursion

        if (debug) { cat("binary operator:",
                         deparse(term[[1]]),",",
                         deparse(term[[2]]),",",
                         deparse(term[[3]]),"\n")
                         cat("term 2: ", deparse(f2), "\n")
                         cat("term 3: ", deparse(f3), "\n")
                         cat("\n")
        }
        c(f2, f3)
    }#--fbx

    fbx_terms <- fbx(term);
    if (debug) cat("fbx(term): ", deparse(fbx_terms),"\n\n");
    if (!is.null(fbx_terms)) fbx_terms = reformulas::sumTerms(c(0,fbx_terms));#--add in zero intercept
    return(fbx_terms);
}#--getREs

hasSmooths<-function(formula){
  return(reformulas::anySpecial(formula, specials=.valid_smooths));
}

#' @title Extract smooth terms from a formula
#' @param term - formula or term from which to extract smooth terms
#' @return list of smooth terms
#' @details Valid smooths are found in rtmbGMACS:::.valid_smooths.
#' Based on [reformulas::findbars_x()].
#'
#' @examplesIf FALSE
#' # example code
#' f = ~ 1 + x + s(x, bs = "tp") + s(y, bs = "gp") + s(z) + (g|p)
#' sms = getSmoothTerms(f);
#' @importFrom reformulas makeOp
#' @export
#'
getSmoothTerms <- function(term,
                            targets = .valid_smooths,
                            debug=FALSE) {

  specials=character(0);
  default.special=NULL;

    ## drop RHS from two-sided formula
    if (length(term) == 3 && identical(term[[1]], quote(`~`))) {
        term <- RHSForm(term, as.form = TRUE)
    }
    ds <- NULL;

    ## base function
    ## defining internally in this way makes debugging slightly
    ## harder, but allows easy propagation of the top-level
    ## arguments down the recursive chain

    fbx <- function(term,target) {
        if (is.name(term) || !is.language(term)) return(NULL)
        if (list(term[[1]]) %in% lapply(specials,as.name)) {
            if (debug) cat("special: ",deparse(term),"\n")
            return(term)
        }
        if (head(term) == as.name(target)) {  ## found s(x)
            if (debug) {
                tt <- sprintf('"%s"', target)
                cat(sprintf("%s term: %s\n", tt, deparse(term)))
            }
            if (is.null(ds)) return(term)
            return(makeOp(term, ds))
        }
        # if (head(term) == as.name("||")) {
        #     if (expand_doublevert_method == "diag_special") {
        #         return(makeOp(makeOp(term[[2]], term[[3]],
        #                              op = quote(`|`)),
        #                       as.name("diag")))
        #     }
        #     ## expand_double_vert_method == "split" if we get here
        #     return(lapply(expandDoubleVert(term), fbx))
        # }
        if (head(term) == as.name("(")) {  ## found (...)
            if (debug) cat("paren term:",deparse(term),"\n")
            return(fbx(term[[2]]))
        }
        stopifnot(is.call(term))
        if (length(term) == 2) {
            ## unary operator, decompose argument
            if (debug) cat("unary operator, target ",target,":",deparse(term[[2]]),"\n")
            return(fbx(term[[2]],target));
        }
        ## binary operator, decompose both arguments
        f2 <- fbx(term[[2]],target=target);
        f3 <- fbx(term[[3]],target=target);

        if (debug) { cat("binary operator, target ",target,":",deparse(term[[2]]),",",
                         deparse(term[[3]]),"\n")
                         cat("term 2: ", deparse(f2), "\n")
                         cat("term 3: ", deparse(f3), "\n")
        }
        c(f2, f3)
    }

    fbx_terms <- list(); k=0;
    for (i in seq_along(targets)){
      if (debug) cat("checking target",i,"'",targets[i],"'\n")
      fbx_res <- fbx(term,targets[i]); j=0;
      for (j in seq_along(fbx_res)) {
        fbx_terms[k+j] = fbx_res[j];
        if (debug) {
          if (!is.null(fbx_terms[k+j])) {
            cat("fbx(term,",targets[i],"): ", deparse(fbx_terms[k+j]),"\n\n");
          } else {
            cat("fbx(term,",targets[k+j],"): is NULL\n\n")
          }
        }
      }
      k = k+j;
    }
    #expandAllGrpVar(fbx_term)
  return(fbx_terms);
}

#' @title: Create "better" names for fixed effects parameters.
#' @return dataframe with original (`orig`) and revised (`revd`) names.
#' @details
#' The original parameter names (given by `colnames(X)`)for those involving factors are combinations of
#' the factor names and the factor levels, with interactions indicated by a colon (":").
#' The revised names for these drop the factor names and split interactions with a "_").
#'
#' @importFrom stringr str_split_1 str_remove
#' @import from tibble tibble
#' @export
#'
createParamNamesForFEs<-function(fe_form,X,debug=TRUE){
  trms = terms(fe_form);
  vars= rownames(attr(trms,"factors"));
  rgx = paste0("^",vars,collapse="|");
  dfrNames <- tibble::tibble(orig=colnames(X),revd="");
  extract <- function(o){
      if (stringr::str_detect(o,"\\(")){
        r = o;
      } else {
        r = stringr::str_remove(o,rgx);
      }
    return(r);
  }
  for (i in 1:nrow(dfrNames)){
    o = stringr::str_split_1(dfrNames$orig[i],stringr::fixed(":"));
    if (debug) cat("o: '",paste(o,collapse=" "),"'\n",sep="");
    if (length(o)==1) {
      dfrNames$revd[i] = extract(o);
    } else {
      r = "";
      for (j in 1:(length(o)-1)){
        r = paste0(r,extract(o[j]),"*");
      }
      j = length(o);
      r = paste0(r,extract(o[j]));
      dfrNames$revd[i] = r;
    }
  }
  return(dfrNames);
}

createParamsInfo<-function(formula, fr, ranOK=TRUE, type="",
                           contrasts = NULL, sparse=FALSE, old_smooths = NULL){
  ##--keep only non-smooths FIXED EFFECTS (remove random effects and smooths) from formula and process using getReXtrms
  fe_form <- formula;
  reformulas::RHSForm(fe_form) <- reformulas::nobars(reformulas::noSpecials(reformulas::RHSForm(reformulas:::noSpecials(fe_form,specials="s"))));
  res_FEs <- getXReTrms(fe_form, fr, ranOK=ranOK, type=type,
                         contrasts = contrasts, sparse=sparse, old_smooths = old_smooths);
  attr(res_FEs$X,"term_names")<-createParamNamesForFEs(fe_form,res_FEs$X);
  View(attr(res_FEs$X,"term_names"));
  ###--only need to keep X, terms, and offset? just X?

  ##--keep only RANDOM EFFECTS (remove fixed effects and smooths) and process using getReXtrms
  has_re <- hasREs(formula);
  if (has_re){
    re_form <- formula;
    reformulas::RHSForm(re_form) <- getREs(re_form);
    res_REs <- getXReTrms(re_form, fr, ranOK=ranOK, type=type,
                           contrasts = contrasts, sparse=sparse, old_smooths = old_smooths);
  }

  ##--keep only SMOOTHS (remove fixed and random effects) and process using getReXtrms
  has_sms <- hasSmooths(formula);
  if (has_sms){
    sm_form <- formula;
    reformulas::RHSForm(sm_form) <- reformulas::RHSForm(makeRHSFormula(getSmoothTerms(sm_form)));
    res_SMs <- getXReTrms(sm_form, fr, ranOK=ranOK, type=type,
                           contrasts = contrasts, sparse=sparse, old_smooths = old_smooths);
  }

  ##--TODO: decide what to return (combine info from above results?)
  return(NULL);#--fill in
}

##--Based on glmmTMB: create FE, smooth, and RE terms from formula----
##' @param formula current formula, containing both fixed & random effects
##' @param fr full model frame
##' @param ranOK random effects allowed here?
##' @param type label for model type
##' @param contrasts a list of contrasts (see ?glmmTMB)
##' @param sparse (logical) return sparse model matrix?
##' @param old_smooths smooth information from a prior model fit (for prediction)
##' @return a list composed of
##' \item{X}{design matrix for fixed effects}
##' \item{Z}{design matrix for random effects}
##' \item{reTrms}{output from \code{\link[reformulas]{mkReTrms}}, possibly augmented with information about \code{mgcv}-style smooth terms}
##' \item{ss}{splitform of the formula}
##' \item{aa}{additional arguments, used to obtain rank}
##' \item{terms}{terms for the fixed effects}
##' \item{offset}{offset vector, or vector of zeros if offset not specified}
##' \item{reXterms}{terms for the model matrix in each RE term}
##' \item{formula}{input formula}
##'
##' @importFrom stats model.matrix contrasts
##' @importFrom methods new
##' @importFrom mgcv smoothCon smooth2random s PredictMat
##' @importFrom reformulas inForm findbars nobars noSpecials sub_specials addForm findbars_x anySpecial RHSForm RHSForm<- extractForm reOnly no_specials splitForm addForm0 makeOp
##' @importFrom utils head
getXReTrms <- function(formula, fr, ranOK=TRUE, type="",
                       contrasts = NULL, sparse=FALSE, old_smooths = NULL) {

    #--determine if model has REs ----
    has_re <- !is.null(reformulas::findbars_x(formula));
    #--determine if model has smooths ----
    has_smooths <- hasSmooths(formula);

    #--determine fixed-effects model matrix X ----
    ## remove random effect parts and smooths from formula:
    fixedform <- formula;
    reformulas::RHSForm(fixedform) <- reformulas::nobars(reformulas::RHSForm(fixedform));

    terms <- NULL ## make sure it's empty in case we don't set it

    nobs <- nrow(fr)

    ## check for empty fixed form
    ## need to ignore environments when checking!
    ##  ignore.environment= arg only works with closures
    idfun <- function(x,y) {
        environment(x) <- emptyenv()
        environment(y) <- emptyenv()
        return(identical(x,y))
    }

    if (idfun(reformulas::RHSForm(fixedform, as.form=TRUE), ~ 0) ||
        idfun(reformulas::RHSForm(fixedform, as.form=TRUE), ~ -1)) {
        ##--no FEs----
        X <- matrix(ncol=0, nrow=nobs); #--FE model matrix has no columns and nobs rows
        offset <- rep(0,nobs);          #--offsets are 0
    } else {
        ##--has FEs:determine model matrix
        ## check for mgcv-style smooth terms, adjust accordingly ...
        if (has_smooths) {
            ###--extract s() terms----
            smooth_terms <- reformulas::findbars_x(fixedform, default.special = NULL, target = "s")
            ###--*remove* s() terms from fixed formula----
            fixedform <- reformulas::noSpecials(fixedform, specials = "s");

            ## FIXME: could be fragile about eval environments (how
            ##  far up do we have to go with eval.parent? Or do we
            ##  use the environment of the formula?
            if (!is.null(old_smooths)) {
                ## we are predicting, want to use old smooths rather than constructing new ones
                smooth_terms2 <- lapply(old_smooths[lengths(old_smooths)>0],
                                        function(s) {
                                            if (is.null(s)) return(NULL)
                                            X <- mgcv::PredictMat(s$sm, fr)   ## get prediction matrix for new data
                                            ## transform to RE parameterization
                                            if (!is.null(s$re$trans.U)) X <- X%*%s$re$trans.U
                                            X <- t(t(X)*s$re$trans.D)
                                            ## re-order columns according to random effect re-ordering...
                                            X[,s$re$rind] <- X[,s$re$pen.ind!=0]
                                            ## re-order penalization index in same way
                                            pen.ind <- s$re$pen.ind; s$pen.ind[s$re$rind] <- pen.ind[pen.ind>0]
                                            ## start return object...
                                            s_new <- list(re = list(rand=list(), Xf=X[,which(s$re$pen.ind==0),drop=FALSE]))
                                            for (i in 1:length(s$re$rand)) { ## loop over random effect matrices
                                                s_new$re$rand[[i]] <- X[, which(pen.ind==i), drop=FALSE]
                                                attr(s_new$re$rand[[i]], "s.label") <- attr(s$re$rand[[i]], "s.label")
                                            }
                                            names(s_new$re$rand) <- names(s$re$rand)
                                            return(s_new)
                                        })
            } else {
                ## new smooths
                smooth_terms2 <- lapply(smooth_terms,
                    function(tt) {
                        ## ‘smoothCon’ returns a list of smooths because factor ‘by’
                        ## variables result in multiple copies of a smooth, each multiplied
                        ## by the dummy variable associated with one factor level.
                        sm <- eval(bquote(mgcv::smoothCon(.(tt),
                                                          ## don't want intercept etc ...
                                                          absorb.cons = TRUE,
                                                          data=fr)))
                        if (length(sm)>1) stop("can't handle 'by' arguments in smooths yet")
                        sm <- sm[[1]]
                        list(sm = sm,
                             ## FIXME: will vnames ("a vector of names
                             ## to avoid as dummy variable names in the
                             ## random effects form") ever be non-empty?
                             re = mgcv::smooth2random(sm, vnames = "", type = 2))
                    })
            } ## create (new) smooth terms
        } ## has_smooths

        ##--determine FE model matrix----
        ###--if model has smooths, then each smooth is represented by 1 column
        tt <- terms(fixedform)
        ###--following commented out from glmmTMB.R
        # pv <- attr(mf$formula,"predvars")
        # attr(tt, "predvars") <- fix_predvars(pv,tt)
        # mf$formula <- tt
        # terms_fixed <- terms(eval(mf,envir=environment(fixedform)))
        # terms <- list(fixed=terms(terms_fixed))
        terms <- list(fixed=terms(tt));
        if (!sparse) {
            X <- model.matrix(reformulas::noSpecials(fixedform, specials = "offset"),
                                                     fr, contrasts);
        } else {
            X <- Matrix::sparse.model.matrix(reformulas::noSpecials(fixedform, specials = "offset"),
                                                                    fr, contrasts);
            ## FIXME? ?sparse.model.matrix recommends MatrixModels::model.Matrix(*,sparse=TRUE)
            ##  (but we may not need it, and would add another dependency etc.)
        }
        if (has_smooths) {
            if (sparse) warning("smooth terms may not be compatible with sparse X matrices")
            for (si in seq_along(smooth_terms2)) {
                cnm <- colnames(X)  ## need to update cnm after each added term ...
                s <- smooth_terms2[[si]]
                if (ncol(s$re$Xf) == 0) next; #--no columns in smooth, so skip
                ## indices within the beta vector corresponding to this smooth
                beta_ind <- ncol(X) + seq(ncol(s$re$Xf))
                snm <- attr(s$re$rand$Xr, "s.label");#--get label for smooth
                X <- cbind(X, s$re$Xf);              #--append smooth columns to X
                colnames(X) <- c(cnm, paste0(snm, seq.int(ncol(s$re$Xf)))); #--name smooth columns in X
                ## store column indices for smooth in X in smooth object
                smooth_terms2[[si]]$re$beta_ind <- beta_ind;
            }
        }

        ## will be 0-column matrix if fixed formula is empty
        offset <- rep(0,nobs);
        if (reformulas::inForm(fixedform,quote(offset))) {
            ## hate to match offset terms with model frame names
            ##  via deparse, but since that what was presumably done
            ##  internally to get the model frame names in the first place ...
            for (o in reformulas::extractForm(fixedform,quote(offset))) {
                offset_nm <- deparse1(o);
                ## don't think this will happen, but ...
                if (length(offset_nm)>1) {
                    stop("trouble reconstructing offset name")
                }
                offset <- offset + fr[[offset_nm]];
            }
        }
    }#--RHS has FEs (possibly including smooths, which are represented as single columns)

    ## ran-effects model frame (for predvars)
    ## important to COPY formula (and its environment)?
    ranform <- formula

    if (!has_re && !has_smooths) {
        reTrms <- reXterms <- NULL
        Z <- new("dgCMatrix",Dim=c(as.integer(nobs),0L)) ## matrix(0, ncol=0, nrow=nobs)
        aa <- integer(0) #added for rr to get rank
        ss <- integer(0)
    } else {
        ##--set up SMOOTHs and REs----
        ## FIXME: check whether predvars are carried along correctly in terms
        if (!ranOK) stop("no random effects allowed in ", type, " term")
        ## FIXME: could use doublevert_split = FALSE here to preserve
        ##  || -> diag() behaviour if we wanted (with new reformulas version)
        ##  WTS: also, converts ar1(...) (at least) to us(...)
        reformulas::RHSForm(ranform) <- reformulas::subbars(reformulas::RHSForm(reformulas::reOnly(formula)));

        if (has_re) {
#            mf$formula <- ranform
            ###--make reTrms----
            ## no_specials so that mkReTrms can handle it
            reTrms <- mkReTrms(reformulas::no_specials(reformulas::findbars_x(formula)),
                               fr,
                               reorder.terms=FALSE,
                               calc.lambdat=TRUE,   #--WTS: changed to TRUE
                               sparse = TRUE)
        } else {
            ## dummy elements
            reTrms <- list(Ztlist = list(), flist = list(), cnms = list(), theta = list(),
                           Zt=list(), Lind=list(),Gp=list(),lower=list(),Lambdat=list(),nl=list(),ord=list()); #--WTS: added these
        } #--if (has_re)

        ##--split `formula` into `fixedFormula` and reTrmFormulas, reTrmAddArgs lists----
        ###--elements in reTrmFormulas are formulas for RE or smooth terms as "effect|group" or "variable"
        sp_form <- reformulas::splitForm(formula, specials = c(names(.valid_covstruct), "s"))

        ##--need to fill in SMOOTH TERMS as REs in **correct locations**----
        ##--post-process mkReTrms to add smooths (incorporate in mkReTrms?)----
        if (has_smooths) {
            ns <- length(sp_form$reTrmClasses);
            augReTrms <- list(Ztlist = vector("list", ns),
                              flist = vector("list", ns),
                              cnms = vector("list", ns),
                              smooth_info = vector("list", ns));
            barpos    <- which(sp_form$reTrmClasses != "s")
            nonbarpos <- which(sp_form$reTrmClasses == "s")
            for (p in c("Ztlist", "flist", "cnms")) {                     #--TODO: fill out missing items like Lambdat, Lind, theta?
                ## fill in values from traditional (bar-containing) REs
                augReTrms[[p]][barpos] <- reTrms[[p]]
                names(augReTrms[[p]])[barpos] <- names(reTrms[[p]])
            }
            ## This code *would* set one theta value for each smooth      #--TODO (WTS): probably want to do this
            ## glmmTMB doesn't use the theta component anyway --
            ## and these computations break because mkReTrms doesn't know
            ## about specials, assumes ntheta = (nlev*(nlev+1)/2)
            ## for each smooth
            ##
            ## augReTrms$theta <- rep(0, ## default value
            ##                        sum(lengths(reTrms$theta)) +
            ##                        length(nonbarpos))
            ## augReTrms$theta[barpos] <- reTrms$theta

            ## only need one 'dummy' factor for all the smooth terms
            ff <- factor(rep(1, nobs))
            augReTrms$flist <- c(reTrms$flist, list(dummy = ff))
            avec <- rep(NA_integer_, ns)
            avec[barpos] <- attr(reTrms$flist, "assign")
            ## mkReTrms returns more than we need (some is for lme4)
            ##  ... which bits are actually used hereafter?
            avec[nonbarpos] <-  length(augReTrms$flist)
            attr(augReTrms$flist, "assign") <- avec
            ncol_fun <- function(x) if (is.null(x)) 0 else ncol(x)
            b_lens <- vapply(augReTrms$Ztlist, ncol_fun, FUN.VALUE = numeric(1))
            for (i in seq_along(smooth_terms2)) {
                s <- smooth_terms2[[i]]
                pos <- nonbarpos[i]
                Zt <- as(t(s$re$rand$Xr), "dgCMatrix")
                b_lens[pos] <- ncol(Zt)
                b_ind <- sum(b_lens[seq_along(b_lens)<i]) + seq(ncol(Zt))
                ###--store b indices----
                smooth_terms2[[i]]$re$b_ind <- b_ind
                npar <- nrow(Zt)
                ###--store Zt----
                augReTrms$Ztlist[[pos]] <- Zt
                nm <- attr(s$re$rand$Xr, "s.label")
                names(augReTrms$Ztlist)[pos] <- nm
                ###--add to cnms----
                augReTrms$cnms[[pos]] <- paste0("dummy", seq(npar))
                names(augReTrms$cnms)[pos] <- "dummy"
            }
            ##--store smooth info in relevant spots----
            for (i in seq_along(nonbarpos)) {
                augReTrms$smooth_info[[nonbarpos[i]]] <- smooth_terms2[[i]]
            }
            ##--reconstitute Zt, Gp----
            augReTrms$Zt <- do.call(rbind, augReTrms$Zt)
            augReTrms$Gp <- cumsum(c(0, vapply(augReTrms$Ztlist, nrow, 0L)))

            ##--replace `reTrms` with `augReTrms`----
            reTrms <- augReTrms;   #--copy smooth-related terms to reTrms
        }#--has_smooths

        ##--create additional arguments (`aa`) list
        sp_form$reTrmClasses[sp_form$reTrmClasses == "s"] <- "homdiag"
        # FIX ME: migrate this (or something like it) down to reTrms,
        ##    allow for more different covstruct types that have additional arguments
        ##  e.g. phylo(.,tree); fixed(.,Sigma)
        # FIX ME: use NA rather than 0 as a placeholder in aa?
        ## FIXME: make sure that eval() happens in the right environment/
        ##    document potential issues
        ## Changed from getting rank to extracting additional argument for propto
        get_arg <- function(v) {
          if (length(v) == 1) return(NA_real_)
          payload <- v[[2]]
          ## rabbit-hole alert. Try to evaluate payload first in model frame,
          ##  then in formula environment (... then in parent env of
          ##  formula env ... ??)
          res <- tryCatch(eval(payload, envir = fr,
                               enclos = environment(formula)),
                          error = function(e)
                            stop("can't evaluate argument ",
                                 sQuote(deparse(payload)),
                                 call. = FALSE))
          return(res)
        }
        aa <- lapply(sp_form$reTrmAddArgs, get_arg); #--additional arguments

        ##--create reXterms list----
        ## terms for the model matrix in each RE term
        ## this is imperfect: it should really be done in mkReTrms/mkBlist,
        ## where we are generating these terms anyway on the way
        ## to constructing Z, but that's in lme4 so we can't change it
        ## unless absolutely necessary
        termsfun <- function(x) {
            ## this is a little magic: copying lme4:::mkBlist approach
            ff <- eval(substitute( ~ foo, list(foo = x[[2]]))) ## make formula from LHS
            tt <- try(terms(ff, data=fr), silent=TRUE)         ## construct terms
            if (inherits(tt,"try-error")) {
                stop(
                    sprintf("can't evaluate RE term %s: simplify?",
                            sQuote(deparse(ff)))
                )
            }
            tt
        }
        ## HACK: should duplicate 'homdiag' definition, keep it as 's' (or call it 'mgcv_smooth")
        ##  so we can recognize it.
        ## Here, we're using the fact that the ...AddArgs stuff is still in an unevaluated form
        drop_s <- function(f, a) {
            if (identical(a[[1]], as.symbol('s'))) NA else termsfun(f)
        }
        reXterms <- Map(drop_s, sp_form$reTrmFormulas, sp_form$reTrmAddArgs); #--create reXterms list

        for (i in seq_along(sp_form$reTrmAddArgs)) {
          if(sp_form$reTrmClasses[i] == "rr") {
            if (!is.na(aa[i]) & is.na(suppressWarnings(as.numeric( aa[i] )))) {
              stop("non-numeric value for reduced-rank dimension", call. = FALSE);
            }
          }
          else if(sp_form$reTrmClasses[i] == "propto"){
            checkProptoNames(aa = aa[[i]], cnms = reTrms$cnms[[i]], reXtrm = reXterms[[i]]);
          }
        }

        ##--save RE term class info as vector----
        ss <- unlist(sp_form$reTrmClasses);

        ##--save transpose of Zt
        Z <- t(reTrms$Zt);   ## still sparse ...
    }#--has RE and/or smooth terms

    ## if(is.null(rankX.chk <- control[["check.rankX"]]))
    ## rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    ## X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    ## if(is.null(scaleX.chk <- control[["check.scaleX"]]))
    ##     scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    ## X <- checkScaleX(X, kind=scaleX.chk)

    ## list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
    ##      wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))

    return(namedList(X, Z, reTrms, ss, aa, terms, offset, reXterms, formula));
} #--getXReTerms

###--From glmmTMB: utils.R----
## generate a list with names equal to values
## See also: \code{tibble::lst}, \code{Hmisc::llist}
namedList <- function (...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L)))
        nm <- snm
    if (any(nonames <- nm == ""))
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}

## reassign predvars to have term vars in the right order,
##  but with 'predvars' values inserted where appropriate
fix_predvars <- function(pv,tt) {
    if (length(tt)==3) {
        ## convert two-sided to one-sided formula
        tt <- RHSForm(tt, as.form=TRUE)
    }
    ## ugh, deparsing again ...
    tt_vars <- vapply(attr(tt, "variables"), deparse1, character(1))[-1]
    ## remove terminal paren - e.g. match term poly(x, 2) to
    ##   predvar poly(x, 2, <stuff>)
    ## beginning of string, including open-paren, colon
    ##  but not *first* comma nor arg ...
    ##  could possibly try init_regexp <- "^([^,]+).*" ?
    init_regexp <- "^([(^:_.[:alnum:]]+).*"
    tt_vars_short <- gsub(init_regexp,"\\1",tt_vars)
    if (is.null(pv) || length(tt_vars)==0) return(NULL)
    new_pv <- quote(list())
    ## maybe multiple variables per pv term ... [-1] ignores head
    ## FIXME: test for really long predvar strings ????
    pv_strings <- vapply(pv,deparse1,FUN.VALUE=character(1))[-1]
    pv_strings <- gsub(init_regexp,"\\1",pv_strings)
    for (i in seq_along(tt_vars)) {
        w <- match(tt_vars_short[[i]],pv_strings)
        if (!is.na(w)) {
            new_pv[[i+1]] <- pv[[w+1]]
        } else {
            ## insert symbol from term vars
            new_pv[[i+1]] <- as.symbol(tt_vars[[i]])
        }
    }
    return(new_pv)
}



##--From glmmTMB::glmmTMB.R: Calculate random effect structure----
###--NOT SURE THE FOLLOWING IS USEFUL!!----
##' Calculates number of random effects, number of parameters,
##' block size and number of blocks.  Mostly for internal use.
##' @param reTrms random-effects terms list
##' @param ss a vector of character strings indicating a valid covariance structure (one for each RE term).
##' Must be one of \code{names(glmmTMB:::.valid_covstruct)};
##' default is to use an unstructured  variance-covariance
##' matrix (\code{"us"}) for all blocks).
##' @param reXterms terms objects corresponding to each RE term
##' @param fr model frame
##' @param aa additional arguments (i.e. rank, or var-cov matrix)
##' @inheritParams glmmTMBControl
##' @return a list
##' \item{blockNumTheta}{number of variance covariance parameters per term}
##' \item{blockSize}{size (dimension) of one block}
##' \item{blockReps}{number of times the blocks are repeated (levels)}
##' \item{covCode}{structure code}
##' \item{simCode}{simulation code; should we "zero" (set to zero/ignore), "fix" (set to existing parameter values), "random" (draw new random deviations)?}
##' \item{fullCor}{logical vector (compute/store full correlation matrix?)}
##' @examples
##' data(sleepstudy, package="lme4")
##' rt <- lme4::lFormula(Reaction~Days+(1|Subject)+(0+Days|Subject),
##'                     sleepstudy)$reTrms
##' rt2 <- lme4::lFormula(Reaction~Days+(Days|Subject),
##'                     sleepstudy)$reTrms
##' getReStruc(rt)
##' getReStruc(rt2)
##' @importFrom stats setNames dist .getXlevels
##' @export
getReStruc <- function(reTrms, ss=NULL, aa=NULL, reXterms=NULL, fr=NULL, full_cor=NULL) {

    ## information from ReTrms is contained in cnms, flist elements
    ## cnms: list of column-name vectors per term
    ## flist: data frame of grouping variables (factors)
    ##   'assign' attribute gives match between RE terms and factors
    if (is.null(reTrms)) return(list())

    ## Get info on sizes of RE components

    assign <- attr(reTrms$flist,"assign")
    nreps <- vapply(assign,
                    function(i) length(levels(reTrms$flist[[i]])),
                    0)
    blksize <- diff(reTrms$Gp) / nreps
    ## figure out number of parameters from block size + structure type

    if (is.null(ss)) {
        ss <- rep("us",length(blksize))
    }

    if (is.null(aa)) {
        aa <- rep(NA, length(blksize))
    }

    getRank <- function(cov_name, a) {
        if (cov_name != "rr") return(0)
        if (is.na(a)) return(2) #default rank is 2 [FIXME: don't hard-code here; specify upstream]
        return(a)
    }

    blkrank <- mapply(getRank, ss, aa)

    parFun <- function(struc, blksize, blkrank) {
        switch(as.character(struc),
               "diag" = blksize, # (heterogenous) diag
               "us" = blksize * (blksize+1) / 2,
               "cs" = blksize + 1,
               "ar1" = 2,
               "hetar1" = blksize + 1,
               "ou" = 2,
               "exp" = 2,
               "gau" = 2,
               "mat" = 3,
               "toep" = 2 * blksize - 1,
               "rr" = blksize * blkrank - (blkrank - 1) * blkrank / 2, #rr
               "homdiag" = 1,  ## (homogeneous) diag
               "propto" = blksize * (blksize+1) / 2 + 1, #propto (same as us, plus one extra for proportional param)
               "homcs" = 2,
               "homtoep" = blksize,
               stop(sprintf("undefined number of parameters for covstruct '%s'", struc))
               )
    }
    blockNumTheta <- mapply(parFun, ss, blksize, blkrank, SIMPLIFY=FALSE)

    covCode <- .valid_covstruct[ss];

    ## set simulation code to 'random' for all RE by default
    simCode <- rep(.valid_simcode[["random"]], length(ss))

    ##--replacement for stats::.getXLevels, which fails if `Terms` has no "response" attribute
    .getXlevels<-function (Terms, m) {
        xvars <- vapply(attr(Terms, "variables"), stats:::deparse2, "")[-1L]
        if ((!is.null(attr(Terms, "response")) && (yvar <- attr(Terms, "response")) > 0))
            xvars <- xvars[-yvar]
        if (length(xvars)) {
            xlev <- lapply(m[xvars], function(x) if (is.factor(x))
                levels(x)
            else if (is.character(x))
                levels(as.factor(x)))
            xlev[!vapply(xlev, is.null, NA)]
        }
    }

    ans <- list()
    for (i in seq_along(ss)) {
        tmp <- list(blockReps = nreps[i],
                    blockSize = blksize[i],
                    blockNumTheta = blockNumTheta[[i]],
                    blockCode = covCode[i],
                    simCode = simCode[i],
                    fullCor = as.integer(full_cor[i])
                    )
        if(ss[i] %in% c("ar1", "hetar1")) {
            ## FIXME: Keep this warning ?
            if (any(reTrms$cnms[[i]][1] == "(Intercept)") )
                warning(paste0(ss[i], "() not meaningful with intercept"))
            if (length(.getXlevels(reXterms[[i]],fr))!=1) {
                stop(paste0(ss[i], "() expects a single, factor variable as the time component"))
            }
        } else if(ss[i] == "ou") {
            times <- parseNumLevels(reTrms$cnms[[i]])
            if (ncol(times) != 1)
                stop("'ou' structure is for 1D coordinates only.")
            if (is.unsorted(times, strictly=TRUE))
                stop("'ou' is for strictly sorted times only.")
            tmp$times <- drop(times)
        } else if(ss[i] %in% c("exp", "gau", "mat")){
            coords <- parseNumLevels(reTrms$cnms[[i]])
            tmp$dist <- as.matrix( dist(coords) )
        }
        ans[[i]] <- tmp
    }
    names(ans) <- names(reTrms$Ztlist)

    return(ans)
}

##--From glmmTMB::glmmTMB.R----
###--.noDispersionFamilies: names of families with no dispersion parameter----
.noDispersionFamilies <- c("binomial", "poisson", "truncated_poisson", "bell")

###--.extraParamFamilies: list with number of additional/shape parameters (default = 0)
.extraParamFamilies <- list('1' = c('t', 'tweedie', 'nbinom12', 'skewnormal'),
                            '2' = 'ordbeta')

###--find_psi: function to identify how many extra parameters stat family `f` has
find_psi <- function(f) {
    for (i in seq_along(.extraParamFamilies)) {
        if (f %in% .extraParamFamilies[[i]]) return(as.numeric(i))
    }
    return(0)
}

###--usesDispersion: function to identify whether stat family `f` has a dispersion parameter
## BMB: why not just sigma(f)!=1.0 ... ? (redundant with sigma.glmmTMB)
usesDispersion <- function(f) {
    is.na(match(f, .noDispersionFamilies))
}

###--.classicDispersionFamilies: names of stat families with classic dispersion parameter----
.classicDispersionFamilies <- c("gaussian","Gamma","t")

##--stripReTrms: function to select only desired pieces from results of getXReTrms----
stripReTrms <- function(xrt, whichReTrms = c("cnms","flist","smooth_info"), which="terms") {
    whichReTrms <- intersect(whichReTrms, names(xrt$reTrms)) ## not all models will have smooth_info
    c(xrt$reTrms[whichReTrms],setNames(xrt[which],which))
}

