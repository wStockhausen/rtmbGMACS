##--Based on glmmTMB:getXReTerms create RE terms from formula----
##' @param formula formula, possibly containing random effects
##' @param model_frame full model frame
##' @param contrasts a list of contrasts (see ?glmmTMB)
##' @param sparse (logical) return sparse model matrix?
##' @return a list composed of
##' \item{Z}{design matrix for random effects}
##' \item{reTrms}{output from \code{mkReTrms}}
##' \item{ss}{splitform of the formula}
##' \item{aa}{additional arguments, used to obtain rank}
##' \item{terms}{terms for the fixed effects}
##' \item{reXterms}{terms for the model matrix in each RE term}
##' \item{formula}{input formula}
##' \item{model_frame}{input model frame}
##' \item{model_matrix}{combination of input model frame and Z matrix}
##'
##' @importFrom stats model.matrix contrasts
##' @importFrom reformulas inForm findbars nobars noSpecials sub_specials addForm findbars_x anySpecial RHSForm RHSForm<- extractForm reOnly no_specials splitForm addForm0 makeOp
##'
##' @export
##'
##' hacked from getXReTrms
createParamsInfo_REs <- function(formula,
                                 model_frame,
                                 contrasts=NULL,
                                 sparse=TRUE) {

  reform <- getREs(formula);
  nobs <- nrow(model_frame);
  if (is.null(reform)){
    #--no random effects, return NULL (TODO: better idea??)
    #--but for now:
    ##--glmmTMB created these if no fixedform:
    Z <-Matrix::Matrix(numeric(0),nrow=nobs,ncol=0,sparse=sparse);
    offset <- rep(0,nobs);          #--offsets are 0
    p <- RTMB::AD(numeric(0));
    re_form     = reform;
    model_matrix = dplyr::bind_cols(model_frame,as.matrix(Z));
    #if (!has_re && !has_smooths) { ##--TODO: anything we need to keep here?
        reTrms <- reXterms <- NULL
        Z <- new("dgCMatrix",Dim=c(as.integer(nobs),0L)) ## matrix(0, ncol=0, nrow=nobs)
        aa <- integer(0) #added for rr to get rank
        ss <- integer(0)
    #}
    ## dummy elements TODO: NULL above--is this better?  Think this can be better organized
    reTrms <- list(Ztlist = list(), flist = list(), cnms = list(), theta = list(),
                   Zt=list(), Lind=list(),Gp=list(),lower=list(),Lambdat=list(),nl=list(),ord=list()); #--WTS: added these
    return(namedList(p, Z, re_form, contrasts, offset, formula, model_frame, model_matrix));
  }

  terms <- NULL ## make sure it's empty in case we don't set it

  ###--make reTrms----
  ## no_specials so that mkReTrms can handle it TODO: think maybe want to include specials to allow cov structures to be evaluated in mkReTrms
  #  bars = reformulas::no_specials(reformulas::findbars_x(formula));
  bars = reformulas::findbars_x(reform,specials=findValidCovStructs(),default.special="diag"); # <-allows covariance structures info to be passed
  reTrms <- mkReTrms(bars,
                     model_frame,
                     reorder.terms=FALSE,
                     sparse = sparse);

  ##--split `formula` into `fixedFormula` and reTrmFormulas, reTrmAddArgs lists----
  ###--elements in reTrmFormulas are formulas for RE or smooth terms as "effect|group" or "variable"
  #sp_form <- reformulas::splitForm(formula, specials = findValidCovStructs())             #--covariance types extracted here!!----
  sp_form <- splitForm_RE(formula); #--TODO: change defaultTerm to "diag"?

  ##--create additional arguments (`aa`) list
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

  ## if(is.null(rankX.chk <- control[["check.rankX"]]))
  ## rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
  ## X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
  ## if(is.null(scaleX.chk <- control[["check.scaleX"]]))
  ##     scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
  ## X <- checkScaleX(X, kind=scaleX.chk)

  ## list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
  ##      wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))

  return(namedList(Z, reTrms, ss, aa, terms, offset, reXterms, formula));
} #--getXReTerms

