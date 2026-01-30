##' From the result of \code{\link{findbars_}} applied to a model formula and
##' and the evaluation frame, create the model matrix, etc. associated with
##' random-effects terms.  See the description of the returned value for a
##' detailed list.
##'
##' @title Create list of structures needed for models with random effects
##' @param sp_form a list of parsed random-effects terms with covariate structure info created using splitForm_RE(formula) [splitForm_RE()]
##' @param model_frame a model frame in which to evaluate these terms
##' @param drop.unused.levels (logical) drop unused factor levels?
##' @param reorder.terms arrange random effects terms in decreasing order of number of groups (factor levels)?
##' @param reorder.vars arrange columns of individual random effects terms in alphabetical order?
##' @param sparse (logical) set up sparse model matrices?
##' @return a list with components
##' \item{Zt}{transpose of the sparse model matrix for the random effects}
##'  \item{Ztlist}{list of components of the transpose of the
##'    random-effects model matrix, separated by random-effects term}
##' \item{Lambdat}{transpose of the sparse relative covariance factor}
##' \item{Lind}{an integer vector of indices determining the mapping of the
##'     elements of the \code{theta} to the \code{"x"} slot of \code{Lambdat}}
##' \item{theta}{initial values of the covariance parameters}
##' \item{lower}{lower bounds on the covariance parameters}
##' \item{flist}{list of grouping factors used in the random-effects terms}
##' \item{cnms}{a list of column names of the random effects according to
##'     the grouping factors}
##' \item{Gp}{a vector indexing the association of
##'    elements of the conditional mode vector
##'    with random-effect terms; if \code{nb} is the vector of numbers
##'    of conditional modes per term (i.e. number of groups times number
##'  of effects per group), \code{Gp} is \code{c(0,cumsum(nb))}
##'     (and conversely \code{nb} is \code{diff(Gp)})}
##' \item{nl}{names of the terms (in the same order as \code{Zt},
##'     i.e. reflecting the \code{reorder.terms} argument)}
##' \item{ord}{an integer vector giving the relationship between the order of the terms in the formula and the terms in the final object (which are ordered by the number of levels in the grouping variable, if \code{reorder.terms} is TRUE)}
##' @details \code{Lambdat}, \code{Lind}, \code{theta}, \code{lower} are likely to
##' be useful only for \code{lme4}; the other terms can be generally useful for
##' constructing mixed-effect models
##' @importFrom Matrix sparseMatrix drop0
## (no methods found in package 'Matrix' for rbind ... ???)
##' @importMethodsFrom Matrix coerce t diag
##' @references \insertRef{lme4}{reformulas})
##' @export
##' @examplesIf FALSE
##' ## (silly/impractical formula, for illustration only)
##' form <- mpg ~ 1 + (1|gear) + (factor(cyl)|gear) + (1 + hp | carb)
##' fr <- model.frame(subbars(form), data = mtcars)
##' rterms <- mkReTrms(findbars(form), fr)
##' names(rterms)
##' ## block sizes (latent variables per block) of each term
##' (nperblock <- lengths(rterms$cnms))
##' ## latent variables per term
##' (nperterm <- diff(rterms$Gp))
##' with(rterms, identical(unname(nl*nperblock), nperterm))
##' ## illustrate reordering of terms
##' dd <- expand.grid(a = 1:7, b = 1:3, c = 1:5, d = 1:9)
##' dd$y <- 1
##' form2 <- y ~ 1 + (1|a) + (1|b) + (1|c) + (1|d)
##' rterms2 <- mkReTrms(findbars(form2), dd, reorder.terms = TRUE)
##' ## reorder elements into original formula order
##' with(rterms2, cnms[order(ord)])
##' ## reorder splitForm output to match mkReTrms components
##' ss <- splitForm(form2)
##' ss$reTrmFormulas[rterms2$ord]
mkReTrms <- function(sp_form,
                     model_frame,
                     drop.unused.levels=TRUE,
                     reorder.terms=TRUE,
                     reorder.vars=FALSE,
                     sparse = TRUE) {
  fr = model_frame; #--TODO: replace `fr` with `model_frame` in code below
  # names(bars) <- reformulas:::barnames(reformulas::no_specials(bars));
  # term.names <- vapply(bars, deparse1, "");
  term.names = names(sp_form$reTrms);

  ##--get component blocks----
  blist <- lapply(sp_form$reTrms, mkBlist, fr, drop.unused.levels,
                  reorder.vars = reorder.vars, sparse = sparse)
  nl <- vapply(blist, `[[`, 0L, "ngroups")   # no. of groups (factor levels) per term
                                    # (in lmer jss:  \ell_i)

  ##--order terms stably by decreasing number of levels in the factor----
  ord <- seq_along(nl)
  if (reorder.terms) {
      if (any(diff(nl) > 0)) {
          ord <- rev(order(nl))
          blist      <- blist     [ord]
          nl         <- nl        [ord]
          term.names <- term.names[ord]
      }
  }
  ##--create Ztlist and Zt----
  Ztlist <- lapply(blist, `[[`, "sm"); #--extract sm's into a list
  Zt <- do.call(rbind, Ztlist)  ## eq. 7, JSS lmer paper
  names(Ztlist) <- term.names
  q <- nrow(Zt)

  ## Create and install Lambdat, Lind, etc.  This must be done after
  ## any potential reordering of the terms.
  colnms <- lapply(blist, `[[`, "colnms")   # list of column names of the
                                            # model matrix per term
  nc <- lengths(colnms)                     # no. of columns per term
                                            # (in lmer jss:  p_i)
  nb <- nc * nl                         # no. of random effects per term
                                        # (in lmer jss:  q_i)
  ## eq. 5, JSS lmer paper
  if (sum(nb) != q) {
      stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)",
                   sum(nb),q))
  }
  boff <- cumsum(c(0L, nb))             # offsets into b

  ll <- list(Zt = drop0(Zt),                    #--TODO: do zero columns need to be dropped
             Gp = unname(c(0L, cumsum(nb))));

  # massage the factor list
  fl <- lapply(blist, `[[`, "ff")
  # check for repeated factors
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  } else asgn <- seq_along(fl)
  names(fl) <- ufn;
  attr(fl, "assign") <- asgn;

  ll$flist  <- fl;
  ll$colnms <- colnms;
  ll$Ztlist <- Ztlist
  ll$nl     <- nl;
  ll$ord    <- ord;
  ll$blist  <- blist;
  ll
} ## {mkReTrms}
