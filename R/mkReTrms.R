##' From the result of \code{\link{findbars}} applied to a model formula and
##' and the evaluation frame, create the model matrix, etc. associated with
##' random-effects terms.  See the description of the returned value for a
##' detailed list.
##'
##' @title Create list of structures needed for models with random effects
##' @param bars a list of parsed random-effects terms with covariate structure info
##' @param model_frame a model frame in which to evaluate these terms
##' @param drop.unused.levels (logical) drop unused factor levels?
##' @param reorder.terms arrange random effects terms in decreasing order of number of groups (factor levels)?
##' @param reorder.vars arrange columns of individual random effects terms in alphabetical order?
##' @param calc.lambdat (logical) compute \code{Lambdat} and \code{Lind} components? (At present these components
##' are needed for \code{lme4} machinery but not for \code{glmmTMB}, and may be large in some cases; see Bates \emph{et al.} 2015
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
##' @importFrom Rdpack reprompt
##' @family utilities
##' @references \insertRef{lme4}{reformulas})
##' @export
##' @examples
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
mkReTrms <- function(bars,
                     model_frame,
                     drop.unused.levels=TRUE,
                     reorder.terms=TRUE,
                     reorder.vars=FALSE,
                     calc.lambdat = TRUE,
                     sparse = NULL) {
  fr = model_frame; #--TODO: replace `fr` with `model_frame` in code below
  # if (!length(bars))
  #   stop("No random effects terms specified in formula",call.=FALSE);
  # stopifnot(is.list(bars), vapply(bars, is.language, NA),
  #           inherits(fr, "data.frame"));
  names(bars) <- reformulas:::barnames(reformulas::no_specials(bars));
  term.names <- vapply(bars, deparse1, "");

  ##--get component blocks----
  blist <- lapply(bars, mkBlist, fr, drop.unused.levels,
                  reorder.vars = reorder.vars, sparse = sparse)
  nl <- vapply(blist, `[[`, 0L, "nlevels")   # no. of levels per term
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
  if (calc.lambdat)
      nth <- as.integer((nc * (nc+1))/2)    # no. of parameters per term       <=number of theta parameters for each RE term, doesn't account for cov. struct.
                                            # (in lmer jss:  ??)
      nb <- nc * nl                         # no. of random effects per term
                                            # (in lmer jss:  q_i)
  ## eq. 5, JSS lmer paper
  if (sum(nb) != q) {
      stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)",
                   sum(nb),q))
  }
  boff <- cumsum(c(0L, nb))             # offsets into b
  if (calc.lambdat) thoff <- cumsum(c(0L, nth))           # offsets into theta
  ## FIXME: should this be done with cBind and avoid the transpose
  ## operator?  In other words should Lambdat be generated directly
  ## instead of generating Lambda first then transposing?
  if (calc.lambdat) {                                                           #--creating Lambdat here----
      mk_b <-function(i) {
          mm <- matrix(seq_len(nb[i]), ncol = nc[i],
                       byrow = TRUE)
          dd <- diag(nc[i])
          ltri <- lower.tri(dd, diag = TRUE)
          ii <- row(dd)[ltri]
          jj <- col(dd)[ltri]
          ## unused: dd[cbind(ii, jj)] <- seq_along(ii)
          data.frame(i = as.vector(mm[, ii]) + boff[i],
                     j = as.vector(mm[, jj]) + boff[i],
                     x = as.double(rep.int(seq_along(ii),
                                           rep.int(nl[i], length(ii))) +
                                   thoff[i]))
      }
      Lambdat <- t(do.call(sparseMatrix,
                           do.call(rbind,
                                   lapply(seq_along(blist), mk_b)))); #--matrix with indices to non-zero elements as values
      Lind <- as.integer(Lambdat@x); #--Lind can "map" theta into Lambdat: Lamdbat@x <- theta[Lind];
      #--NOTE: might need to re-create Lambdat as adsparse? or convert it to adsparse once Lind is obtained
  } else {
      Lambdat <- Lind <- NULL
  }
  thet <- NULL
  if (calc.lambdat)  thet <- numeric(sum(nth))
  ll <- list(Zt = drop0(Zt), theta = thet, Lind = Lind,
             Gp = unname(c(0L, cumsum(nb))))
  ## lower bounds on theta elements are 0 if on diagonal, else -Inf
  ll$lower <- -Inf * (thet + 1)
  if (calc.lambdat) {
      ll$lower[unique(diag(Lambdat))] <- 0
      Lambdat@x[] <- ll$theta[ll$Lind]  # initialize elements of Lambdat
  }
  ll$theta[] <- is.finite(ll$lower) # initial values of theta are 0 off-diagonal, 1 on
  ll$Lambdat <- Lambdat
  # massage the factor list
  fl <- lapply(blist, `[[`, "ff")
  # check for repeated factors
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  } else asgn <- seq_along(fl)
  names(fl) <- ufn
  ## DON'T need fl to be a data.frame ...
  ## fl <- do.call(data.frame, c(fl, check.names = FALSE))
  attr(fl, "assign") <- asgn
  ll$flist <- fl
  ll$colnms <- colnms
  ll$Ztlist <- Ztlist
  ll$nl <- nl
  ll$ord <- ord
  ll
} ## {mkReTrms}
