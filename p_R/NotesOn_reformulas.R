#--decomposing formula with random effects
##--uses "reformulas" R package
require(reformulas);
##--random effects are terms in a formula object represented as
###--a formula followed by a vertical bar or double vertical bars followed by a formula for grouping variables,
###--all enclosed in parentheses.
f1<-y~x+(1+z|g);
f2<-y~x+(1+z||g);

#--`expandDoubleVerts` expands RE terms with double vertical bars into separate, *indepndent* RE terms
reformulas::expandDoubleVerts(f1); #--yields: y ~ x + (1 + z | g) (no change)
reformulas::expandDoubleVerts(f2); #--yields: y ~ x + ((1 | g) + (0 + z | g))

#--`RHSForm` extracts righthand side of a formula
reformulas::RHSForm(f1,as.form=FALSE); #--yields: x+(1+z|g) as language object
reformulas::RHSForm(f1,as.form=TRUE);  #--yields: ~x+(1+z|g) as a formula object w/ "no" response (attr(*, "response")= num 0)

#--`subbars` recursively replaces | and || operators with + on individual terms
##--this provides a formula suitable for using with `model.frame`
subbars(f1); #--yields: y ~ x + (1 + z + g)
subbars(f2); #--yields: y ~ x + (1 + z + g)

#--`nobars` removes random effects terms from a mixed effects formula
##--this provides a formula with fixed effects only
nobars(f1); #--yields: y ~ x
nobars(f2); #--yields: y ~ x

#--`splitForm` splits a mixed model formula into
splitForm(f1); #--yields the a list with the following elements:
# $fixedFormula = y ~ x
# $reTrmFormulas$reTrmFormulas[[1]] = 1 + z | g
# $reTrmAddArgs$reTrmAddArgs[[1]] = us()
# $reTrmClasses[1] = "us"
splitForm(f2); #--yields the a list with the following elements:
# $fixedFormula = y ~ x
# $reTrmFormulas$reTrmFormulas[[1]] = 1 + z | g
# $reTrmAddArgs$reTrmAddArgs[[1]] = diag()
# $reTrmClasses[1] = "diag"

#--`mkReTrms` create the model matrix, etc. associated with random-effects terms.
# mkReTrms(
#   bars,  :a list of parsed random-effects terms
#   fr,    :a model frame in which to evaluate these terms
#   drop.unused.levels = TRUE,
#   reorder.terms = TRUE,     :arrange random effects terms in decreasing order of number of groups (factor levels)
#   reorder.vars = FALSE,     :arrange columns of individual random effects terms in alphabetical order?
#   calc.lambdat = TRUE,      :compute Lambdat and Lind components? (see Bates et al. 2015)
#   sparse = NULL             :create sparse model matrices
# )
# Value: a list with components:
#   Zt      - transpose of the sparse model matrix for the random effects
#   Ztlist  - list of components of the transpose of the random-effects model matrix, separated by random-effects term
#   Lambdat - transpose of the sparse relative covariance factor
#   Lind    - an integer vector of indices determining the mapping of the elements of the theta to the "x" slot of Lambdat
#   theta   - initial values of the covariance parameters
#   lower   - lower bounds on the covariance parameters
#   flist   - list of grouping factors used in the random-effects terms
#   cnms    - a list of column names of the random effects according to the grouping factors
#   Gp      - a vector indexing the association of elements of the conditional mode vector with random-effect terms;
#               if nb is the vector of numbers of conditional modes per term (i.e. number of groups times number of effects per group),
#               Gp is c(0,cumsum(nb)) (and conversely nb is diff(Gp))
#   nl      - names of the terms (in the same order as Zt, i.e. reflecting the reorder.terms argument)
#   ord     - an integer vector giving the relationship between the order of the terms in the formula
#               and the terms in the final object (which are ordered by the number of levels in the grouping variable, if reorder.terms is TRUE)

nobars(f1); #--yields: y ~ x
nobars(f2); #--yields: y ~ x


