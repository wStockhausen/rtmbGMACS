#--test parameterization using glmmTMB approach to specifications
require(glmmTMB);

dirPrj = file.path(getwd(),"../.."); #--rstudioapi::getActiveProject();
source(file.path(dirPrj,"R","DimensionsUtilities.R"),local=TRUE);
source(file.path(dirPrj,"R","DimensionsFunctions.R"),local=TRUE);
source(file.path(dirPrj,"R","ragged_array.R"),local=TRUE);

#--define some dimensions
dZs = list( MALE=  list(IMMATURE=list(`NEW SHELL`=seq(25,180,5),
                                      `OLD SHELL`=seq(25,180,5)),
                          MATURE=list(`NEW SHELL`=seq(25,180,5),
                                      `OLD SHELL`=seq(25,180,5))),
          FEMALE=list(IMMATURE=list(`NEW SHELL`=seq(25,130,5),
                                    `OLD SHELL`=seq(25,130,5)),
                        MATURE=list(`NEW SHELL`=seq(25,130,5),
                                    `OLD SHELL`=seq(25,130,5))));
attr(dZs,"dmnms")<-c("x","m","s","z");
dfrZs = createSparseDimsMap(z=dZs);
dfrYZs = createSparseDimsMap(y=c(`2020`=2020,`2021`=2021,`2022`=2022,`2023`=2023),
                             z=dZs);

#--glmmTMB specification----
##--ziformula:----
###   a one-sided (i.e., no response variable) formula for zero-inflation combining fixed and random effects: 
###   the default ~0 specifies no zero-inflation. Specifying ~. sets the zero-inflation formula identical to 
###   the right-hand side of formula (i.e., the conditional effects formula); terms can also be added or subtracted. 
###   When using ~. as the zero-inflation formula in models where the conditional effects formula contains an offset term, 
###   the offset term will automatically be dropped. The zero-inflation model uses a logit link.
##--dispformula:----	
###   a one-sided formula for dispersion combining fixed and random effects: 
###   the default ~1 specifies the standard dispersion given any family. The argument is ignored for families 
###   that do not have a dispersion parameter. For an explanation of the dispersion parameter for each family, see sigma. 
###   The dispersion model uses a log link. In Gaussian mixed models, dispformula=~0 fixes the residual variance to be 0 
###   (actually a small non-zero value), forcing variance into the random effects. 
###   The precise value can be controlled via control=glmmTMBControl(zero_dispval=...); 
###   the default value is sqrt(.Machine$double.eps).
##--weights: ----
###   weights, as in glm. Not automatically scaled to have sum 1.
##--contrasts: ----
###   an optional list, e.g., list(fac1="contr.sum"). See the contrasts.arg of model.matrix.default.
##--start:----
###   starting values, expressed as a list with possible components 
###   beta, betazi, betadisp (fixed-effect parameters for conditional, zero-inflation, dispersion models); 
###   b, bzi, bdisp (conditional modes for conditional, zero-inflation, and dispersion models); 
###   theta, thetazi, thetadisp (random-effect parameters, on the standard deviation/Cholesky scale, for conditional, z-i, and disp models); 
###   psi (extra family parameters, e.g., shape for Tweedie models).
##--map: ----
###   a list specifying which parameter values should be fixed to a constant value rather than estimated. 
###   `map` should be a named list containing factors corresponding to a subset of the internal parameter names (see start parameter). 
###    Distinct factor values are fitted as separate parameter values, NA values are held fixed: e.g., 
###    map=list(beta=factor(c(1,2,3,NA))) would fit the first three fixed-effect parameters of the conditional model and 
###    fix the fourth parameter to its starting value. 
###    In general, users will probably want to use start to specify non-default starting values for fixed parameters. 
###    See MakeADFun for more details.
##--priors:----
###   a data frame of priors, in a similar format to that accepted by the brms package; see priors
##
##--Note:---- 
###   By default, vector-valued random effects are fitted with unstructured (general symmetric positive definite) variance-covariance matrices. 
###   Structured variance-covariance matrices can be specified in the form struc(terms|group), where struc is one of
###   diag (diagonal, heterogeneous variance)
###   ar1 (autoregressive order-1, homogeneous variance)
###   cs (compound symmetric, heterogeneous variance)
###   ou (* Ornstein-Uhlenbeck, homogeneous variance)
###   exp (* exponential autocorrelation)
###   gau (* Gaussian autocorrelation)
###   mat (* Matérn process correlation)
###   toep (* Toeplitz)
###   rr (reduced-rank/factor-analytic model)
###   homdiag (diagonal, homogeneous variance)
###   propto (* proportional to user-specified variance-covariance matrix)
##
##--Note:----
###   Smooths taken from the mgcv package can be included in glmmTMB formulas using s; 
###   these terms will appear as additional components in both the fixed and the random-effects terms. 
###   This functionality is experimental for now. We recommend using REML=TRUE. 
###   See s for details of specifying smooths (and smooth2random and the appendix of Wood (2004) for technical details).

#--create test formula
f = sparse_idx~1 + x + s(z,k=6) + (1+m|y);
glF = glmmTMB::glmmTMB(f,
                       dfrYZs |>  dplyr::mutate(z=as.numeric(z)),
                       ziformula=~0,      #--formula (1-sided) for zero inflation 
                       dispformula=~1,    #--formula (1-sided) for dispersion, combining fixed and REs
                       weights=NULL,
                       contrasts=NULL,
                       priors=NULL,
                       map=NULL,
                       doFit=FALSE
                       )

