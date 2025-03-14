#'
#' @title Compare a gmacs selectivity function in R and AD contexts
#' @description Function to compare a gmacs selectivity function in R and AD contexts.
#'
#' @param fcn - name of gmacs selectivity function
#' @param z - vector of sizes at which to evaluate selectivity function
#' @param ... - named parameters passed on to selectivity function
#' @param map_ - named list pasted to [TMB::MakeADFun()] identifying fixed parameters
#' @param title - plot title
#' @param verbose - flag (T/F) to print/plot extra information
#'
#' @return nothing
#'
#' @details This function:
#' \enumerate{
#'   \item calculates `fcn` using normal R evaluation
#'   \item adds noise to the "data" created in 1.
#'   \item creates an RTMB objective function to fit the data
#'   \item fits the data using nlminb
#'   \item prints out the estimated parameters (if verbose=TRUE)
#'   \item returns a ggplot2 with a comparison of the R, AD, observed and estimated values
#' }
#'
#' See testing/testIndivSelFunctions1.R for example usage.
#'
#' @export
#'
compareSelFun1<-function(fcn,z,...,map_=list(),title=fcn,verbose=FALSE){
  if (verbose) {
    cat("\n#----------------------------------------------------------\n")
    cat("#  testing",title,"\n");
  }
  params = list(...);
  if (verbose) print(params);

  dZ = z[2] - z[1];

  #--calculate the function using normal R evaluation (can pass ... here)----
  if (fcn=="const_sel")##--const_sel----
    selr = const_sel(z,...,verbose=verbose);#--constant

  if (fcn=="asclogistic")##--asclogistic----
    selr = asclogistic(z,...,verbose=verbose);#--ascending logistic
  if (fcn=="asclogistic1")
    selr = asclogistic1(z,...,verbose=verbose);#--alternative ascending logistic
  if (fcn=="asclogistic5095")
    selr = asclogistic5095(z,...,verbose=verbose);#--alternative ascending logistic
  if (fcn=="asclogistic50D95")
    selr = asclogistic50D95(z,...,verbose=verbose);#--alternative ascending logistic

  if (fcn=="ascnormal1")##--ascnormal1----
    selr = ascnormal1(z,...,dZ=dZ,verbose=verbose);#--ascending normal
  if (fcn=="ascnormal2")##--ascnormal2----
    selr = ascnormal2(z,...,dZ=dZ,verbose=verbose);#--ascending normal
  if (fcn=="ascnormal2a")##--ascnormal2a----
    selr = ascnormal2a(z,...,dZ=dZ,verbose=verbose);#--ascending normal
  if (fcn=="ascnormal2b")##--ascnormal2b----
    selr = ascnormal2b(z,...,dZ=dZ,verbose=verbose);#--ascending normal
  if (fcn=="ascnormal3")##--ascnormal3----
    selr = ascnormal3(z,...,dZ=dZ,verbose=verbose);#--ascending normal

  if (fcn=="dbllogistic")##--dbllogistic----
    selr = dbllogistic(z,...,verbose=verbose);#--double logistic
  if (fcn=="dbllogistic5095")##--dbllogistic5095----
    selr = dbllogistic5095(z,...,verbose=verbose);#--alternative double logistic

    if (fcn=="dblnormal4")##--dblnormal4----
      selr = dblnormal4(z,...,dZ=dZ,verbose=verbose);#--4-parameter double normal
    if (fcn=="dblnormal4a")##--dblnormal4a----
      selr = dblnormal4a(z,...,dZ=dZ,verbose=verbose);#--4-parameter double normal
    if (fcn=="dblnormal6")##--dblnormal6----
      selr = dblnormal6(z,...,dZ=dZ,verbose=verbose);#--6-parameter double normal

  if (fcn=="stackedLogistic1")##--stackedLogistic1----
    selr = stackedLogistic1(z,...,verbose=verbose);

  if (fcn=="selSpline")##--selSpline----
    selr = selSpline(z,...,verbose=verbose);
  if (fcn=="selSplineClmpd")##--selSplineClmpd----
    selr = selSplineClmpd(z,...,verbose=verbose);
  if (fcn=="selSplineClmpdRight")##--selSplineClmpdRight----
    selr = selSplineClmpdRight(z,...,verbose=verbose);
  if (fcn=="selSplineClmpdLeft")##--selSplineClmpdLeft----
    selr = selSplineClmpdLeft(z,...,verbose=verbose);

  print(z);
  print(selr);
  p = ggplot(tibble::tibble(z=z,s=selr),aes(x=z,y=s)) +
        geom_point() + geom_line() +
        scale_y_continuous(limits=c(min(c(0,selr)),NA)) +
        labs(subtitle="normal R");
  if (verbose) print(p);

  # #--make RTMB tape
  # if (verbose) cat("making tape.\n")
  # f<-function(p){fr(z,p,verbose=verbose)};#--"reduce" function to parameters-only for AD taping
  # F = MakeTape(f,numeric(length(unlist(params))));
  # if (verbose) show(F);     #--just FYI
  # selad = F(params);        #--evaluate using RTMB AD methods
  # if (verbose) print(F$jacobian(params));#--print the jacobian (for fun)

  #--create data to test fitting----
  set.seed(111);
  data = list(fcn=fcn,z=z,obs=selr+rnorm(length(selr),0,0.05));

  #--define simple objective function----
  if (verbose) cat("setting up objective function.\n")
  objfn<-function(pars){
    z<-data$z;
    obs<-OBS(data$obs);

    dZ = z[2]-z[1];

    if (data$fcn=="const_sel")##--const_sel----
      prd_sel = const_sel(z,pars$pCnst,verbose=verbose);

    if (data$fcn=="asclogistic")##--asclogistic----
      prd_sel = asclogistic(z,pars$pZ50,pars$pSlp,pars$pRefZ,verbose=verbose);
    if (data$fcn=="asclogistic1")
      prd_sel = asclogistic1(z,pars$pZ50,pars$pWdZ,pars$pRefZ,verbose=verbose);
    if (fcn=="asclogistic5095")
      prd_sel = asclogistic5095(z,pars$pZ50,pars$pZ95,pars$pRefZ,verbose=verbose);#--alternative ascending logistic
    if (fcn=="asclogistic50D95")
      prd_sel = asclogistic50D95(z,pars$pZ50,pars$pZ9550,pars$pRefZ,verbose=verbose);#--alternative ascending logistic

    if (fcn=="ascnormal1")##--ascnormal1----
      prd_sel = ascnormal1(z,pars$pAscZ1,pars$pAscWdZ,dZ,verbose=verbose);#--ascending normal
    if (fcn=="ascnormal2")##--ascnormal2----
      prd_sel = ascnormal2(z,pars$pAscZ1,pars$pAscRefS,pars$pRefZ,dZ,verbose=verbose);#--ascending normal
    if (fcn=="ascnormal2a")##--ascnormal2a----
      prd_sel = ascnormal2a(z,pars$pAscZ1,pars$pZatRefS,pars$pRefS,dZ,verbose=verbose);#--ascending normal
    if (fcn=="ascnormal2b")##--ascnormal2b----
      prd_sel = ascnormal2b(z,pars$pAscZ1,pars$pDZ2RefS,pars$pRefS,dZ,verbose=verbose);#--ascending normal
    if (fcn=="ascnormal3")##--ascnormal3----
      prd_sel = ascnormal3(z,pars$pDZ1,pars$pSatZ2,pars$pMxZ1,pars$pRefZ2,dZ,verbose=verbose);#--ascending normal

    if (fcn=="dbllogistic")##--dbllogistic----
      prd_sel = dbllogistic(z,pars$pAscZ50,pars$pAscSlp,pars$pDscZ50,pars$pDscSlp,pars$pRefZ,verbose=verbose);#--double logistic
    if (fcn=="dbllogistic5095")##--dbllogistic5095----
      prd_sel = dbllogistic5095(z,pars$pAscZ50,pars$pAscZ95,pars$pDscZ95,pars$pDscZ50,pars$pRefZ,verbose=verbose);#--alternative double logistic

    if (fcn=="dblnormal4")##--dblnormal4----
      prd_sel = dblnormal4(z,pars$pAscZ1,pars$pAscWd,pars$pDscDZ,pars$pDscWd,dZ,verbose=verbose);#--4-parameter double normal
    if (fcn=="dblnormal4a")##--dblnormal4a----
      prd_sel = dblnormal4a(z,pars$pAscZ1,pars$pAscWd,pars$pDscSclDZ,pars$pDscWd,pars$pRefZ,dZ,verbose=verbose);#--4-parameter double normal
    if (fcn=="dblnormal6")##--dblnormal6----
      prd_sel = dblnormal6(z,pars$pAscZ1,pars$pAscWd,pars$pDscZ1,pars$pDscWd,pars$pAscFlr,pars$pDscFlr,dZ,verbose=verbose);#--6-parameter double normal

    if (fcn=="stackedLogistic1")##--stackedLogistic1----
      prd_sel = stackedLogistic1(z,pars$pMnZ1,pars$pSdZ1,pars$pMnZ2,pars$pSdZ2,pars$pOmga,verbose=verbose);

    if (fcn=="selSpline")##--selSpline----
      prd_sel = selSpline(z,pars$params,pars$knots,verbose=verbose);
    if (fcn=="selSplineClmpd")##--selSplineClmpd----
      prd_sel = selSplineClmpd(z,pars$params,pars$knots,verbose=verbose);
    if (fcn=="selSplineClmpdRight")##--selSplineClmpdRight----
      prd_sel = selSplineClmpdRight(z,pars$params,pars$knots,verbose=verbose);
    if (fcn=="selSplineClmpdLeft")##--selSplineClmpdLeft----
      prd_sel = selSplineClmpdLeft(z,pars$params,pars$knots,verbose=verbose);

    ad_pred_sel = AD(1)*prd_sel;#--just to get a different object for ADREPORTING
    REPORT(prd_sel);

    obs %~% dnorm(prd_sel,1); #--adds -log-likelihood to hidden variable `.nll`

    ADREPORT(ad_pred_sel);

    return(.nll);
  }

  #--create RTMB model----
  if (verbose) cat("MakeADFun'ing\n");
  if (verbose) print(params);
  mdl = MakeADFun(objfn,parameters=params,map=map_,silent=verbose);
  if (verbose) {
    print(mdl$par);
    print(mdl$fn());
    print(mdl$gr());
  }
  #--optimize the model----
  #----(won't exactly reproduce `params`, but should be close)
  if (verbose) cat("Optimizing model\n");
  opt = nlminb(mdl$par,mdl$fn,mdl$gr,hessian=mdl$he);
  if (verbose) {
    print(opt$par);
    print(RTMB::sdreport(mdl));
  }
  #--plot the results----
  if (verbose) cat("Plotting results\n");
  dfr = dplyr::bind_cols(tibble::as_tibble(data),
                         `true R`  = selr,
                         est = mdl$report()$prd_sel) |>
        tidyr::pivot_longer(c("true R","est","obs"),
                            names_to="type",values_to="value") |>
        dplyr::mutate(type=factor(type));
  p = ggplot(dfr |> dplyr::filter(type=="true R"),aes(x=z,y=value)) +
        geom_point(data=dfr |> dplyr::filter(type=="obs"),colour="red") +
        geom_point(colour="blue") +
        geom_line(colour="blue") +
        geom_line(data=dfr |> dplyr::filter(type=="est"),colour="green",linetype=3) +
        scale_y_continuous(limits=c(min(c(0,dfr$value)),NA)) +
        labs(subtitle=title,x="size",y="selectivity") +
        wtsPlots::getStdTheme();
  print(p);
  if (verbose) cat("#----------------------------------------------------------\n")
  return(list(obj=mdl,p=p));
}

