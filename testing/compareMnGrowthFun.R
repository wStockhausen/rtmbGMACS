#'
#' @title Compare a gmacs growth function in R and AD contexts
#' @description Function to compare a gmacs growth function in R and AD contexts.
#'
#' @param fr - gmacs growth function
#' @param z - vector of sizes at which to evaluate growth function
#' @param params - vector of "true" parameter values for function
#' @param consts - vector of constants
#' @param scale - gamma distribution scale parameter to add variability
#' @param title - graph title
#'
#' @return nothing
#'
#' @details This function:
#' \enumerate{
#'   \item calculates `fr` using normal R evaluation
#'   \item creates an RTMB tape, evaluates it in an AD context, and prints the jacobian
#'   \item adds variability to the "data" created in 1.
#'   \item creates an RTMB objective function to fit the data
#'   \item fits the data using nlminb
#'   \item prints out the estimated parameters
#'   \item plots a comparison of the R, AD, observed and estimated values
#' }
#'
#' See testing/testGrowthFunctions.R for example usage.
#'
#' @export
#'
compareGrowthFun<-function(fr,z,params,consts,scale,title){
  cat("\n#----------------------------------------------------------\n")
  cat("#  testing",title,"\n");
  gr = fr(z,params,consts);#--evaluate using normal R evaluation
  mi = gr - z;#--mean molt increment
  # ggplot(tibble::tibble(x=z,y=gr),aes(x=x,y=y)) +
  #   geom_point() + geom_line() +
  #   scale_y_continuous(limits=c(min(c(0,gr)),NA));

  #--make RTMB tape
  f<-function(p){fr(z,p,consts)};#--"reduce" function to parameters-only for AD taping
  F = MakeTape(f,numeric(length(params)));
  show(F);                  #--just FYI
  gad = F(params);          #--evaluate using RTMB AD methods
  print(F$jacobian(params));#--print the jacobian (for fun)

  #--test in simple objective function
  set.seed(111);
  data = list(z=z,
              obs=z+rgamma(length(mi),mi,scale=scale), #--post-molt size
              scale=scale);
  objfn<-function(ps){
    getAll(data,ps);
    obs<-OBS(obs);#--observed post-molt size

    prd_g = fr(z,pars,consts);
    REPORT(prd_g);

    obs_mi = obs-z;#--"observed" molt increment
    obs_mi %~% dgamma(prd_g-z,scale=scale); #--adds -log-likelihood to hidden variable `.nll`

    sim_obs = z + obs_mi;
    REPORT(sim_obs);

  }
  #--create model----
  mdl = MakeADFun(objfn,parameters=list(pars=params))
  print(mdl$par);
  print(mdl$fn(params));
  print(mdl$gr(params));
  #--optimize the model----
  #----(won't exactly reproduce `params`, but should be close)
  opt = nlminb(mdl$par,mdl$fn,mdl$gr,hessian=mdl$he);
  print(opt$par);
  print(RTMB::sdreport(mdl));
  dfr = dplyr::bind_cols(tibble::as_tibble(data),
                         `true R`  = gr,
                         `true AD` = gad,
                         est = mdl$report()$prd_g) |>
        tidyr::pivot_longer(c("true R","true AD","est","obs"),
                            names_to="type",values_to="value") |>
        dplyr::mutate(type=factor(type));
  p = ggplot(dfr |> dplyr::filter(type!="obs"),aes(x=z,y=value,colour=type,fill=type)) +
        geom_point(data=dfr |> dplyr::filter(type=="obs")) +
        geom_line() +
        scale_y_continuous(limits=c(min(c(0,dfr$value)),NA)) +
        labs(subtitle=title) +
        wtsPlots::getStdTheme();
  print(p);
  cat("#----------------------------------------------------------\n")
  return(invisible(mdl));
}
