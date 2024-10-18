#'
#' @title Compare a gmacs selectivity function in R and AD contexts
#' @description Function to compare a gmacs selectivity function in R and AD contexts.
#'
#' @param fr - gmacs selectivity function
#' @param z - vector of sizes at which to evaluate selectivity function
#' @param params - vector of "true" parameter values for function
#' @param ref - reference value
#'
#' @return nothing
#'
#' @details This function:
#' \enumerate{
#'   \item calculates `fr` using normal R evaluation
#'   \item creates an RTMB tape, evaluates it in an AD context, and prints the jacobian
#'   \item adds noise to the "data" created in 1.
#'   \item creates an RTMB objective function to fit the data
#'   \item fits the data using nlminb
#'   \item prints out the estimated parameters
#'   \item plots a comparison of the R, AD, observed and estimated values
#' }
#'
#' See testing/testSelFunctions.R for example usage.
#'
#' @export
#'
compareSelFun<-function(fr,z,params,ref,title){
  cat("\n#----------------------------------------------------------\n")
  cat("#  testing",title,"\n");
  selr = fr(z,params,ref);#--evaluate using normal R evaluation
  # ggplot(tibble::tibble(z=z,s=selr),aes(x=z,y=s)) +
  #   geom_point() + geom_line() +
  #   scale_y_continuous(limits=c(min(c(0,selr)),NA));

  #--make RTMB tape
  f<-function(p){fr(z,p,ref)};#--"reduce" function to parameters-only for AD taping
  F = MakeTape(f,numeric(length(params)));
  show(F);                  #--just FYI
  selad = F(params);        #--evaluate using RTMB AD methods
  print(F$jacobian(params));#--print the jacobian (for fun)

  #--test in simple objective function
  set.seed(111);
  data = list(z=z,obs=selr+rnorm(length(selr),0,0.05));
  objfn<-function(ps){
    getAll(data,ps);
    obs<-OBS(obs);

    prd_sel = fr(z,pars,ref);
    REPORT(prd_sel);

    obs %~% dnorm(prd_sel,1); #--adds -log-likelihood to hidden variable `.nll`
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
                         `true R`  = selr,
                         `true AD` = selad,
                         est = mdl$report()$prd_sel) |>
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
  return(invisible(NULL));
}
