#'
#' @title Compare a gmacs selectivity function in R and AD contexts
#' @description Function to compare a gmacs selectivity function in R and AD contexts.
#'
#' @param fr - gmacs selectivity function
#' @param z - vector of sizes at which to evaluate selectivity function
#' @param params - vector of "true" parameter values for function
#' @param consts - vector of constants
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
compareSelFun<-function(fr,z,params,consts,title){
  cat("\n#----------------------------------------------------------\n")
  cat("#  testing",title,"\n");
  selr = fr(z,params,consts);#--evaluate using normal R evaluation
  pg = ggplot(tibble::tibble(z=z,s=selr),aes(x=z,y=s)) +
        geom_point() + geom_line() +
        scale_y_continuous(limits=c(min(c(0,selr)),NA)) +
        labs(subtitle="normal R");
  print(pg);

  #--make RTMB tape
  f<-function(p){fr(z,p,consts)};#--"reduce" function to parameters-only for AD taping
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

    prd_sel = fr(z,pars,consts);
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
  p = ggplot(dfr |> dplyr::filter(type=="true AD"),aes(x=z,y=value)) +
        geom_point(data=dfr |> dplyr::filter(type=="obs"),colour="red") +
        geom_point(data=dfr |> dplyr::filter(type=="true R"),colour="blue") +
        geom_line(colour="green") +
        geom_line(data=dfr |> dplyr::filter(type=="est"),colour="green",linetype=3) +
        scale_y_continuous(limits=c(min(c(0,dfr$value)),NA)) +
        labs(subtitle=title) +
        wtsPlots::getStdTheme();
  print(p);
  cat("#----------------------------------------------------------\n")
  return(invisible(NULL));
}

#' #' @param params - vector of {knots,values} for spline function
#' #' @param sd - sd of observation noise
#' compareSplineSelFun<-function(fr,z,params,ref,title,dfr){
#'   cat("\n#----------------------------------------------------------\n")
#'   cat("#  testing",title,"\n");
#'   x = rep(seq(25,175,25),1);
#'   y = sin(2*pi*x/100);
#'   params = c(x,y);
#'   #--create "true" spline
#'   spl<-function(z,params,ref=NULL){
#'     np = length(params)/2;
#'     xk = params[1:np];   #--knots
#'     yk = params[1:np+np];#--values at knots
#'     st = splinefun(xk,yk,method="natural");#--true spline
#'     xp = seq(from=dplyr::first(xk),to=dplyr::last(xk),length=100);
#'     dfrxp = dplyr::bind_rows(
#'               tibble::tibble(x=xp,y=st(xp),name="fine"),
#'               tibble::tibble(x=z, y=st(z), name="at z's")
#'             );
#'     print(ggplot(dfrxp,aes(x=x,y=y,colour=name)) + geom_line());
#'     return(st);
#'   }
#'   splf = spl(z,c(x,y),FALSE);
#'   ys = splf(z,0)
#'
#'   selr = fr(z,params,ref);#--evaluate using normal R evaluation
#'   # ggplot(tibble::tibble(z=z,s=selr),aes(x=z,y=s)) +
#'   #   geom_point() + geom_line() +
#'   #   scale_y_continuous(limits=c(min(c(0,selr)),NA));
#'
#'   #--make RTMB tape
#'   f<-function(p){splf(z,p,ref)};#--"reduce" function to parameters-only for AD taping
#'   F = MakeTape(f,numeric(length(params)));
#'   show(F);                  #--just FYI
#'   selad = F(params);        #--evaluate using RTMB AD methods
#'   print(F$jacobian(params));#--print the jacobian (for fun)
#'
#'   #--test in simple objective function
#'   set.seed(111);
#'   data = list(z=z,obs=selr+rnorm(length(selr),0,0.05));
#'   objfn<-function(ps){
#'     getAll(data,ps);
#'     obs<-OBS(obs);
#'
#'     prd_sel = fr(z,pars,ref);
#'     REPORT(prd_sel);
#'
#'     obs %~% dnorm(prd_sel,1); #--adds -log-likelihood to hidden variable `.nll`
#'   }
#'   #--create model----
#'   mdl = MakeADFun(objfn,parameters=list(pars=params))
#'   print(mdl$par);
#'   print(mdl$fn(params));
#'   print(mdl$gr(params));
#'   #--optimize the model----
#'   #----(won't exactly reproduce `params`, but should be close)
#'   opt = nlminb(mdl$par,mdl$fn,mdl$gr,hessian=mdl$he);
#'   print(opt$par);
#'   print(RTMB::sdreport(mdl));
#'   dfr = dplyr::bind_cols(tibble::as_tibble(data),
#'                          `true R`  = selr,
#'                          `true AD` = selad,
#'                          est = mdl$report()$prd_sel) |>
#'         tidyr::pivot_longer(c("true R","true AD","est","obs"),
#'                             names_to="type",values_to="value") |>
#'         dplyr::mutate(type=factor(type));
#'   p = ggplot(dfr |> dplyr::filter(type=="true AD"),aes(x=z,y=value)) +
#'         geom_point(data=dfr |> dplyr::filter(type=="obs"),colour="red") +
#'         geom_point(data=dfr |> dplyr::filter(type=="true R"),colour="blue") +
#'         geom_line(colour="green") +
#'         geom_line(data=dfr |> dplyr::filter(type=="est"),colour="green",linetype=3) +
#'         scale_y_continuous(limits=c(min(c(0,dfr$value)),NA)) +
#'         labs(subtitle=title) +
#'         wtsPlots::getStdTheme();
#'   print(p);
#'   cat("#----------------------------------------------------------\n")
#'   return(invisible(NULL));
#' }
