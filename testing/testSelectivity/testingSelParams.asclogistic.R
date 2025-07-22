#--test estimability of asclogistic function parameters by fitting to
#--simulated data using RTMB
require(RTMB);
require(ggplot2);
source(file.path(here::here(),"R/SelectivityFunctions.R"));

verbose=TRUE;
z = seq(25,180,by=5);

z50 = 100.00; #--size at which selectivity = 0.50 (logit-scale mean)
slp =   0.05; #--slope at z50
refZ  = 125;  #--reference *size*
params = list(z50=z50,slp=slp);      #--individual parameters supplied as list objects
params = list(p=c(z50=z50,slp=slp)); #--parameters supplied in a vector
map_   = NULL;

#--reference selectivity----
selr = asclogistic(z,pZ50=z50,pSlp=slp,refZ=refZ,debug=verbose);

  #--create data to fit----
  set.seed(111);
  data = list(z=z,obs=selr+rnorm(length(selr),0,0.05),
              refZ=refZ);

  #--create model/objective function----
  objfn<-function(pars){
    z<-data$z;
    dZ = z[2]-z[1];
    obs<-OBS(data$obs);

#    params = c(pars$z50,pars$slp);
    params = pars$p;
    cat("params = \n"); cat(str(params),"\n\n");
    print(params);

    prd_sel = asclogistic(z,params=params,refZ=data$refZ,debug=verbose);
#    prd_sel = asclogistic(z,pZ50=pars$z50,pSlp=pars$slp,refZ=data$refZ,debug=verbose);
    cat("prd_sel = \n"); cat(str(prd_sel),"\n\n");
    print(prd);

    ADREPORT(pars$z50);
    ADREPORT(pars$slp);
    ADREPORT(pars$refZ);

    REPORT(prd_sel);

    obs %~% dnorm(prd_sel,1); #--adds -log-likelihood to hidden variable `.nll`
  }
  objfn(params);
  selad = selr;

  #--MakeADFun the model----
  if (verbose) cat("MakeADFun'ing\n");
  mdl = MakeADFun(objfn,parameters=params,map=map_);
  if (verbose) {
    print(mdl$par);
    print(mdl$fn());
    print(mdl$gr());
  }
  #--optimize the model----
  #----(won't exactly reproduce `params`, but should be close)
  opt = nlminb(mdl$par,mdl$fn,mdl$gr,hessian=mdl$he);
  if (verbose) {
    print(opt$par);
    print(RTMB::sdreport(mdl));
    source(file.path(here::here(),"R/getDFR_sdreport.R"))
    print(getDFR_sdreport(mdl));
  }
  #----plot the model results
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
        geom_vline(xintercept=refZ,linetype=3) +
        scale_y_continuous(limits=c(min(c(0,dfr$value)),NA)) +
        labs(subtitle="",x="size",y="selectivity") +
        wtsPlots::getStdTheme();
  print(p);

