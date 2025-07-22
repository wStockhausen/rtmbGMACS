#--test estimability of ascnormal1 function parameters by fitting to
#--simulated data using RTMB
require(RTMB);
require(ggplot2);
source(file.path(here::here(),"R/SelectivityFunctions.R"));

verbose=TRUE;
z = seq(25,180,by=5);

pMdZ = 105.0; #--size at which selectivity reaches 1
pWdZ =  20.0; #--width
refZ  = 125;  #--reference *size*
params = list(pWdZ=pWdZ,pMdZ=pMdZ);      #--individual parameters supplied as list objects
params = list(p=c(pWdZ=pWdZ,pMdZ=pMdZ)); #--parameters supplied in a vector
map_   = NULL;

#--reference selectivity----
selr = ascnormal1(z,ascWdZ=pWdZ,ascMdZ=pMdZ,refZ=refZ,debug=verbose);

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

    prd_sel = ascnormal1(z,params=params,refZ=data$refZ,debug=verbose);
    cat("prd_sel = \n"); cat(str(prd_sel),"\n\n");
    print(prd_sel);

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

