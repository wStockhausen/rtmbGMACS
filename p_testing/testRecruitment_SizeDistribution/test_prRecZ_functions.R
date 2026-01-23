prRecZ2<-function(mnZ,wdZ,zMn,zMx,zBs,dZ){  #--now prRecZ1 in R/calcRecruitment_SizeDistribution.R
  prZ = AD(array(0,dim=length(zBs)));
  if (RTMB:::ad_context()){
    nzMn = RTMB:::getValues(zMn);#--zMn is fixed, so this should not be a problem
    nzMx = RTMB:::getValues(zMx);#--zMx is fixed, so this should not be a problem
  } else {
    nzMn = zMn;
    nzMx = zMx;
  }
  idx = which((nzMn<=zBs)&(zBs<=nzMx));
  shp = (mnZ^2)/(wdZ^2);
  scl = (wdZ^2)/mnZ;
  prZmn = pgamma(zMn-dZ/2,shape=shp,scale=scl);
  prZmx = pgamma(zMx+dZ/2,shape=shp,scale=scl);
  prZ[idx] = pgamma(zBs[idx]+dZ/2,shape=shp,scale=scl) -
                 pgamma(zBs[idx]-dZ/2,shape=shp,scale=scl);
  # cat("prZ before scaling\n")
  # print(prZ);
  # print(prZmx- prZmn);
  #--normalize and apply squarewave window so sum  across zBs is 1 for each population category
  prZ[idx] = prZ[idx]/(prZmx- prZmn);
  # cat("prZ after scaling\n")
  # print(prZ);
  # cat("sum(prZ):",sum(prZ),"\n");
  return(prZ);
}

verbose=TRUE;
z = seq(25,80,by=5);
dZ = z[2]-z[1];

mnZ = 45.0;  #--mean of recruitment size distribution
wdZ = 15.0; #--width of recruitment size distribution
zMn = 30;  #--min size
zMx = 65;  #--min size
params = list(mnZ=mnZ,wdZ=wdZ,zMn=zMn,zMx=zMx);
map_   = list(zMn=factor(NA),zMx=factor(NA));

prR = prRecZ2(mnZ,wdZ,zMn,zMx,z,dZ);
sum(prR);

#--create data to fit----
set.seed(111);
data = list(z=z,obs=prR+rnorm(length(prR),0,0.05));

#--create model/objective function----
objfn<-function(pars){
  z<-data$z;
  dZ = z[2]-z[1];
  obs<-OBS(data$obs);

  #--either----
  # getAll(pars);
  # prd_sel = asclogistic(z,z50,slp,refZ,verbose);
  # ADREPORT(z50);
  # ADREPORT(slp);
  # ADREPORT(refZ);
  #--or----
  prdR = prRecZ2(pars$mnZ,pars$wdZ,pars$zMn,pars$zMx,z,dZ);
  ADREPORT(pars$mnZ);
  ADREPORT(pars$wdZ);
  ADREPORT(pars$zMn);
  ADREPORT(pars$zMx);

  REPORT(prdR);

  obs %~% dnorm(prdR,1); #--adds -log-likelihood to hidden variable `.nll`
}
objfn(params);

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
  source(file.path(rstudioapi::getActiveProject(),"R/getDFR_sdreport.R"))
  print(getDFR_sdreport(mdl));
}
  dfr = dplyr::bind_cols(tibble::as_tibble(data),
                         `true R`  = prR,
                         est = mdl$report()$prdR) |>
        tidyr::pivot_longer(c("true R","est","obs"),
                            names_to="type",values_to="value") |>
        dplyr::mutate(type=factor(type));
  p = ggplot(dfr |> dplyr::filter(type=="true R"),aes(x=z,y=value)) +
        geom_point(data=dfr |> dplyr::filter(type=="obs"),colour="red") +
        geom_point(colour="blue") +
        geom_line(colour="blue") +
        geom_line(data=dfr |> dplyr::filter(type=="est"),colour="green",linetype=3) +
        geom_vline(xintercept=zMn,linestyle=3) +
        geom_vline(xintercept=zMx,linestyle=3) +
        scale_y_continuous(limits=c(min(c(0,dfr$value)),NA)) +
        labs(subtitle="",x="size",y="pr(recruitment)") +
        wtsPlots::getStdTheme();
  print(p);

