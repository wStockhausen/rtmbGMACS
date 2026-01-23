#--estimate allometry
##--NOTE: TYPE = "pwrLaw2" converges successfully; "pwrLaw1" does not!
require(ggplot2);
dirThs = dirname(rstudioapi::getSourceEditorContext()$path);

#--get allometric data----
dfrZW  = wtsUtilities::getObj(file.path(dirThs,"rda_AllometryData.RData")) |>
           dplyr::rename(obs=w) |> tibble::rownames_to_column(var="obs_id") |>
           dplyr::mutate(dplyr::across(dplyr::everything(),as.character));
plt = ggplot(dfrZW,aes(x=as.numeric(z),y=as.numeric(obs),colour=x)) + geom_point();
print(plt);
##--remove some clearly bad values
dfrZWp = dfrZW |>
          dplyr::filter(!((as.numeric(z)>125)&(as.numeric(obs)<100))) |>
          dplyr::filter(!((as.numeric(z)<100)&(as.numeric(obs)>400))) |>
          dplyr::filter(!((as.numeric(z)< 40)&(as.numeric(obs)> 75))) |>
          dplyr::filter(obs_id!=3667) |>
          dplyr::filter(as.numeric(obs)>0);
plt = ggplot(dfrZWp,aes(x=as.numeric(z),y=as.numeric(obs),colour=x)) + geom_point() +
        coord_cartesian(xlim=c(0,60),ylim=c(0,100))
print(plt);
dfrZW = dfrZWp; rm(dfrZWp);

##--map functions/parameters to data----
dfrZWppp = dfrZW |>
             dplyr::filter(dplyr::between(as.numeric(y),2015,2024)) |>  #--keep only data w/in model time frame5
             dplyr::arrange(obs_id) |>
             dplyr::mutate(m=ifelse(m=="all",NA,m),     #--convert maturity state "all" to NA
                           p=ifelse(x=="female",NA,p));#--convert post-molt age for females to NA

#--parameters
##--TYPE=="pwrLaw1" (w = pA * (RTMB::AD(z)^pB);)
params1 = list(pA=NA,pB=NA);
params1[["pA"]] = c(0.00027,0.000562,0.000441);  #--male, imm female,mat female
params1[["pB"]] = c(3.022134,2.8169282,2.898686);#--male, imm female,mat female
##--TYPE=="pwrLaw2" (w = exp(pLnA+pLnS+pB*ln(RTMB::AD(z)/pZ0));)
params2 = list(pLnA=NA,pB=NA);
pLnS    = log(1.0); #--same scale as pA
pZ0     = 80;   #--reference size (mm CW) for individual crab weight
params2[["pB"]]   = params1[["pB"]];#--male, imm female,mat female
params2[["pLnA"]] = log(c(0.00027,0.000562,0.000441))-pLnS + params2[["pB"]]*log(pZ0);  #--male, imm female,mat female
##--check conversion
exp(params2[["pLnA"]]);
params1[["pA"]]*(pZ0)^(params1[["pB"]]);

#--source objective function----
source(file.path(dirThs,"objfun_EstimateAllometryAlt.R"));

#--run objective function in R environment for testing----
if (FALSE) {
  testing=TRUE;
  TYPE = "pwrLaw2";
  if (TYPE=="pwrLaw1") params = params1;
  if (TYPE=="pwrLaw2") params = params2;
  NLLTYPE = "dlnorm";
  objfun<-objfun_EstimateAllometryAlt(params);
  objfun$nll;
  sum(objfun$dfrZWs$nll);
  ggplot(objfun$dfrZWs,aes(x=prd,y=obs,colour=x)) + geom_point() + geom_abline(slope=1);
}

#--run in RTMB environment----
if (FALSE){
  require(RTMB);
  testing=FALSE;
  TYPE = "pwrLaw2";
  if (TYPE=="pwrLaw1") params = params1;
  if (TYPE=="pwrLaw2") params = params2;
  NLLTYPE = "dnorm";
  objfun<-MakeADFun(objfun_EstimateAllometryAlt,
                    params,
                    random=NULL,
                    map=list(),
                    hessian=TRUE,
                    LaplaceNonZeroGradient=TRUE,
                    ADreport=FALSE,
                    ridge.correct=FALSE, #--Note: apparently EXPERIMENTAL
                    silent=FALSE);
  opt <- nlminb(objfun$par, objfun$fn, objfun$gr,objfun$he);
  sdrep <- sdreport(objfun);
  print(sdrep);
  if (TYPE=="pwrLaw2"){
    est = sdrep$value;
    sdr = sdrep$sd;
    cat("pA: \n  est      sd\n");
    for (i in 1:length(est)){
      cat(est[i],sdr[i],"\n");
    }
  }

  par_best = objfun$env$last.par.best;
  par_best;
  objfun$fn(par_best);
  objfun$gr(par_best);
  objfun$he(par_best);
  rep = objfun$report(par_best);
  dfrZWs = rep$dfrZWs;
  plt = ggplot() +
          geom_point(data=dfrZWs,mapping=aes(x=obs,y=prd,colour=x,fill=x)) +
          geom_abline(slope=1) + facet_grid(x~.);
  print(plt);
  plt = ggplot() +
          geom_point(data=dfrZWs,mapping=aes(x=z,y=nll,colour=x,fill=x)) + facet_grid(x~.);
  print(plt);
}

