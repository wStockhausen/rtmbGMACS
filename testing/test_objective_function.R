#--test the gmacs objective function
require(RTMB);
require(ggplot2);
require(rtmbGMACS);

#--set up inputs list----
inputs = list(testing=TRUE,verbose=TRUE);

##--set up model dimensions----
#----"y","s","r","x","m","a","p","z"
dims = list();
dims$y_ = 2020:2024;                           #--years
dims$s_ = 1:4;                                 #--seasons
dims$r_ = "EBS";                               #--regions
dims$x_ = c("male","female");                  #--sexes
dims$m_ = c("immature","mature");              #--maturity states
dims$p_ = c("new shell","old shell");          #--post-molt ages
dims$zc = seq(25,185,5); n_zcs = length(zc);   #--size bin cutpts
dims$zb =  0.5*(zc[2:n_zcs]+zc[1:(n_zcs-1)]);  #--size bin centers
dims$f_ = c("TCF","NMFS");                     #--fleets
for (dim in names(dims)){
  dimp = dims[[dim]];
  names(dimp) = as.character(dimp);
  dims[[dim]] = dimp;
}

dims$dmsYS   = createSparseDimsMap(y=dims$y_,s=dims$s_);
dims$dmsN    = createSparseDimsMap(x=dims$x_,m=dims$m_,p=dims$p_,z=dims$zb);
dims$dmsFlts = createSparseDimsMap(f=dims$f_);

inputs$dims = dims;

##--set up options----
options = list();
options$initN = "noPop"; #--add options: "unfished", "equilibrium", "estimated"
inputs$options = options;

##--set up data----
data = list();
inputs$data = data;

##--set up information regarding model functions and parameters----
param_info = list();
tblPrcs = tibble::tibble();
###--initial abundance----
if (options$initN=="noPop"){
  param_info[["initN"]] = NULL;#--no parameters because option = "noPop"
} else
if (options$initN=="unfished"){
  stop("not implemented yet.")
} else
if (options$initN=="equilibrium"){
  stop("not implemented yet.")
} else
if (options$initN=="estimated"){
  #--reference state
  idxRef = getSubset.DimsMap(dmsN,x=dims$x_["male"],m=dims$m_["mature"],
                             p=dims$p_["new shell"],z=dims$zb["102.5"])[[1]];
  pNdevs = vector("numeric",nrow(dmsN));
} else

###--recruitment----
#####--pRec_lnRbar: ln-scale mean recruitment, millions----
pRec_lnRbar = dplyr::bind_rows(
                tibble::tibble(name="pRec_lnRbar",id=1,init=log(5),lower=log(2),upper=log(7),phs=1,RE=FALSE,prior=NA,p1=NA,p2=NA),
                tibble::tibble(name="pRec_lnRbar",id=2,init=log(3),lower=log(2),upper=log(7),phs=1,RE=FALSE,prior=NA,p1=NA,p2=NA)
              );
####--rec devs----
pRec_lnRdvs = dplyr::bind_rows(
                 tibble::tibble(name="pRec_lnRdvs",id=1,init=0,lower=-5.0,upper=5.0,phs=1,RE=TRUE,prior=NA,p1=NA,p2=NA),
                 tibble::tibble(name="pRec_lnRdvs",id=2,init=0,lower=-5.0,upper=5.0,phs=1,RE=TRUE,prior=NA,p1=NA,p2=NA)
               );

####--sex ratio----
fcn = "logistic";
#--mean
pRec_LgtSexRatio = tibble::tibble(name="pRec_SexRatio",id=1,init=logit(0.5),lower=logit(0.4),upper=logit(0.6),phs=1,RE=FALSE,prior=NA,p1=NA,p2=NA);
####--recruitment size distribution----
fcn = "gamma"; #--gamma truncated at lowest size bin cutpoint
#--mean = shape*scale =  shape/rate
pRec_MnZ = tibble::tibble(name="pRec_MnZ",id=1,init=40,lower=30,upper=50,phs=1,RE=FALSE,prior=NA,p1=NA,p2=NA);
#--sdv  = mean*scale  =  mean/rate
pRec_SdZ = tibble::tibble(name="pRec_MnZ",id=1,init=40,lower=30,upper=50,phs=1,RE=FALSE,prior=NA,p1=NA,p2=NA);

inputs$param_info = params; rm(params);

###--natural mortality----
###--growth----
####--molting----
####--growth|molt----
####--combined----
###--maturity
options$terminal_molt = TRUE;
###--selectivity----
###--fishing mortality----
###--survey catchability----
###--ageing
###--movement


#--define parameters list
param_list = list();
if (inputs$testing) {
  param_list[["dummy"]] = 0.0;
} else {

}

#--make objective function
obj = MakeADFun(obj_fun,param_list,random=NULL,map=list(),silent=FALSE);

