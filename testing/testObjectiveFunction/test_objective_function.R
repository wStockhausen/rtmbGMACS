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
zc = seq(25,185,5); n_zcs = length(dims$zc);
dims$zc = zc;                                  #--size bin cutpts
dims$zb =  0.5*(zc[2:nzcs]+zc[1:(nzcs-1)]);    #--size bin centers
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
###--options: initial N----
options$initN = "zeroPop";  #--add options: "unfished", "equilibrium", "estimated"
options$initN = "inputPop"; #--add options: "unfished", "equilibrium", "estimated"
###--options: Recruitment----
#opts_rec = ??
###--options: Natural Mortality----
#opts_nm = ??
###--options: Growth----
#opts_grw = ??
###--options: Selectivity/Retention----
#opts_sel = ??
###--options: Fishery characteristics----
#opts_fsh = ??
###--options: Survey characteristics----
#opts_srv = ??
###--options: Movement----
#opts_mov = ??

###--non-data likelihood components----

###--add to `inputs`
inputs$options = options;

##--set up data----
data = list();
inputs$data = data;

##--set up information regarding model functions and parameters----
param_info = list();
tblPrcs = tibble::tibble();
###--initial abundance----
if (options$initN=="zeroPop"){
  param_info[["pInitN"]] = vector("numeric",nrow(dmsN));#--no parameters because option = "zeroPop"
} else
if (options$initN=="inputPop"){
  param_info[["pInitN"]] = NULL;#--no parameters because option = "zeroPop"
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
} else {
  stop(paste0("Initial N (options$initN) option '",options$initN,"' not recognized."));
}

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
opts_maturity = list()
opts_maturity$postmolt      = TRUE;#--maturity determined by postmolt size
opts_maturity$terminal_molt = TRUE;#--undergoes terminal molt
options$maturity = opts_maturity;
###--selectivity----
###--fishing mortality----
###--survey catchability----
###--ageing
###--movement

inputs$options = options;

#--define parameters list
param_list = list();
if (inputs$testing) {
  param_list[["pDummy"]] = 0.0;
}
param_list[["pInitN"]] =

#--make objective function
##--`inputs` list is implicit
obj = MakeADFun(obj_fun,param_list,random=NULL,map=list(),silent=FALSE);

