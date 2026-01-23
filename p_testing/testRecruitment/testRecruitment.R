#--test the gmacs objective function
require(RTMB);
require(ggplot2);
require(rtmbGMACS);

#--set up inputs list----
inputs = list(testing=TRUE,verbose=TRUE);

##--set up model dimensions----
#----"y","s","r","x","m","a","p","z"
dims = list();
dims$y = 2020:2024;                           #--years
dims$s = 1:4;                                 #--seasons
dims$r = "EBS";                               #--regions
dims$x = c("male","female");                  #--sexes
dims$m = c("immature","mature");              #--maturity states
dims$a = c("new_shell","old_shell");          #--post-molt ages
dims$zc = seq(25,185,5);                       #--size bin cutpts
zc = dims$zc; nzcs = length(zc);
dims$zb =  0.5*(zc[2:nzcs]+zc[1:(nzcs-1)]);  #--size bin centers
dims$f = c("TCF","NMFS");                     #--fleets
for (dim in names(dims)){
  dimp = dims[[dim]];
  names(dimp) = as.character(dimp);
  dims[[dim]] = dimp;
}

dims$dmsYS   = createSparseDimsMap(y=dims$y,s=dims$s);
dims$dmsN    = createSparseDimsMap(r=dims$r,x=dims$x,m=dims$m,a=dims$a,z=dims$zb);
dims$dmsFlts = createSparseDimsMap(f=dims$f);

inputs$dims = dims;

##--set up options----
options = list();
###--options: initial N----
#opts_InitN = ??
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

###--add to `inputs`
inputs$options = options;

##--set up data----
data = list();
inputs$data = data;

###--data likelihood components----
###--non-data likelihood components----

##--set up information regarding model functions and parameters----
param_info = list();
###--initial abundance----
###--recruitment----
fn  ="inputRecruitmentSpecifications.txt";
param_info$recruitment = readParamInfo_Recruitment(fn,TRUE);
params_rec = extractRecruitmentParameters(inputs$dims,param_info$recruitment);
###--natural mortality----
###--growth----
####--molting----
####--growth|molt----
####--combined----
###--maturity
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
##--`inputs` list is implicit
obj = MakeADFun(obj_fun,param_list,random=NULL,map=list(),silent=FALSE);

