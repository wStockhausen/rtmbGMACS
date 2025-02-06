#--test the gmacs objective function
require(RTMB);
require(ggplot2);
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path)
if (FALSE){
  require(rtmbGMACS);
} else {
  source(file.path(dirPrj,"R","DimensionsFunctions.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Dataframe.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R","readParamInfo_InitialN.R"))
  source(file.path(dirPrj,"R","extractParameters_InitialN.R"))
}

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
##--set up data----
###--data likelihood components----
###--non-data likelihood components----

##--set up information regarding model functions and parameters----
params_info = list(); #--list with input parameter information
params_list = list(); #--list with extracted parameter information
parameters_ = list(); #--list defining parameters (MakeADFun input `parameters`)
params_map_ = list(); #--list defining how to collect and fix parameters (MakeADFun input `map`)
REs          = vector("character",length=0); #--character vector identifying random effect parameters (MakeADFun input `random`)

###--initial abundance----
fn  ="inputSpecs_InitialN.txt";
params_info$initN = readParamInfo_InitialN(fn,TRUE);
params_list$initN = extractParameters_InitialN(inputs$dims,params_info$initN);
parameters_ = assignParameters_IniitialN(parameters_,params_list$initN);
params_map_ = mapParameters_InitialN(params_map_,params_list$initN);

###--recruitment----
# fn  ="inputRecruitmentSpecifications.txt";
# param_info$recruitment = readParamInfo_Recruitment(fn,TRUE);
# params_rec = extractParameters_Recruitment(inputs$dims,param_info$recruitment);

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
source(file.path(dirThs,"obj_fun_iniitialN.R"))
##--`inputs` list is implicit
obj = MakeADFun(obj_fun,param_list,random=NULL,map=list(),silent=FALSE);

