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
}

#--set up inputs list----
inputs = list(testing=TRUE,verbose=TRUE);

##--set up model dimensions----
#----"y","s","r","x","m","a","p","z"
dims = list();
dims$y_ = 1975:2024;                           #--years
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

dims$dms_ys     = createSparseDimsMap(y=dims$y_,s=dims$s_);
dims$dms_yrxmsz = createSparseDimsMap(y=dims$y_,r=dims$r_,x=dims$x_,m=dims$m_,p=dims$p_,z=dims$zb);
dims$dms_f      = createSparseDimsMap(f=dims$f_);

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
# fn  ="inputSpecs_InitialN.txt";
# params_info$initN = readParamInfo_InitialN(fn,TRUE);
# params_list$initN = extractParameters_InitialN(inputs$dims,params_info$initN);
# parameters_ = assignParameters_IniitialN(parameters_,params_list$initN);
# params_map_ = mapParameters_InitialN(params_map_,params_list$initN);

###--recruitment----
# fn  ="inputRecruitmentSpecifications.txt";
# param_info$recruitment = readParamInfo_Recruitment(fn,TRUE);
# params_rec = extractParameters_Recruitment(inputs$dims,param_info$recruitment);

###--natural mortality----
###--growth & ageing----
####--molting----
####--growth|molt----
####--combined----
###--maturity
###--selectivity----
fn  ="inputSpecs_SelFcns.txt";
params_info$initN = readParamInfo_SelFcns(fn,TRUE);
params_list$initN = extractParameters_SelFcns(inputs$dims,params_info$initN);
parameters_ = assignParameters_SelFcns(parameters_,params_list$initN);
params_map_ = mapParameters_SelFcns(params_map_,params_list$initN);
###--fishing mortality----
###--survey catchability----
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

