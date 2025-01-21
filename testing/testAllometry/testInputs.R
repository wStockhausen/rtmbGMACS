#--test the gmacs objective function
require(RTMB);
require(ggplot2);
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path)
if (FALSE){
  require(rtmbGMACS);
} else {
  source(file.path(dirPrj,"R","DimensionsFunctions.R"))
  source(file.path(dirPrj,"R","DimensionsUtilities.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Dataframe.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Transforms.R"))
  source(file.path(dirPrj,"R","readParamInfo_Allometry.R"))
  source(file.path(dirPrj,"R","extractParameters_Allometry.R"))
}

#--set up inputs list----
inputs = list(testing=TRUE,verbose=TRUE);

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
dims = list();
dims$y_ = 1975:2024;                           #--years
dims$s_ = 1:4;                                 #--seasons
dims$r_ = "EBS";                               #--regions
dims$x_ = c("male","female");                  #--sexes
dims$m_ = c("immature","mature");              #--maturity states
dims$a_ = c("new_shell","old_shell");          #--post-molt ages
dims$zc = seq(24.5,184.5,5);
zc = dims$zc; n_zcs = length(zc);              #--size bin cutpts
dims$zb =  0.5*(zc[2:n_zcs]+zc[1:(n_zcs-1)]);  #--size bin centers
dims$f_ = c("TCF","NMFS");                     #--fleets
for (dim in names(dims)){
  dimp = dims[[dim]];
  names(dimp) = as.character(dimp);
  dims[[dim]] = dimp;
}

dims$dms_ys     = createSparseDimsMap(y=dims$y_,s=dims$s_);
dims$dms_yrxmaz = createSparseDimsMap(y=dims$y_,r=dims$r_,x=dims$x_,m=dims$m_,a=dims$a_,z=dims$zb);
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

###--allometry----
conn = file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.txt");
res  = readParamInfo_Allometry(conn,TRUE);
lstAllom = extractParameters_Allometry(res,dims$dms_yrxmaz);
###--initial abundance----
###--recruitment----
###--natural mortality----
###--growth & ageing----
####--molting----
####--growth|molt----
####--combined----
###--maturity
###--selectivity----
###--fishing mortality----
###--survey catchability----
###--movement

#--define parameters list
par_list = list();
if (inputs$testing) {
  param_list[["dummy"]] = 0.0;
} else {

}

#--make objective function
source(file.path(dirThs,"obj_fun_Allometry.R"))
##--`inputs` list is implicit
obj = MakeADFun(obj_fun,param_list,random=NULL,map=list(),silent=FALSE);

