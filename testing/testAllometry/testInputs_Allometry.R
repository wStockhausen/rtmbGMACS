#--test the gmacs objective function
require(RTMB);
require(ggplot2);
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
setwd(dirThs);
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
  source(file.path(dirPrj,"R","extractParamInfo_Allometry.R"))
}

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirThs,"..","r_setupModelDimensions.TestA.R"))
dims = setupModelDims();

###--allometry----
conn = file.path("./inputSpecs_Allometry.data-vertical.txt");
res      = readParamInfo_Allometry(conn,TRUE);
lstAllom = extractParamInfo_Allometry(res,dims$dmsYSN);

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
conn = file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.data.txt");
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
}

#--make objective function
source(file.path(dirThs,"obj_fun_Allometry.R"))
##--`inputs` list is implicit
obj = MakeADFun(obj_fun,param_list,random=NULL,map=list(),silent=FALSE);

