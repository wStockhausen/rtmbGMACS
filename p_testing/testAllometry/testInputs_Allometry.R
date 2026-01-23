#--test allometry specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R","readParamInfo_Allometry.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Allometry.R"))
  source(file.path(dirPrj,"R","calcAllometry.R"))
}

inputs = list();

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "data-vertical";
if (type=="data-vertical"){
  ###--allometry with data-vertical----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn     = file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.data-vertical.txt");
  res      = readParamInfo_Allometry(conn,TRUE);
  lstAllom = extractParamInfo_Allometry(res,dims,FALSE);
  params = list(pAllom_FPs=lstAllom$params);#--"FP" for "fixed" parameters
} else
if (type=="data-horizontal"){
  ###--allometry with data-horizontal----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn     = file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.data-horizontal.txt");
  res      = readParamInfo_Allometry(conn,TRUE);
  lstAllom = extractParamInfo_Allometry(res,dims);
  params = list(pAllom_FPs=lstAllom$params);#--"FP" for "fixed" parameters
} else
if (type=="function"){
  ###--allometry with function----
  dims = setupModelDims(zcs=seq(55.5,104.5,5));
  conn     = file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.function.txt");
  res      = readParamInfo_Allometry(conn,FALSE);
  lstAllom = extractParamInfo_Allometry(res,dims,FALSE);
  params = list(pAllom_MPs=lstAllom$MPs$params);
  if (!is.null(lstAllom$OPs$params)) params[["pAllom_OPs"]]=lstAllom$OPs$params;
  if (!is.null(lstAllom$OPs$params)) params[["pAllom_DPs"]]=lstAllom$DPs$params;
  if (!is.null(lstAllom$OPs$params)) params[["pAllom_REs"]]=lstAllom$REs$params;
}
inputs$dims     = dims;
inputs$lstAllom = lstAllom;#--add lstAllom to inputs

#--test wAtZ function----
source(file.path(dirPrj,"R/calcAllometry.R"));
wAtZ = calcAllometry(inputs$dims,inputs$lstAllom,params,FALSE,loopIC_=TRUE); #--slower
wAtZ = calcAllometry(inputs$dims,inputs$lstAllom,params,FALSE,loopIC_=FALSE);#--faster

#--test wAtZ in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate weight-at-size----
  info = inputs$lstAllom;
  wAtZ = calcAllometry(dims,info,params,verbose);
  REPORT(wAtZ);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=list(),silent=FALSE);

