#--test growth specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfo_Growth.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Growth.R"))
  source(file.path(dirPrj,"R","calcGrowth.R"))
}

inputs = list();

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "data-vertical";
if (type=="data-vertical"){
  ###--growth with data-vertical----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth/inputSpecs_Growth.data-vertical.txt");
  res   = readParamInfo_Growth(conn,TRUE);
  lstNM = extractParamInfo_Growth(res,dims,FALSE);
  params = list(pNM_FPs=lstNM$params);#--"FP" for "fixed" parameters
} else
if (type=="data-horizontal"){
  ###--growth with data-horizontal----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth/inputSpecs_Growth.data-horizontal.txt");
  res   = readParamInfo_Growth(conn,TRUE);
  lstNM = extractParamInfo_Growth(res,dims);
  params = list(pNM_FPs=lstNM$params);#--"FP" for "fixed" parameters
} else
if (type=="function"){
  ###--growth with function----
  dims = setupModelDims(zcs=seq(55.5,104.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth/inputSpecs_Growth.function.txt");
  res   = readParamInfo_Growth(conn,FALSE);
  lstNM = extractParamInfo_Growth(res,dims,FALSE);
  params = list(pNM_MPs=lstNM$MPs$params);
  if (!is.null(lstNM$OPs$params)) params[["pNM_OPs"]]=lstNM$OPs$params;
  if (!is.null(lstNM$OPs$params)) params[["pNM_DPs"]]=lstNM$DPs$params;
  if (!is.null(lstNM$OPs$params)) params[["pNM_REs"]]=lstNM$REs$params;
}
inputs$dims  = dims;
inputs$lstNM = lstNM;#--add lstNM to inputs

#--test wAtZ function----
source(file.path(dirPrj,"R/calcGrowth.R"));
M = calcGrowth(inputs$dims,inputs$lstNM,params,FALSE,loopIC_=TRUE); #--slower
M = calcGrowth(inputs$dims,inputs$lstNM,params,FALSE,loopIC_=FALSE);#--faster

#--test M in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate weight-at-size----
  info = inputs$lstNM;
  M = calcGrowth(dims,info,params,verbose);
  REPORT(M);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=list(),silent=FALSE);

