#--test Growth_PrMolt specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfo_Growth_PrMolt.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Growth_PrMolt.R"))
  source(file.path(dirPrj,"R","calcGrowth_PrMolt.R"))
}

inputs = list();

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "data-vertical";
if (type=="data-vertical"){
  ###--growth with data-vertical----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth_PrMolt/inputSpecs_Growth_PrMolt.data-vertical.txt");
  res   = readParamInfo_Growth_PrMolt(conn,TRUE);
  lstPrMolt = extractParamInfo_Growth_PrMolt(res,dims,FALSE);
  params = list(pPrMolt_FPs=lstPrMolt$params);#--"FP" for "fixed" parameters
} else
if (type=="data-horizontal"){
  ###--probability of molting with data-horizontal----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth_PrMolt/inputSpecs_Growth_PrMolt.data-horizontal.txt");
  res   = readParamInfo_Growth_PrMolt(conn,TRUE);
  lstPrMolt = extractParamInfo_Growth_PrMolt(res,dims);
  params = list(pPrMolt_FPs=lstPrMolt$params);#--"FP" for "fixed" parameters
} else
if (type=="function"){
  ###--probability of molting with function----
  dims = setupModelDims(zcs=seq(55.5,104.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth_PrMolt/inputSpecs_Growth_PrMolt.function.txt");
  res   = readParamInfo_Growth_PrMolt(conn,FALSE);
  lstPrMolt = extractParamInfo_Growth_PrMolt(res,dims,FALSE);
  params = list(pPrMolt_MPs=lstPrMolt$MPs$params);
  if (!is.null(lstPrMolt$OPs$params)) params[["pPrMolt_OPs"]]=lstPrMolt$OPs$params;
  if (!is.null(lstPrMolt$DPs$params)) params[["pPrMolt_DPs"]]=lstPrMolt$DPs$params;
  if (!is.null(lstPrMolt$REs$params)) params[["pPrMolt_REs"]]=lstPrMolt$REs$params;
}
inputs$dims  = dims;
inputs$lstPrMolt = lstPrMolt;#--add lstPrMolt to inputs

#--test function for probability of molting----
source(file.path(dirPrj,"R/calcGrowth_PrMolt.R"));
prM = calcGrowth_PrMolt(inputs$dims,inputs$lstPrMolt,params,FALSE,loopIC_=TRUE); #--slower
prM = calcGrowth_PrMolt(inputs$dims,inputs$lstPrMolt,params,FALSE,loopIC_=FALSE);#--faster

#--test M in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  REPORT(dims);
  #--calculate weight-at-size----
  info = inputs$lstPrMolt;
  prM = calcGrowth_PrMolt(dims,info,params,verbose);
  REPORT(prM);

  #--calculate objective function value
  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

verbose=TRUE;
params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=list(),silent=FALSE);
rep = obj$report();
dfrPrM = dplyr::bind_cols(rep$dims$dmsC,prM=rep$prM[1,1,]) |>
         dplyr::mutate(z=as.numeric(as.character(z)),
                       category=paste(x,m,p));
View(dfrPrM)
ggplot(dfrPrM,aes(x=z,y=prM,colour=category,shape=p)) +
  geom_line() + geom_point() +
  labs(x="size (mm CW)",y="probability of moting",shape="post-molt age") +
  wtsPlots::getStdTheme()
