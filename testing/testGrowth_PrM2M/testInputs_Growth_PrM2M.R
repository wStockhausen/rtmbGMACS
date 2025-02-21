#--test Growth_PrM2M specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","MiscFunctions.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Dataframe.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Transforms.R"))
  source(file.path(dirPrj,"R","readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R","readParamInfo_Growth_PrM2M.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType2.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Growth_PrM2M.R"))
  source(file.path(dirPrj,"R","calcGrowth_PrM2M.R"))
}


##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "data-vertical";
if (type=="data-vertical"){
  ###--growth with data-vertical----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth_PrM2M/inputSpecs_Growth_PrM2M.data-vertical.txt");
  res   = readParamInfo_Growth_PrM2M(conn,TRUE);
  lstPrM2M = extractParamInfo_Growth_PrM2M(res,dims,FALSE);
  params = list(pPrM2M_FPs=lstPrM2M$params);#--"FP" for "fixed" parameters
} else
if (type=="data-horizontal"){
  ###--growth with data-horizontal----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth_PrM2M/inputSpecs_Growth_PrM2M.data-horizontal.txt");
  res   = readParamInfo_Growth_PrM2M(conn,TRUE);
  lstPrM2M = extractParamInfo_Growth_PrM2M(res,dims);
  params = list(pPrM2M_FPs=lstPrM2M$params);#--"FP" for "fixed" parameters
} else
if (type=="function"){
  ###--growth with function----
  dims = setupModelDims(zcs=seq(55.5,104.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth_PrM2M/inputSpecs_Growth_PrM2M.function.txt");
  res   = readParamInfo_Growth_PrM2M(conn,FALSE);
  lstPrM2M = extractParamInfo_Growth_PrM2M(res,dims,FALSE);
  params = list(pPrM2M_MPs=lstPrM2M$MPs$params);
  if (!is.null(lstPrM2M$OPs$params)) params[["pPrM2M_OPs"]]=lstPrM2M$OPs$params;
  if (!is.null(lstPrM2M$OPs$params)) params[["pPrM2M_DPs"]]=lstPrM2M$DPs$params;
  if (!is.null(lstPrM2M$OPs$params)) params[["pPrM2M_REs"]]=lstPrM2M$REs$params;
}
inputs = list();
inputs$dims  = dims;
inputs$lstPrM2M = lstPrM2M;#--add lstPrM2M to inputs

#--test probability of maturing function----
source(file.path(dirPrj,"R/calcGrowth_PrM2M.R"));
prM2M = calcGrowth_PrM2M(inputs$dims,inputs$lstPrM2M,params,TRUE);

#--test prM2M in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  REPORT(dims);

  #--calculate probability of maturing----
  info = inputs$lstPrM2M;
  prM2M = calcGrowth_PrM2M(dims,info,params,verbose);
  REPORT(prM2M);

  #--calculate objective function value
  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=list(),silent=FALSE);
rep = obj$report();
dfrPrM2M = dplyr::bind_cols(rep$dims$dmsC,prM=rep$prM2M$`2020`$`1`) |>
             dplyr::mutate(z=as.numeric(as.character(z)),
                           category=paste(x,m,p));
View(dfrPrM2M)
ggplot(dfrPrM2M,aes(x=z,y=prM,colour=category,shape=p)) +
  geom_line() + geom_point() +
  labs(x="size (mm CW)",y="probability of moting",shape="post-molt age") +
  wtsPlots::getStdTheme()

