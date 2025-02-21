#--test Growth_PrGr specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfo_Growth_PrGr.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType2.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Growth_PrGr.R"))
  source(file.path(dirPrj,"R","calcGrowth_PrGr.R"))
}


##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "data-vertical";
if (type=="data-vertical"){
  ###--growth with data-vertical----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth_PrGr/inputSpecs_Growth_PrGr.data-vertical.txt");
  res   = readParamInfo_Growth_PrGr(conn,TRUE);
  lstPrGr = extractParamInfo_Growth_PrGr(res,dims,FALSE);
  params = list(pPrGr_FPs=lstPrGr$params);#--"FP" for "fixed" parameters
} else
if (type=="data-horizontal"){
  ###--growth with data-horizontal----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth_PrGr/inputSpecs_Growth_PrGr.data-horizontal.txt");
  res   = readParamInfo_Growth_PrGr(conn,TRUE);
  lstPrGr = extractParamInfo_Growth_PrGr(res,dims);
  params = list(pPrGr_FPs=lstPrGr$params);#--"FP" for "fixed" parameters
} else
if (type=="function"){
  ###--growth with function----
  dims = setupModelDims(zcs=seq(55.5,104.5,5));
  conn  = file.path(dirPrj,"testing/testGrowth_PrGr/inputSpecs_Growth_PrGr.function.txt");
  res   = readParamInfo_Growth_PrGr(conn,FALSE);
  lstPrGr = extractParamInfo_Growth_PrGr(res,dims,FALSE);
  params = list(pPrGr_MPs=lstPrGr$MPs$params);
  if (!is.null(lstPrGr$OPs$params)) params[["pPrGr_OPs"]]=lstPrGr$OPs$params;
  if (!is.null(lstPrGr$OPs$params)) params[["pPrGr_DPs"]]=lstPrGr$DPs$params;
  if (!is.null(lstPrGr$OPs$params)) params[["pPrGr_REs"]]=lstPrGr$REs$params;
}
inputs = list();
inputs$dims  = dims;
inputs$lstPrGr = lstPrGr;#--add lstPrGr to inputs

#--test probability of maturing function----
source(file.path(dirPrj,"R/calcGrowth_PrGr.R"));
PrGr = calcGrowth_PrGr(inputs$dims,inputs$lstPrGr,params,TRUE);

#--test PrGr in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  REPORT(dims);

  #--calculate probability of maturing----
  info = inputs$lstPrGr;
  PrGr = calcGrowth_PrGr(dims,info,params,verbose);
  REPORT(PrGr);

  #--calculate objective function value
  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=list(),silent=FALSE);
rep = obj$report();
tmPrGr = rep$PrGr[["2020"]][["1"]];
rownames(tmPrGr) = as.character(1:dims$nCs);
colnames(tmPrGr) = as.character(1:dims$nCs);
dfrPrGr=tibble::as_tibble(tmPrGr) |>
         dplyr::mutate(row_id=dplyr::row_number()) |>
         tidyr::pivot_longer(-row_id,names_to="col_id",values_to="value") |>
         dplyr::mutate(col_id=as.numeric(col_id),
                       value=ifelse(value>0,value,NA));
View(dfrPrGr)
ggplot(dfrPrGr,aes(x=col_id,y=row_id,fill=value)) +
  geom_raster() +
  labs(x="size (mm CW)",y="probability of moting",shape="post-molt age") +
  wtsPlots::getStdTheme()

