#--test selectivity specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","MiscFunctions.R"))
  source(file.path(dirPrj,"R","readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R","readParamInfo_Selectivity.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Selectivity.R"))
  source(file.path(dirPrj,"R","Functions_Selectivity.R"))
  source(file.path(dirPrj,"R","calcSelectivity.R"))
}

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "pre-specified-vertical";
if (type=="pre-specified-vertical"){
  ###--selectivity functions with data-vertical type input----
  dims = setupModelDims(zcs=seq(25.5,80.5,5));
  conn   = file.path(dirPrj,"testing/testSelectivity/inputSpecs_Selectivity.pre-specified-vertical.txt");
  res    = readParamInfo_Selectivity(conn,TRUE);
  lstSel = extractParamInfo_Selectivity(res,dims,FALSE);
  params = list(pSel_FPs=lstSel$params);#--"FP" for "fixed" parameters
  map = lstSel$map;
} else
if (type=="pre-specified-horizontal"){
  ###--selectivity functions with data-horizontal type input----
  dims   = setupModelDims(zcs=seq(24.5,184.5,5));
  conn   = file.path(dirPrj,"testing/testSelectivity/inputSpecs_Selectivity.pre-specified-horizontal.txt");
  res    = readParamInfo_Selectivity(conn,TRUE);
  lstSel = extractParamInfo_Selectivity(res,dims);
  params = list(pSel_FPs=lstSel$params);#--"FP" for "fixed" parameters
  map = lstSel$map;
} else
if (type=="function"){
  ###--selectivity functions with function type input----
  dims = setupModelDims(zcs=seq(25.5,84.5,5));
  conn   = file.path(dirPrj,"testing/testSelectivity/inputSpecs_Selectivity.function.txt");
  res    = readParamInfo_Selectivity(conn,FALSE);
  lstSel = extractParamInfo_Selectivity(res,dims,FALSE);
  params = list(pSel_MPs=lstSel$MPs$params);
  map = lstSel$map;
  if (!is.null(lstSel$OPs$params)) params[["pSel_OPs"]]=lstSel$OPs$params;
  if (!is.null(lstSel$DPs$params)) params[["pSel_DPs"]]=lstSel$DPs$params;
  if (!is.null(lstSel$REs$params)) params[["pSel_REs"]]=lstSel$REs$params;
}
inputs        = list();
inputs$dims   = dims;
inputs$lstSel = lstSel;#--add lstSel to inputs

#--test calcSelectivity function----
source(file.path(dirPrj,"R/calcSelectivity.R"));
lstSelVals = calcSelectivity(inputs$dims,inputs$lstSel,params,FALSE);
plotSelVals<-function(if_=1,iy_=1,is_=1){
  dfrSel = dims$dmsC |>
            dplyr::mutate(value=(lstSelVals[[if_]])[iy_,is_,],
                          z=as.numeric(as.character(z)));
  ggplot(dfrSel, aes(x=as.numeric(z),y=value)) +
    geom_line() + geom_hline(yintercept=0,linetype=3,colour="white") +
    facet_grid(m+p~x) +
    labs(x="size",y="selectivity functions") +
    wtsPlots::getStdTheme();
}
plotSelVals(1,1,1)
plotSelVals(2,1,1)
plotSelVals(1,4,1)
plotSelVals(2,4,1)
plotSelVals(3,2,1)

#--test selectivity functions in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate selectivity functions----
  info = inputs$lstSel;
  lstSel = calcSelectivity(dims,info,params,verbose);
  REPORT(lstSel);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

verbose=FALSE;
params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=map,silent=FALSE);
rep = obj$report();
lstSelVals = rep$lstSel;
plotSelVals(1,1,1)
plotSelVals(2,1,1)
plotSelVals(1,4,1)
plotSelVals(2,4,1)

