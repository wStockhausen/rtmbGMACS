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
  source(file.path(dirPrj,"R","SelectivityFunctions1.R"))
  #source(file.path(dirPrj,"R","calcSelectivity.R"))
}

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "data-vertical";
if (type=="data-vertical"){
  ###--recruitment time series with data-vertical----
  dims = setupModelDims(zcs=seq(25.5,80.5,5));
  conn   = file.path(dirPrj,"testing/testSelectivity/inputSpecs_Selectivity.data-vertical.txt");
  res    = readParamInfo_Selectivity(conn,TRUE);
  lstSel = extractParamInfo_Selectivity(res,dims,FALSE);
  params = list(pSel_FPs=lstSel$params);#--"FP" for "fixed" parameters
  map = lstSel$map;
} else
if (type=="data-horizontal"){
  ###--recruitment time series with data-horizontal----
  dims   = setupModelDims(zcs=seq(24.5,184.5,5));
  conn   = file.path(dirPrj,"testing/testSelectivity/inputSpecs_Selectivity.data-horizontal.txt");
  res    = readParamInfo_Selectivity(conn,TRUE);
  lstSel = extractParamInfo_Selectivity(res,dims);
  params = list(pSel_FPs=lstSel$params);#--"FP" for "fixed" parameters
  map = lstSel$map;
} else
if (type=="function"){
  ###--recruitment time series with function----
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
inputs          = list();
inputs$dims     = dims;
inputs$lstSel = lstSel;#--add lstSel to inputs

#--test calcSelectivity function----
source(file.path(dirPrj,"R/calcSelectivity.R"));
lstSelVals = calcSelectivity(inputs$dims,inputs$lstSel,params,TRUE);
dfrSel = tibble::as_tibble(arrSel[1:5,1,1]) |>
         dplyr::mutate(row_id=dplyr::row_number()) |>
         dplyr::mutate(y=as.numeric(dims$y[row_id]));
ggplot(dfrSel, aes(x=y,y=value)) +
  geom_line() + geom_hline(yintercept=0,linetype=3,colour="white") +
  labs(x="year",y="selectivity functions") +
  wtsPlots::getStdTheme()

#--test arrSel in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate recruitment time seriess----
  info = inputs$lstSel;
  lstSel = calcSelectivity(dims,info,params,verbose);
  REPORT(arrSel);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

verbose=FALSE;
params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=map,silent=FALSE);
rep = obj$report();
sum(rep$arrSel[1,1,]);#--should be = 1
dfrSel = tibble::as_tibble(rep$arrSel[1,1,]) |>
         dplyr::mutate(row_id=dplyr::row_number()) |>
         dplyr::mutate(r=dims$dmsC$r[row_id],
                       x=dims$dmsC$x[row_id],
                       m=dims$dmsC$m[row_id],
                       p=dims$dmsC$p[row_id],
                       z=as.numeric(as.character(dims$dmsC$z[row_id])),
                       mp=paste(m,p),
                       rxmp_from=paste(r,x,m,p));
ggplot(dfrSel, aes(x=z,y=value)) +
  geom_line() + geom_hline(yintercept=0,linetype=3,colour="white") +
  facet_grid(x~mp) +
  labs(x="size (mm CW)",y="selectivity functions") +
  wtsPlots::getStdTheme()
