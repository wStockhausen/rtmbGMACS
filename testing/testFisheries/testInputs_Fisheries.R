#--test fisheries specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfo_Fisheries.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Fisheries.R"))
  source(file.path(dirPrj,"R","calcFisheries.R"))
}

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "data-vertical";
if (type=="data-vertical"){
  ###--recruitment time series with data-vertical----
  dims = setupModelDims(zcs=seq(25.5,80.5,5));
  conn   = file.path(dirPrj,"testing/testFisheries/inputSpecs_Fisheries.data-vertical.txt");
  res    = readParamInfo_Fisheries(conn,TRUE);
  lstFsh = extractParamInfo_Fisheries(res,dims,FALSE);
  params = list(pSel_FPs=lstFsh$params);#--"FP" for "fixed" parameters
  map = lstFsh$map;
} else
if (type=="data-horizontal"){
  ###--recruitment time series with data-horizontal----
  dims   = setupModelDims(zcs=seq(24.5,184.5,5));
  conn   = file.path(dirPrj,"testing/testFisheries/inputSpecs_Fisheries.data-horizontal.txt");
  res    = readParamInfo_Fisheries(conn,TRUE);
  lstFsh = extractParamInfo_Fisheries(res,dims);
  params = list(pSel_FPs=lstFsh$params);#--"FP" for "fixed" parameters
  map = lstFsh$map;
} else
if (type=="function"){
  ###--recruitment time series with function----
  dims = setupModelDims(zcs=seq(25.5,84.5,5));
  conn   = file.path(dirPrj,"testing/testFisheries/inputSpecs_Fisheries.function.txt");
  res    = readParamInfo_Fisheries(conn,FALSE);
  lstFsh = extractParamInfo_Fisheries(res,dims,FALSE);
  params = list(pSel_MPs=lstFsh$MPs$params);
  map = lstFsh$map;
  if (!is.null(lstFsh$OPs$params)) params[["pSel_OPs"]]=lstFsh$OPs$params;
  if (!is.null(lstFsh$DPs$params)) params[["pSel_DPs"]]=lstFsh$DPs$params;
  if (!is.null(lstFsh$REs$params)) params[["pSel_REs"]]=lstFsh$REs$params;
}
inputs          = list();
inputs$dims     = dims;
inputs$lstFsh = lstFsh;#--add lstFsh to inputs

#--test calcFisheries function----
source(file.path(dirPrj,"R/calcFisheries.R"));
lstFshVals = calcFisheries(inputs$dims,inputs$lstFsh,params,TRUE);
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
  #--calculate fisheries time series----
  info = inputs$lstFsh;
  lstFsh = calcFisheries(dims,info,params,verbose);
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
