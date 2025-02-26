#--test recruitment category proportions specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfo_Recruitment_TimeSeries.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Recruitment_TimeSeries.R"))
  source(file.path(dirPrj,"R","calcRecruitment_TimeSeries.R"))
}

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "data-vertical";
if (type=="data-vertical"){
  ###--recruitment time series with data-vertical----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn     = file.path(dirPrj,"testing/testRecruitment_TimeSeries/inputSpecs_Recruitment_TimeSeries.data-vertical.txt");
  res      = readParamInfo_Recruitment_TimeSeries(conn,TRUE);
  lstRecTS = extractParamInfo_Recruitment_TimeSeries(res,dims,FALSE);
  params = list(pRecTS_FPs=lstRecTS$params);#--"FP" for "fixed" parameters
} else
if (type=="data-horizontal"){
  ###--recruitment time series with data-horizontal----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn     = file.path(dirPrj,"testing/testRecruitment_TimeSeries/inputSpecs_Recruitment_TimeSeries.data-horizontal.txt");
  res      = readParamInfo_Recruitment_TimeSeries(conn,TRUE);
  lstRecTS = extractParamInfo_Recruitment_TimeSeries(res,dims);
  params = list(pRecTS_FPs=lstRecTS$params);#--"FP" for "fixed" parameters
} else
if (type=="function"){
  ###--recruitment time series with function----
  dims = setupModelDims(zcs=seq(25.5,84.5,5));
  conn     = file.path(dirPrj,"testing/testRecruitment_TimeSeries/inputSpecs_Recruitment_TimeSeries.function.txt");
  res      = readParamInfo_Recruitment_TimeSeries(conn,FALSE);
  lstRecTS = extractParamInfo_Recruitment_TimeSeries(res,dims,FALSE);
  params = list(pRecTS_MPs=lstRecTS$MPs$params);
  if (!is.null(lstRecTS$OPs$params)) params[["pRecTS_OPs"]]=lstRecTS$OPs$params;
  if (!is.null(lstRecTS$DPs$params)) params[["pRecTS_DPs"]]=lstRecTS$DPs$params;
  if (!is.null(lstRecTS$REs$params)) params[["pRecTS_REs"]]=lstRecTS$REs$params;
}
inputs          = list();
inputs$dims     = dims;
inputs$lstRecTS = lstRecTS;#--add lstRecTS to inputs

#--test arrRecTS function----
source(file.path(dirPrj,"R/calcRecruitment_TimeSeries.R"));
arrRecTS = calcRecruitment_TimeSeries(inputs$dims,inputs$lstRecTS,params,TRUE,loopIC_=TRUE); #--slower
arrRecTS = calcRecruitment_TimeSeries(inputs$dims,inputs$lstRecTS,params,FALSE,loopIC_=FALSE);#--faster
dfrRecTS = tibble::as_tibble(arrRecTS[1,1,]) |>
         dplyr::mutate(row_id=dplyr::row_number()) |>
         dplyr::mutate(r=dims$dmsC$r[row_id],
                       x=dims$dmsC$x[row_id],
                       m=dims$dmsC$m[row_id],
                       p=dims$dmsC$p[row_id],
                       z=as.numeric(as.character(dims$dmsC$z[row_id])),
                       mp=paste(m,p),
                       rxmp_from=paste(r,x,m,p));
ggplot(dfrRecTS, aes(x=z,y=value)) +
  geom_line() + geom_hline(yintercept=0,linetype=3,colour="white") +
  facet_grid(x~mp) +
  labs(x="size (mm CW)",y="recruitment time series") +
  wtsPlots::getStdTheme()

#--test arrRecTS in RTMB objective function----
map = list(pRecTS_MPs=factor(c(1,2,NA,NA),levels=c(1,2,3,4)));
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate recruitment time seriess----
  info = inputs$lstRecTS;
  arrRecTS = calcRecruitment_TimeSeries(dims,info,params,verbose,loopIC_=FALSE);
  REPORT(arrRecTS);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

verbose=FALSE;
params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=map,silent=FALSE);
rep = obj$report();
sum(rep$arrRecTS[1,1,]);#--should be = 1
dfrRecTS = tibble::as_tibble(rep$arrRecTS[1,1,]) |>
         dplyr::mutate(row_id=dplyr::row_number()) |>
         dplyr::mutate(r=dims$dmsC$r[row_id],
                       x=dims$dmsC$x[row_id],
                       m=dims$dmsC$m[row_id],
                       p=dims$dmsC$p[row_id],
                       z=as.numeric(as.character(dims$dmsC$z[row_id])),
                       mp=paste(m,p),
                       rxmp_from=paste(r,x,m,p));
ggplot(dfrRecTS, aes(x=z,y=value)) +
  geom_line() + geom_hline(yintercept=0,linetype=3,colour="white") +
  facet_grid(x~mp) +
  labs(x="size (mm CW)",y="recruitment category proportions") +
  wtsPlots::getStdTheme()
