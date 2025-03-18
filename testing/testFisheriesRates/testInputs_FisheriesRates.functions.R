#--test "function" fishery rates specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfo_FisheriesRates.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_FisheriesRates.R"))
  source(file.path(dirPrj,"R","calcFisheriesRates.R"))
}

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))
dims = setupModelDims(zcs=seq(25.5,84.5,5));

##--set up model inputs
inputs = list(dims=dims);
params = list();
map = list();

##--get selectivity information
source(file.path(dirPrj,"R","readParamInfo_Selectivity.R"))
source(file.path(dirPrj,"R","extractParamInfo_Selectivity.R"))
source(file.path(dirPrj,"R","calcSelectivity.R"))
conn   = file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FshSelectivity.function.txt");
res    = readParamInfo_Selectivity(conn,FALSE);
lstSel = extractParamInfo_Selectivity(res,dims,FALSE);
params[["pSel_MPs"]]=lstSel$MPs$params;
map = addList(map,lstSel$map);
if (!is.null(lstSel$OPs$params)) params[["pSel_OPs"]]=lstSel$OPs$params;
if (!is.null(lstSel$DPs$params)) params[["pSel_DPs"]]=lstSel$DPs$params;
if (!is.null(lstSel$REs$params)) params[["pSel_REs"]]=lstSel$REs$params;
inputs$lstSel = lstSel;#--add lstSel to inputs

##--get fishery capture rates information
source(file.path(dirPrj,"R","readParamInfo_FisheriesCaptureRates.R"))
source(file.path(dirPrj,"R","extractParamInfo_FisheriesCaptureRates.R"))
source(file.path(dirPrj,"R","calcFisheriesCaptureRates.R"))
conn    = file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FisheriesCaptureRates.function.txt");
res     = readParamInfo_FisheriesCaptureRates(conn,FALSE);
lstFCRs = extractParamInfo_FisheriesCaptureRates(res,dims,FALSE);
params[["pFCRs_MPs"]]=lstFCRs$MPs$params;
map = addList(map,lstFCRs$map);
if (!is.null(lstFCRs$OPs$params)) params[["pFCRs_OPs"]]=lstFCRs$OPs$params;
if (!is.null(lstFCRs$DPs$params)) params[["pFCRs_DPs"]]=lstFCRs$DPs$params;
if (!is.null(lstFCRs$REs$params)) params[["pFCRs_REs"]]=lstFCRs$REs$params;
inputs$lstFCRs = lstFCRs;#--add lstFCRs to inputs

##--get fishery handling mortality rates information
source(file.path(dirPrj,"R","readParamInfo_FisheriesHandlingMortalityRates.R"))
source(file.path(dirPrj,"R","extractParamInfo_FisheriesHandlingMortalityRates.R"))
source(file.path(dirPrj,"R","calcFisheriesHandlingMortalityRates.R"))
conn    = file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FisheriesHandlingMortalityRates.function.txt");
res     = readParamInfo_FisheriesHandlingMortalityRates(conn,FALSE);
lstHMRs = extractParamInfo_FisheriesHandlingMortalityRates(res,dims,FALSE);
params[["pHMRs_MPs"]]=lstHMRs$MPs$params;
map = addList(map,lstHMRs$map);
if (!is.null(lstHMRs$OPs$params)) params[["pHMRs_OPs"]]=lstHMRs$OPs$params;
if (!is.null(lstHMRs$DPs$params)) params[["pHMRs_DPs"]]=lstHMRs$DPs$params;
if (!is.null(lstHMRs$REs$params)) params[["pHMRs_REs"]]=lstHMRs$REs$params;
inputs$lstHMRs = lstHMRs;#--add lstHMRs to inputs

###--fishery rates inputs with function format----
conn   = file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FisheriesRates.txt");
res    = readInfo_FisheriesRates(conn,FALSE);
lstSrv = extractInfo_FisheriesRates(res,dims,FALSE);
params[["pFsh_MPs"]]=lstFsh$MPs$params;
map = addList(map,lstFsh$map);
if (!is.null(lstFsh$OPs$params)) params[["pFsh_OPs"]]=lstFsh$OPs$params;
if (!is.null(lstFsh$DPs$params)) params[["pFsh_DPs"]]=lstFsh$DPs$params;
if (!is.null(lstFsh$REs$params)) params[["pFsh_REs"]]=lstFsh$REs$params;
inputs$lstFsh = lstFsh;#--add lstFsh to inputs

#--test calcFisheriesRates function----
source(file.path(dirPrj,"R/calcFisheriesRates.R"));
#--inputs list elements: dims, lstFsh,lstSels,lstFCRs,lstHMRs
lstFshVals = calcFisheriesRates(inputs,params,TRUE);

#--plot results----
plotSrvVals<-function(if_=1,iy_=1,is_=1){
  dfrSrv = dims$dmsC |>
            dplyr::mutate(value=(lstSrvVals[[if_]])[iy_,is_,],
                          z=as.numeric(as.character(z)));
  ggplot(dfrSrv, aes(x=as.numeric(z),y=value)) +
    geom_line() + geom_hline(yintercept=0,linetype=3,colour="grey") +
    facet_grid(m+p~x) +
    labs(x="size",y="survey catchability") +
    wtsPlots::getStdTheme();
}
plotSrvVals(1,1,1)
plotSrvVals(2,1,1)
plotSrvVals(1,2,1)
plotSrvVals(2,2,1)

#--test in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate selectivity functions----
  info = inputs$lstSel;
  lstSelVals = calcSelectivity(dims,info,params,verbose);
  REPORT(lstSelVals);
  #--calculate survey catchability time series----
  info = inputs$lstSrv;
  lstSrvVals = calcFisheriesRates(dims,info,params,lstSelVals,verbose);
  REPORT(lstSrvVals);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

verbose=FALSE;
params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=map,silent=FALSE);
rep = obj$report();
lstSrvVals = rep$lstSrv;
plotSrvVals(1,1,1)
plotSrvVals(2,1,1)
plotSrvVals(1,2,1)
plotSrvVals(2,2,1)

