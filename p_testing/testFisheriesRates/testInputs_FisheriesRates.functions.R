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
  source(file.path(dirPrj,"R","Functions_Selectivity.R"))
  source(file.path(dirPrj,"R","readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","readParamInfo_Selectivity.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Selectivity.R"))
  source(file.path(dirPrj,"R","calcSelectivity.R"))
  source(file.path(dirPrj,"R","readParamInfo_FisheriesCaptureRates.R"))
  source(file.path(dirPrj,"R","extractParamInfo_FisheriesCaptureRates.R"))
  source(file.path(dirPrj,"R","calcFisheriesCaptureRates.R"))
  source(file.path(dirPrj,"R","readParamInfo_FisheriesHandlingMortalityRates.R"))
  source(file.path(dirPrj,"R","extractParamInfo_FisheriesHandlingMortalityRates.R"))
  source(file.path(dirPrj,"R","calcFisheriesHandlingMortalityRates.R"))
  source(file.path(dirPrj,"R","readInfoType1.R"))
  source(file.path(dirPrj,"R","extractInfoType1.R"))
  source(file.path(dirPrj,"R","readInfo_FisheriesRates.R"))
  source(file.path(dirPrj,"R","extractInfo_FisheriesRates.R"))
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
conn    = file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FshSelectivity.function.txt");
res     = readParamInfo_Selectivity(conn,FALSE);
lstSels = extractParamInfo_Selectivity(res,dims,FALSE);
params[["pSel_MPs"]]=lstSels$MPs$params;
map = addList(map,lstSels$map);
if (!is.null(lstSels$OPs$params)) params[["pSel_OPs"]]=lstSels$OPs$params;
if (!is.null(lstSels$DPs$params)) params[["pSel_DPs"]]=lstSels$DPs$params;
if (!is.null(lstSels$REs$params)) params[["pSel_REs"]]=lstSels$REs$params;
inputs$lstSels = lstSels;#--add lstSel to inputs

##--get fishery capture rates information
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
conn    = file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FisheriesHandlingMortalityRates.function.txt");
res     = readParamInfo_FisheriesHandlingMortalityRates(conn,FALSE);
lstHMRs = extractParamInfo_FisheriesHandlingMortalityRates(res,dims,FALSE);
params[["pHMRs_MPs"]]=lstHMRs$MPs$params;
map = addList(map,lstHMRs$map);
if (!is.null(lstHMRs$OPs$params)) params[["pHMRs_OPs"]]=lstHMRs$OPs$params;
if (!is.null(lstHMRs$DPs$params)) params[["pHMRs_DPs"]]=lstHMRs$DPs$params;
if (!is.null(lstHMRs$REs$params)) params[["pHMRs_REs"]]=lstHMRs$REs$params;
inputs$lstHMRs = lstHMRs;#--add lstHMRs to inputs

###--get fishery rates info----
conn   = file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FisheriesRates.txt");
res    = readInfo_FisheriesRates(conn,FALSE);
lstFsh = extractInfo_FisheriesRates(res,dims,FALSE);
inputs$lstFsh = lstFsh;#--add lstFsh to inputs

#--test calcFisheriesRates function----
lstArrSels = calcSelectivity(dims,lstSels,params,FALSE);
lstArrFCRs = calcFisheriesCaptureRates(dims,lstFCRs,params,FALSE);
lstArrHMRs = calcFisheriesHandlingMortalityRates(dims,lstHMRs,params,FALSE);
source(file.path(dirPrj,"R/calcFisheriesRates.R"));
lstFshVals = calcFisheriesRates(dims,inputs$lstFsh,lstArrFCRs,lstArrHMRs,lstArrSels,TRUE);

#--convert list/arrays to dataframe----
getDFR3<-function(dims,lstArr){
  lstDfrs = list();
  nms = names(lstArr);
  for (il_ in 1:length(nms)){
    #--testing: il_=1;
    nm = nms[il_];
    arrYSC = lstArr[[nm]];
    dfrp = dims$dmsYSC;
    dfrp$value = as.vector(aperm(arrYSC,perm=c(3,2,1)));
    lstDfrs[[nm]] = dfrp |> dplyr::mutate(flt=nm);
    rm(dfrp);
  }
  return(dplyr::bind_rows(lstDfrs));
}
dfrFCRs = getDFR3(dims,lstArrFCRs) ; View(dfrFCRs);
dfrHMRs = getDFR3(dims,lstArrHMRs) ; View(dfrHMRs);
dfrSels = getDFR3(dims,lstArrSels) ; View(dfrSels);

getDFR4<-function(dims,lstArr,
                  types=c("FCR","RMR","DMR","TFMR","cap_rt","hm_rt","cap_sel","ret_sel")){
  lstDfrs = list();
  nms = names(lstArr);
  for (il_ in 1:length(nms)){
    #--testing: il_=1;
    nm = nms[il_];
    arrYSC = lstArr[[nm]];
    dfrp = tibble::tibble(type=types) |>
             dplyr::cross_join(dims$dmsYSC);
    dfrp$value = as.vector(aperm(arrYSC,perm=c(4,3,2,1)));
    lstDfrs[[nm]] = dfrp |> dplyr::mutate(flt=nm);
    rm(dfrp);
  }
  return(dplyr::bind_rows(lstDfrs));
}
dfrFsh = getDFR4(dims,lstFshVals) |> dplyr::select(!sparse_idx) |>
           tidyr::pivot_wider(names_from="type",values_from="value");
View(dfrFsh);

#--test in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate selectivity functions, fully-selected capture rates, handling mortality rates----
  lstArrSels = calcSelectivity(dims,inputs$lstSels,params,FALSE);
  lstArrFCRs = calcFisheriesCaptureRates(dims,inputs$lstFCRs,params,FALSE);
  lstArrHMRs = calcFisheriesHandlingMortalityRates(dims,inputs$lstHMRs,params,FALSE);
  REPORT(lstArrSels);
  REPORT(lstArrFCRs);
  REPORT(lstArrHMRs);
  #--calculate fisheries mortality rates----
  lstSrvVals = calcFisheriesRates(dims,inputs$lstFsh,lstArrFCRs,lstArrHMRs,lstArrSels,TRUE);
  REPORT(lstFshVals);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

verbose=FALSE;
params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=map,silent=FALSE);
rep = obj$report();
dfrFshVals = getDFR4(dims,rep$lstFshVals) |> dplyr::select(!sparse_idx) |>
           tidyr::pivot_wider(names_from="type",values_from="value");
View(dfrFshVals);


