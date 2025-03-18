#--test "function" fishery handling mortality rates specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfo_FisheriesHandlingMortalityRates.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_FisheriesHandlingMortalityRates.R"))
  source(file.path(dirPrj,"R","calcFisheriesHandlingMortalityRates.R"))
}

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "function";
if (type=="function"){
  dims = setupModelDims(zcs=seq(25.5,84.5,5));
  params = list();
  map = list();
  ###--fisheries handling mortality rates inputs with function format----
  conn   = file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FisheriesHandlingMortalityRates.function.txt");
  res    = readParamInfo_FisheriesHandlingMortalityRates(conn,FALSE);
  lstHMRs = extractParamInfo_FisheriesHandlingMortalityRates(res,dims,FALSE);
  params[["pHMRs_MPs"]]=lstHMRs$MPs$params;
  map = addList(map,lstHMRs$map);
  if (!is.null(lstHMRs$OPs$params)) params[["pHMRs_OPs"]]=lstHMRs$OPs$params;
  if (!is.null(lstHMRs$DPs$params)) params[["pHMRs_DPs"]]=lstHMRs$DPs$params;
  if (!is.null(lstHMRs$REs$params)) params[["pHMRs_REs"]]=lstHMRs$REs$params;
}
inputs         = list();
inputs$dims    = dims;
inputs$lstHMRs = lstHMRs;#--add lstHMRs to inputs

#--test calcFisheriesHandlingMortalityRates function----
source(file.path(dirPrj,"R/calcFisheriesHandlingMortalityRates.R"));
lstHMRsVals = calcFisheriesHandlingMortalityRates(inputs$dims,inputs$lstHMRs,params,verbose=TRUE);

#--convert list/arrays to dataframe----
getDFR<-function(dims,lstArr){
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
dfrHMRs = getDFR(dims,lstHMRsVals);

#--plot results----
plotHMRsVals<-function(if_=1,iy_=1,is_=1){
  dfrHMRs = dims$dmsC |>
            dplyr::mutate(value=(lstHMRsVals[[if_]])[iy_,is_,],
                          z=as.numeric(as.character(z)));
  ggplot(dfrHMRs, aes(x=as.numeric(z),y=value)) +
    geom_line() + geom_hline(yintercept=0,linetype=3,colour="grey") +
    facet_grid(m+p~x) +
    labs(x="size",y="fisheries handling mortality rates") +
    wtsPlots::getStdTheme();
}
plotHMRsVals(1,1,1)
plotHMRsVals(2,1,1)
plotHMRsVals(1,2,1)
plotHMRsVals(2,2,1)

#--test in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate fisheries capture rates time series----
  info = inputs$lstHMRs;
  lstHMRsVals = calcFisheriesHandlingMortalityRates(dims,info,params,verbose=verbose);
  REPORT(lstHMRsVals);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

verbose=FALSE;
params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=map,silent=FALSE);
rep = obj$report();
lstHMRsVals = rep$lstHMRs;
dfrHMRs = getDFR(dims,lstHMRsVals);
plotHMRsVals(1,1,1)
plotHMRsVals(2,1,1)
plotHMRsVals(1,2,1)
plotHMRsVals(2,2,1)

