#--test "function" fishery capture rates specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfo_FisheriesCaptureRates.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_FisheriesCaptureRates.R"))
  source(file.path(dirPrj,"R","calcFisheriesCaptureRates.R"))
}

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "function";
if (type=="function"){
  dims = setupModelDims(zcs=seq(25.5,84.5,5));
  params = list();
  map = list();
  ###--fisheries capture rates inputs with function format----
  conn   = file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FisheriesCaptureRates.function.txt");
  res    = readParamInfo_FisheriesCaptureRates(conn,FALSE);
  lstFCRs = extractParamInfo_FisheriesCaptureRates(res,dims,FALSE);
  params[["pFCRs_MPs"]]=lstFCRs$MPs$params;
  map = addList(map,lstFCRs$map);
  if (!is.null(lstFCRs$OPs$params)) params[["pFCRs_OPs"]]=lstFCRs$OPs$params;
  if (!is.null(lstFCRs$DPs$params)) params[["pFCRs_DPs"]]=lstFCRs$DPs$params;
  if (!is.null(lstFCRs$REs$params)) params[["pFCRs_REs"]]=lstFCRs$REs$params;
}
inputs         = list();
inputs$dims    = dims;
inputs$lstFCRs = lstFCRs;#--add lstFCRs to inputs

#--test calcFisheriesCaptureRates function----
source(file.path(dirPrj,"R/calcFisheriesCaptureRates.R"));
lstFCRsVals = calcFisheriesCaptureRates(inputs$dims,inputs$lstFCRs,params,verbose=TRUE);

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
dfrFCRs = getDFR(dims,lstFCRsVals);

#--plot results----
plotFCRsVals<-function(if_=1,iy_=1,is_=1){
  dfrFCRs = dims$dmsC |>
            dplyr::mutate(value=(lstFCRsVals[[if_]])[iy_,is_,],
                          z=as.numeric(as.character(z)));
  ggplot(dfrFCRs, aes(x=as.numeric(z),y=value)) +
    geom_line() + geom_hline(yintercept=0,linetype=3,colour="grey") +
    facet_grid(m+p~x) +
    labs(x="size",y="fisehries capture rates") +
    wtsPlots::getStdTheme();
}
plotFCRsVals(1,1,1)
plotFCRsVals(2,1,1)
plotFCRsVals(1,2,1)
plotFCRsVals(2,2,1)

#--test in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate fisheries capture rates time series----
  info = inputs$lstFCRs;
  lstFCRsVals = calcFisheriesCaptureRates(dims,info,params,verbose=verbose);
  REPORT(lstFCRsVals);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

verbose=FALSE;
params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=map,silent=FALSE);
rep = obj$report();
lstFCRsVals = rep$lstFCRs;
dfrFCRs = getDFR(dims,lstFCRsVals);
plotFCRsVals(1,1,1)
plotFCRsVals(2,1,1)
plotFCRsVals(1,2,1)
plotFCRsVals(2,2,1)

