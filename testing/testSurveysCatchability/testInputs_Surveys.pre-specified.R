#--test "pre-specified" surveys catchability specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfo_SurveysCatchability.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_SurveysCatchability.R"))
  source(file.path(dirPrj,"R","Functions_Catchability.R"))
  source(file.path(dirPrj,"R","calcSurveysCatchability.R"))
}

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "pre-specified-horizontal";
if (type=="pre-specified-vertical"){
  ###--survey catchability with pre-specified-vertical format----
  dims = setupModelDims(zcs=seq(25.5,80.5,5));
  conn   = file.path(dirPrj,"testing/testSurveysCatchability/inputSpecs_SurveysCatchability.pre-specified-vertical.txt");
  res    = readParamInfo_SurveysCatchability(conn,TRUE);
  lstSrv = extractParamInfo_SurveysCatchability(res,dims,FALSE);
  params = list(pSrvQ_FPs=lstSrv$params);#--"FP" for "fixed" parameters
  map = lstSrv$map;
  lstSel = NULL;
} else
if (type=="pre-specified-horizontal"){
  ###--survey catchability inputs with pre-specified-horizontal format----
  dims   = setupModelDims(zcs=seq(25.5,80.5,5));
  conn   = file.path(dirPrj,"testing/testSurveysCatchability/inputSpecs_SurveysCatchability.pre-specified-horizontal.txt");
  res    = readParamInfo_SurveysCatchability(conn,FALSE);
  lstSrv = extractParamInfo_SurveysCatchability(res,dims,FALSE);
  params = list(pSrvQ_FPs=lstSrv$params);#--"FP" for "fixed" parameters
  map = lstSrv$map;
  lstSel = NULL;
}
inputs        = list();
inputs$dims   = dims;
inputs$lstSel = lstSel;#--add lstSel to inputs (but NULL if type=="pre-specified")
inputs$lstSrv = lstSrv;#--add lstSrv to inputs

#--test calcSurveys function----
lstSelVals=NULL;#--don't need selectivities object
source(file.path(dirPrj,"R/calcSurveysCatchability.R"));
lstSrvVals = calcSurveysCatchability(inputs$dims,inputs$lstSrv,params,lstSelVals,TRUE);

#--plot results----
plotSrvVals<-function(if_=1,iy_=1,is_=1){
  dfrSrv = dims$dmsC |>
            dplyr::mutate(value=(lstSrvVals[[if_]])[iy_,is_,],
                          z=as.numeric(as.character(z)));
  ggplot(dfrSrv, aes(x=as.numeric(z),y=value)) +
    geom_hline(yintercept=0,linetype=3,colour="grey") +
    geom_line() +
    facet_grid(m+p~x) +
    labs(x="size",y="survey catchability") +
    wtsPlots::getStdTheme();
}
plotSrvVals(1,1,1)
plotSrvVals(2,1,1)
plotSrvVals(1,4,1)
plotSrvVals(2,4,1)

#--test calcSurveysCatchability in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate fisheries time series----
  info = inputs$lstSrv;
  lstSrv = calcSurveysCatchability(dims,info,params,input$lstSel,verbose);
  REPORT(lstSrv);

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
plotSrvVals(1,4,1)
plotSrvVals(2,4,1)
