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
  source(file.path(dirPrj,"R","readParamInfo_Recruitment_CategoryProportions.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Recruitment_CategoryProportions.R"))
  source(file.path(dirPrj,"R","calcRecruitment_CategoryProportions.R"))
}

inputs = list();

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "data-vertical";
if (type=="data-vertical"){
  ###--recruitment category proportionss with data-vertical----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn     = file.path(dirPrj,"testing/testRecruitment_CategoryProportions/inputSpecs_Recruitment_CategoryProportions.data-vertical.txt");
  res      = readParamInfo_Recruitment_CategoryProportions(conn,TRUE);
  lstRecCRs = extractParamInfo_Recruitment_CategoryProportions(res,dims,FALSE);
  params = list(pRecCRs_FPs=lstRecCRs$params);#--"FP" for "fixed" parameters
} else
if (type=="data-horizontal"){
  ###--recruitment category proportionss with data-horizontal----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn     = file.path(dirPrj,"testing/testRecruitment_CategoryProportions/inputSpecs_Recruitment_CategoryProportions.data-horizontal.txt");
  res      = readParamInfo_Recruitment_CategoryProportions(conn,TRUE);
  lstRecCRs = extractParamInfo_Recruitment_CategoryProportions(res,dims);
  params = list(pRecCRs_FPs=lstRecCRs$params);#--"FP" for "fixed" parameters
} else
if (type=="function"){
  ###--recruitment category proportionss with function----
  dims = setupModelDims(zcs=seq(25.5,84.5,5));
  conn     = file.path(dirPrj,"testing/testRecruitment_CategoryProportions/inputSpecs_Recruitment_CategoryProportions.function.txt");
  res      = readParamInfo_Recruitment_CategoryProportions(conn,FALSE);
  lstRecCRs = extractParamInfo_Recruitment_CategoryProportions(res,dims,FALSE);
  params = list(pRecCRs_MPs=lstRecCRs$MPs$params);
  if (!is.null(lstRecCRs$OPs$params)) params[["pRecCRs_OPs"]]=lstRecCRs$OPs$params;
  if (!is.null(lstRecCRs$DPs$params)) params[["pRecCRs_DPs"]]=lstRecCRs$DPs$params;
  if (!is.null(lstRecCRs$REs$params)) params[["pRecCRs_REs"]]=lstRecCRs$REs$params;
}
inputs          = list();
inputs$dims     = dims;
inputs$lstRecCRs = lstRecCRs;#--add lstRecCRs to inputs

#--test prZatR function----
source(file.path(dirPrj,"R/calcRecruitment_CategoryProportions.R"));
prZatR = calcRecruitment_CategoryProportions(inputs$dims,inputs$lstRecCRs,params,FALSE,loopIC_=TRUE); #--slower
prZatR = calcRecruitment_CategoryProportions(inputs$dims,inputs$lstRecCRs,params,FALSE,loopIC_=FALSE);#--faster
dfrPrZ = tibble::as_tibble(prZatR[1,1,]) |>
         dplyr::mutate(row_id=dplyr::row_number()) |>
         dplyr::mutate(r=dims$dmsC$r[row_id],
                       x=dims$dmsC$x[row_id],
                       m=dims$dmsC$m[row_id],
                       p=dims$dmsC$p[row_id],
                       z=as.numeric(as.character(dims$dmsC$z[row_id])),
                       mp=paste(m,p),
                       rxmp_from=paste(r,x,m,p));
ggplot(dfrPrZ, aes(x=z,y=value)) +
  geom_line() + geom_hline(yintercept=0,linetype=3,colour="white") +
  facet_grid(x~mp) +
  labs(x="size (mm CW)",y="recruitment category proportionss") +
  wtsPlots::getStdTheme()

#--test prZatR in RTMB objective function----
map = list(pRecCRs_MPs=factor(c(1,2,NA,NA),levels=c(1,2,3,4)));
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate recruitment category proportionsss----
  info = inputs$lstRecCRs;
  prZatR = calcRecruitment_CategoryProportions(dims,info,params,verbose,loopIC_=FALSE);
  REPORT(prZatR);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

verbose=FALSE;
params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=map,silent=FALSE);
rep = obj$report();
sum(rep$prZatR[1,1,]);#--should be = 1
dfrPrZ = tibble::as_tibble(rep$prZatR[1,1,]) |>
         dplyr::mutate(row_id=dplyr::row_number()) |>
         dplyr::mutate(r=dims$dmsC$r[row_id],
                       x=dims$dmsC$x[row_id],
                       m=dims$dmsC$m[row_id],
                       p=dims$dmsC$p[row_id],
                       z=as.numeric(as.character(dims$dmsC$z[row_id])),
                       mp=paste(m,p),
                       rxmp_from=paste(r,x,m,p));
ggplot(dfrPrZ, aes(x=z,y=value)) +
  geom_line() + geom_hline(yintercept=0,linetype=3,colour="white") +
  facet_grid(x~mp) +
  labs(x="size (mm CW)",y="recruitment category proportions") +
  wtsPlots::getStdTheme()
