#--test recruitment size distribution specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfo_Recruitment_SizeDistribution.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Recruitment_SizeDistribution.R"))
  source(file.path(dirPrj,"R","calcRecruitment_SizeDistribution.R"))
}

inputs = list();

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "data-vertical";
if (type=="data-vertical"){
  ###--recruitment size distribution with data-vertical----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn     = file.path(dirPrj,"testing/testRecruitment_SizeDistribution/inputSpecs_Recruitment_SizeDistribution.data-vertical.txt");
  res      = readParamInfo_Recruitment_SizeDistribution(conn,TRUE);
  lstRecZ = extractParamInfo_Recruitment_SizeDistribution(res,dims,FALSE);
  params = list(pRecZ_FPs=lstRecZ$params);#--"FP" for "fixed" parameters
} else
if (type=="data-horizontal"){
  ###--recruitment size distribution with data-horizontal----
  dims = setupModelDims(zcs=seq(24.5,184.5,5));
  conn     = file.path(dirPrj,"testing/testRecruitment_SizeDistribution/inputSpecs_Recruitment_SizeDistribution.data-horizontal.txt");
  res      = readParamInfo_Recruitment_SizeDistribution(conn,TRUE);
  lstRecZ = extractParamInfo_Recruitment_SizeDistribution(res,dims);
  params = list(pRecZ_FPs=lstRecZ$params);#--"FP" for "fixed" parameters
} else
if (type=="function"){
  ###--recruitment size distribution with function----
  dims = setupModelDims(zcs=seq(25.5,84.5,5));
  conn     = file.path(dirPrj,"testing/testRecruitment_SizeDistribution/inputSpecs_Recruitment_SizeDistribution.function.txt");
  res      = readParamInfo_Recruitment_SizeDistribution(conn,FALSE);
  lstRecZ = extractParamInfo_Recruitment_SizeDistribution(res,dims,FALSE);
  params = list(pRecZ_MPs=lstRecZ$MPs$params);
  if (!is.null(lstRecZ$OPs$params)) params[["pRecZ_OPs"]]=lstRecZ$OPs$params;
  if (!is.null(lstRecZ$DPs$params)) params[["pRecZ_DPs"]]=lstRecZ$DPs$params;
  if (!is.null(lstRecZ$REs$params)) params[["pRecZ_REs"]]=lstRecZ$REs$params;
}
inputs          = list();
inputs$dims     = dims;
inputs$lstRecZ = lstRecZ;#--add lstRecZ to inputs

#--test prZatR function----
source(file.path(dirPrj,"R/calcRecruitment_SizeDistribution.R"));
prZatR = calcRecruitment_SizeDistribution(inputs$dims,inputs$lstRecZ,params,FALSE,loopIC_=TRUE); #--slower
prZatR = calcRecruitment_SizeDistribution(inputs$dims,inputs$lstRecZ,params,FALSE,loopIC_=FALSE);#--faster
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
  labs(x="size (mm CW)",y="recruitment size distribution") +
  wtsPlots::getStdTheme()

#--test prZatR in RTMB objective function----
map = list(pRecZ_MPs=factor(c(1,2,NA,NA),levels=c(1,2,3,4)));
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate recruitment size distributions----
  info = inputs$lstRecZ;
  prZatR = calcRecruitment_SizeDistribution(dims,info,params,verbose,loopIC_=FALSE);
  REPORT(prZatR);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

verbose=FALSE;
params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=map,silent=FALSE);
rep = obj$report();
sum(rep$prZatR[1,1,]);#--should be = number of population categories recruitment occurs in
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
  labs(x="size (mm CW)",y="recruitment size distribution") +
  wtsPlots::getStdTheme()
