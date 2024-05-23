#--RTMB version of Andre's Ex1a
#
require(RTMB);
require(ggplot2);
require(dplyr);
#compiler::enableJIT(0);

dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
FileName <- file.path(dirThs,"ex1a.dat");

# Read in the control parameters
Nyear <- scan(FileName,skip=1,n=1,quiet=T)
Nclass <- scan(FileName,skip=3,n=1,quiet=T)

# Read biological parameters
Length <- scan(FileName,skip=5,n=Nclass,quiet=T);#--size class midpoints
WL <- scan(FileName,skip=7,n=Nclass,quiet=T);    #--weight-at-size
X <- matrix(0,nrow=Nclass,ncol=Nclass)           #--size class transition matrix
for (Iclass in 1:Nclass) X[Iclass,] <- scan(FileName,skip=9+(Iclass-1),n=Nclass,quiet=T)
skip <- 9+Nclass;

# Read initial values for fisehry and survey logistic selectivities
S50 <- scan(FileName,skip=skip+1,n=1,quiet=T);#--fishery selectivity
S95 <- scan(FileName,skip=skip+2,n=1,quiet=T);#--fishery selectivity
SS50 <- scan(FileName,skip=skip+4,n=1,quiet=T);#--survey selectivity
SS95 <- scan(FileName,skip=skip+5,n=1,quiet=T);#--survey selectivity

# Natural mortality
M <- scan(FileName,skip=skip+7,n=1,quiet=T)
skip <- skip + 8

# Catches
CWObs <- rep(0,Nyear)
for (Iyear in 1:Nyear) CWObs[Iyear] <- scan(FileName,skip=skip+Iyear,n=1,quiet=T)
skip <- skip +Nyear+1

# Catch-at-length
CALObs <- matrix(0,nrow=Nyear,ncol=Nclass)
for (Iyear in 1:Nyear) CALObs[Iyear,] <- scan(FileName,skip=skip+Iyear,n=Nclass,quiet=T)
for (Iyear in 1:Nyear) CALObs[Iyear,] <- CALObs[Iyear,]/sum(CALObs[Iyear,])
Neff <- scan(FileName,skip=skip+Nyear+2,n=1,quiet=T)
skip <- skip+Nyear+2+1

# Biomass index
BioIndex <- rep(0,Nyear)
for (Iyear in 1:Nyear) BioIndex[Iyear] <- scan(FileName,skip=skip+Iyear,n=1,quiet=T)
BioSig <- scan(FileName,skip=skip+Nyear+2,n=1,quiet=T)

# Set up fishery selectivity
S = 1.0/(1+exp(-log(19.0)*(Length-S50)/(S95-S50)));
# Set up survey selectivity
SurveyS = 1.0/(1+exp(-log(19.0)*(Length-SS50)/(SS95-SS50)));

source(file.path(dirThs,"rtmb_objfun.R"),chdir=TRUE);
source(file.path(dirThs,"rtmb_getDFR_sdreport.R"),chdir=TRUE);
source(file.path(dirThs,"rtmb_plotResults.R"),chdir=TRUE);

#--parameter values
ParVec <- c(0.0967594569308,
            0.870434495402,0.200399053425,1.05381737431,0.318731122710,-0.918363665392,
            -2.16879783519,-1.80398463200,-1.51453776273,-1.29629862470,-1.11179305214,-0.936759780099,-0.780924556922,-0.652694662549,-0.571948582890,-0.374205217483,-0.370671729651,-0.327256603866,-0.261995846926,-0.231684378053,-0.202194242986,-0.169478399063,-0.218943322854,-0.288907219146,-0.325585728593,-0.296255891488,-0.316494143474,-0.370454824963,-0.429089329265,-0.525589367564,-0.648138027610,
            0.832106344761,0.194604826143,-1.04054310341,0.746052244377,-0.996844601768,-0.518249733505,-0.769018641558,0.356647903296,-0.436211851283,0.186658837889,-0.132747421365,0.720642981419,0.297371523264,-0.671838396088,0.390538452094,0.344467053149,0.356127773876,-0.339423918902,0.561520251858,0.515069864196,-0.507304928860,-0.337304176785,-0.167803906345,0.415482624421,0.00000000000)
nPV = length(ParVec);
#--add selectivity parameter values (so they can be estimated at some point)
ParVec = c(ParVec,log(S50),log(S95),log(SS50),log(SS95),log(S),log(SurveyS));
#--create some indices into ParVec for selectivity parameters
ipS50  = nPV+1; nPV = nPV+1;  #--index of (ln-scale) parameter for S50 (fishery)
ipS95  = nPV+1; nPV = nPV+1;  #--index of (ln-scale) parameter for S95 (fishery)
ipSS50 = nPV+1; nPV = nPV+1;  #--index of (ln-scale) parameter for SS50 (survey)
ipSS95 = nPV+1; nPV = nPV+1;  #--index of (ln-scale) parameter for SS95 (survey)
ipSps   = (nPV+1):(nPV+Nclass); nPV = nPV+Nclass;#--indices for nonparametric fishery selectivity
ipSSps  = (nPV+1):(nPV+Nclass); nPV = nPV+Nclass;#--indices for nonparametric survey selectivity

#-------------------------------------------------------------------------------
#--CASE 1: create parameters list for TMB object, don't estimate selectivity functions
params = list(pars=ParVec); #--define input parameters list
#----create list to turn on/off estimation of individual parameters
mpars1=1:nPV;
#----turn of estimation of all selectivity parameters
mpars1[ipS50]  = NA;#--turn off estimation
mpars1[ipS95]  = NA;#--turn off estimation
mpars1[ipSS50] = NA;#--turn off estimation
mpars1[ipSS95] = NA;#--turn off estimation
mpars1[ipSps]  = NA;#--turn off estimation
mpars1[ipSSps] = NA;#--turn off estimation
map1 = list(pars=factor(mpars1));#--list to group parameters or turn on/off estimation

#--create TMB model object
#----fix selectivity parameters, use parametric representations
useParFishSel = TRUE;
useParSurvSel = TRUE;
objfun(params);#--check objective function using standard R calculations

#--create RTMB object
if (exists("objfun")) {
  #--need to "re-source" the objfun functions if objfun(params) has been called
  rm(objfun,calcLogisticSel,calcFandZ,calcN);
  source(file.path(dirThs,"rtmb_objfun.R"),chdir=TRUE);
}
obj1 <- RTMB::MakeADFun(objfun,
                       params,
                       map=map1,
                       silent=TRUE);
obj1$fn();#--check objective function value at input parameters
obj1$report(obj1$par);#--check reported values (parts of likelihood)

#--estimate the parameters
opt1 = nlminb(obj1$par,obj1$fn,obj1$gr,
             control=list(eval.max=10000,iter.max=1000,trace=100));
obj1$report()

#--get the std devs
sdr1 = RTMB::sdreport(obj1);
dfr1a = getDFR_sdreport(obj1,report=FALSE,model="model 1");
dfr1b = getDFR_sdreport(obj1,report=TRUE, model="model 1");
dfr1b |> dplyr::filter(name %in% c("LikeCatch","LikeCAL","LikeBio","Penal","obj_fun"))
plotResults(dfr1b);

#-------------------------------------------------------------------------------
#--CASE 2: create TMB model object for "fuller" model, estimating parametric selectivity parameters
useParFishSel = TRUE;
useParSurvSel = TRUE;
ParVec2 = 0*ParVec;
ParVec2[ipS50]  = ParVec[ipSS50];#--initial value for S50: mix initial value up a bit
ParVec2[ipS95]  = ParVec[ipSS95];#--initial value for S95: mix initial value up a bit
ParVec2[ipSS50] = ParVec[ipS50];#--initial value for SS50: mix initial value up a bit
ParVec2[ipSS95] = ParVec[ipS95];#--initial value for SS95: mix initial value up a bit
ParVec2[ipSps]  = log(S);       #--initial values for nonparametric fishery selectivity (not estimated)
ParVec2[ipSSps] = log(SurveyS); #--initial values for nonparametric survey selectivity (not estimated)
#----create list to estimate selectivity parameters
mpars2=1:nPV;      #--turn on all parameters
mpars2[ipSps]  = NA;#--turn off estimation of nonparametric fishery selectivity
mpars2[ipSSps] = NA;#--turn off estimation of nonparametric survey selectivity
map2 = list(pars=factor(mpars2));#--"map" to group parameters or turn on/off estimation

#--create TMB object
obj2 <- RTMB::MakeADFun(objfun,
                        parameters=list(pars=ParVec2),
                        map=map2,
                        silent=TRUE);

#--estimate the parameters
opt2 = nlminb(obj2$par,obj2$fn,obj2$gr,hessian=obj2$he,
              control=list(eval.max=10000,iter.max=1000,trace=100));

#--get the sdreport
sdr2 = RTMB::sdreport(obj2);
dfr2a = getDFR_sdreport(obj2,report=FALSE,model="model 2");
dfr2b = getDFR_sdreport(obj2,report=TRUE,model="model 2");
#--look at differences
plotResults(dplyr::bind_rows(dfr1b,dfr2b));

#-------------------------------------------------------------------------------
#--CASE 3: create TMB model object for nonparametric fishery selectivity model
#----does not estimate (parametric) survey selectivity
useParFishSel = FALSE;
useParSurvSel = TRUE;
ParVec3 = 0*ParVec;
ParVec3[ipS50]  = ParVec[ipS50];#--initial value for S50 (not estimated)
ParVec3[ipS95]  = ParVec[ipS95];#--initial value for S95 (not estimated)
ParVec3[ipSS50] = ParVec[ipSS50];#--initial value for SS50 (not estimated)
ParVec3[ipSS95] = ParVec[ipSS95];#--initial value for SS95 (not estimated)
ParVec3[ipSps[Nclass]] = 0;#--initial value for ln-scale fishery selectivity in Nclass (not estimated)
ParVec3[ipSSps] = log(SurveyS); #--initial values for nonparametric survey selectivity (not estimated)
#--above uses 0's for initial values for nonparametric fishery selectivity curve
#----create list to turn off estimation of parametric selectivity parameters and
#----nonparametric survey selectivity
mpars3=1:nPV;
mpars3[ipS50]  = NA;#--turn off estimation
mpars3[ipS95]  = NA;#--turn off estimation
mpars3[ipSS50] = NA;#--turn off estimation
mpars3[ipSS95] = NA;#--turn off estimation
mpars3[ipSps[Nclass]] = NA;#--turn off estimation in Nclass
mpars3[ipSSps] = NA;#--turn off estimation
map3 = list(pars=factor(mpars3));#--list to group parameters or turn on/off estimation

#--create TMB object
obj3 <- RTMB::MakeADFun(objfun,
                        parameters=list(pars=ParVec3),
                        map=map3,
                        silent=TRUE);

#--estimate the parameters
opt3 = nlminb(obj3$par,obj3$fn,obj3$gr,hessian=obj3$he,
              control=list(eval.max=10000,iter.max=1000,trace=1));

#--get the sdreport
sdr3 = RTMB::sdreport(obj3);
dfr3a = getDFR_sdreport(obj3,report=FALSE,model="model 3");
dfr3b = getDFR_sdreport(obj3,report=TRUE,model="model 3");
plotResults(dplyr::bind_rows(dfr1b,dfr2b,dfr3b));
#--look at differences
dfrp = dfr3a |> dplyr::mutate(Par3=ParVec) |>
               dplyr::mutate(diff=Par3-Est,
                             dev =(Par3-Est)/Std);

#-------------------------------------------------------------------------------
#--CASE 4: create TMB model object for nonparametric survey selectivity model
#----does not estimate (parametric) fishery selectivity
useParFishSel = TRUE;
useParSurvSel = FALSE;
ParVec4 = 0*ParVec;
ParVec4[ipS50]  = ParVec[ipS50];#--initial value for S50 (not estimated)
ParVec4[ipS95]  = ParVec[ipS95];#--initial value for S95 (not estimated)
ParVec4[ipSS50] = ParVec[ipSS50];#--initial value for SS50 (not estimated)
ParVec4[ipSS95] = ParVec[ipSS95];#--initial value for SS95 (not estimated)
ParVec4[ipSps]  = log(S);   #--initial values for ln-scale nonparametric fishery selectivity (not estimated)
ParVec4[ipSSps[Nclass]] = 0;#--value for ln-scale nonparametric survey selectivity in Nclass (not estimated)
#--above uses 0's for initial values for nonparametric selectivity curves
#----create list to turn off estimation of parametric selectivity parameters and
#----nonparametric fishery selectivity
mpars4=1:nPV;
mpars4[ipS50]  = NA;#--turn off estimation
mpars4[ipS95]  = NA;#--turn off estimation
mpars4[ipSS50] = NA;#--turn off estimation
mpars4[ipSS95] = NA;#--turn off estimation
mpars4[ipSps]  = NA;#--turn off estimation
mpars4[ipSSps[Nclass]] = NA;#--turn off estimation in Nclass
map4 = list(pars=factor(mpars4));#--list to group parameters or turn on/off estimation

#--create TMB object
obj4 <- RTMB::MakeADFun(objfun,
                        parameters=list(pars=ParVec4),
                        map=map4,
                        silent=TRUE);

#--estimate the parameters
opt4 = nlminb(obj4$par,obj4$fn,obj4$gr,hessian=obj4$he,
              control=list(eval.max=10000,iter.max=1000,trace=1));

#--get the sdreport
sdr4 = RTMB::sdreport(obj4);
dfr4a = getDFR_sdreport(obj4,report=FALSE,model="model 4");
dfr4b = getDFR_sdreport(obj4,report=TRUE,model="model 4");
plotResults(dplyr::bind_rows(dfr1b,dfr2b,dfr3b,dfr4b));
#--look at differences
dfrp = dfr4a |> dplyr::mutate(Par4=ParVec) |>
               dplyr::mutate(diff=Par4-Est,
                             dev =(Par4-Est)/Std);

#-------------------------------------------------------------------------------
#--CASE 5: create TMB model object for nonparametric fishery AND survey selectivity model
#----does not estimate parametric fishery or survey selectivity
#----estimates nonparametric fishery or survey selectivity
useParFishSel = FALSE;
useParSurvSel = FALSE;
ParVec5 = 0*ParVec;
ParVec5[ipS50]  = ParVec[ipS50];#--initial value for S50 (not estimated)
ParVec5[ipS95]  = ParVec[ipS95];#--initial value for S95 (not estimated)
ParVec5[ipSS50] = ParVec[ipSS50];#--initial value for SS50 (not estimated)
ParVec5[ipSS95] = ParVec[ipSS95];#--initial value for SS95 (not estimated)
ParVec5[ipSps[Nclass]]  = 0;#--value for ln-scale nonparametric fishery selectivity in Nclass (not estimated)
ParVec5[ipSSps[Nclass]] = 0;#--value for ln-scale nonparametric survey selectivity in Nclass (not estimated)
#--above uses 0's for initial values for nonparametric selectivity curves
#----create list to turn off estimation of parametric selectivity parameters
mpars5=1:nPV;
mpars5[ipS50]  = NA;#--turn off estimation
mpars5[ipS95]  = NA;#--turn off estimation
mpars5[ipSS50] = NA;#--turn off estimation
mpars5[ipSS95] = NA;#--turn off estimation
mpars5[ipSps[Nclass]]  = NA;#--turn off estimation in Nclass
mpars5[ipSSps[Nclass]] = NA;#--turn off estimation in Nclass
map5 = list(pars=factor(mpars5));#--list to group parameters or turn on/off estimation

#--create TMB object
obj5 <- RTMB::MakeADFun(objfun,
                        parameters=list(pars=ParVec5),
                        map=map5,
                        silent=TRUE);

#--estimate the parameters
opt5 = nlminb(obj5$par,obj5$fn,obj5$gr,hessian=obj5$he,
              control=list(eval.max=10000,iter.max=1000,trace=1));

#--get the sdreport
sdr5 = RTMB::sdreport(obj5);
dfr5a = getDFR_sdreport(obj5,report=FALSE,model="model 5");
dfr5b = getDFR_sdreport(obj5,report=TRUE,model="model 5");
plotResults(dplyr::bind_rows(dfr1b,dfr2b,dfr3b,dfr4b,dfr5b));
#--look at differences
dfrp = dfr5a |> dplyr::mutate(Par5=ParVec) |>
               dplyr::mutate(diff=Par5-Est,
                             dev =(Par5-Est)/Std);


