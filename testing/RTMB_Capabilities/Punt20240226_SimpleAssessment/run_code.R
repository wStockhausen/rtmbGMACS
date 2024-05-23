#--Example where MakeADFun gives error after objective function is first run in R.
#
require(RTMB);
require(ggplot2);
require(dplyr);
compiler::enableJIT(0); ##<--uncommenting removes error at MakeADFun

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
#--CASE 1: create parameters list for RTMB object, don't estimate selectivity functions
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

#--run the code in "R mode"
objfun(params);#--check objective function using standard R calculations

#--create the RTMB object
obj1 <- RTMB::MakeADFun(objfun,
                       params,
                       map=map1,
                       silent=TRUE);


