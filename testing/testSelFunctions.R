#--script to test functions meant for RTMB
require(ggplot2);
require(RTMB);
dirPrj = rstudioapi::getActiveProject();
source(file.path(dirPrj,"R/SelectivityFunctions.R"),local=TRUE);
source(file.path(dirPrj,"testing/compareSelFun.R"),local=TRUE);
source(file.path(dirPrj,"R/getDFR_sdreport.R"),local=TRUE);

#--define size bins
z = seq(25,180,by=5);

#--sel function: ascnormal2a(z,p,refS)----
params = c(50,   #--size at sel = refS
           100); #--size at which ascending limb reaches 1
refS  = 0.5;      #--reference *selectivity*
#--set `f` to function in global environment to be tested
compareSelFun(ascnormal2a,z,params,refS,title="ascnormal2a(z,p,refS)");

#--sel function: ascnormal2(z,p,refZ)----
params = c(0.5,   #--selectivity at size = refZ
           100);  #--size at which ascending limb reaches 1
refZ  = 50;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(ascnormal2,z,params,refZ,title="ascnormal2(z,p,refZ)");

#--sel function: ascnormal(z,p,refZ)----
params = c(25,   #--width of ascending limb
           100);  #--size at which ascending limb reaches 1
refZ  = 50;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(ascnormal,z,params,refZ,title="ascnormal(z,p,refZ)");

#--sel function: dbllogistic5095(z,p,refZ)----
params = c( 50,  #--ascending limb size at 50% selected
           100,  #--ascending limb size at 95% selected
           125,  #--descending limb size at 95% selected
           150); #--descending limb size at 50% selected
refZ  = 0;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(dbllogistic5095,z,params,refZ,title="dbllogistic5095(z,p,refZ)");





