#--script to test growth functions
require(ggplot2);
require(RTMB);
dirPrj = rstudioapi::getActiveProject();
source(file.path(dirPrj,"R/MiscADMBfunctions.R"),local=TRUE);
source(file.path(dirPrj,"R/GrowthFunctions.R"),local=TRUE);
source(file.path(dirPrj,"R/getDFR_sdreport.R"),local=TRUE);
source(file.path(dirPrj,"testing/compareSelFun.R"),local=TRUE);

#--define size bins
z = seq(2.5,182.5,by=5);

#--mean post-molt function: calcMnPostMoltZ1(z,p,c)----
params = c(log(2.0),  #--ln-scale post-molt size for pre-molt size = 1
           1); #--ln-scale slope of mean growth curve
consts  = NULL;  #--reference *size*
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,y=calcMnPostMoltZ1(z,params,consts)),aes(x=z,y=y)) +
  geom_line() + geom_abline(slope=1,linetype=3);
compareGrowthFun(calcMnPostMoltZ1,z,params,consts,scale=0.8,title="calcMnPostMoltZ1(z,p,c)");

#--mean post-molt function: calcMnPostMoltZ2(z,p,c)----
params = c( 30,  #--mean post-molt size for pre-molt size = consts[1]
           160); #--mean post-molt size for pre-molt size = consts[2]
consts = c( 25,  #--pre-molt size corresponding to mean post-molt size given by params[1]
           140); #--pre-molt size corresponding to mean post-molt size given by params[2]
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,y=calcMnPostMoltZ2(z,params,consts)),aes(x=z,y=y)) +
  geom_line() + geom_abline(slope=1,linetype=3);
compareGrowthFun(calcMnPostMoltZ2,z,params,consts,scale=0.8,title="calcMnPostMoltZ2(z,p,c)");

#--mean post-molt function: calcMnPostMoltZ3(z,p,c)----
params = c(30,  #--mean post-molt size for pre-molt size = consts[1]
            1); #--ln-scale slope of mean growth curve
consts = 25;    #--pre-molt size corresponding to mean post-molt size given by params[1]
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,y=calcMnPostMoltZ3(z,params,consts)),aes(x=z,y=y)) +
  geom_line() + geom_abline(slope=1,linetype=3);
mdl = compareGrowthFun(calcMnPostMoltZ3,z,params,consts,scale=0.8,title="calcMnPostMoltZ3(z,p,c)");
mdl$simulate()

