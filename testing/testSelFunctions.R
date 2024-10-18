#--script to test functions meant for RTMB
require(ggplot2);
require(RTMB);
dirPrj = rstudioapi::getActiveProject();
source(file.path(dirPrj,"R/MiscADMBfunctions.R"),local=TRUE);
source(file.path(dirPrj,"R/SelectivityFunctions.R"),local=TRUE);
source(file.path(dirPrj,"R/getDFR_sdreport.R"),local=TRUE);
source(file.path(dirPrj,"testing/compareSelFun.R"),local=TRUE);

#--define size bins
z = seq(25,180,by=5);

#--sel function: dbllogistic5095(z,p,refZ)----
params = c( 50,  #--ascending limb size at 50% selected
           100,  #--ascending limb size at 95% selected
           125,  #--descending limb size at 95% selected
           150); #--descending limb size at 50% selected
refZ  = 0;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(dbllogistic5095,z,params,refZ,title="dbllogistic5095(z,p,refZ)");

#--sel function: dbllogistic(z,p,refZ)----
params = c( 50,  #--ascending limb size at which selectivity = 0.5 (logit-scale mean)
           0.5,  #--ascending limb slope at 50% selected
           125,  #--descending limb size at which selectivity = 0.5 (logit-scale mean)
           0.5); #--descending limb size at 50% selected
refZ  = 0;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(dbllogistic,z,params,refZ,title="dbllogistic(z,p,refZ)");

#--sel function: asclogistic50D95(z,p,refZ)----
params = c( 50,  #--size at which selectivity = 0.5 (logit-scale mean)
            75); #--difference between sizes at 95\% and 50\%-selected
refZ  = 0;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(asclogistic50D95,z,params,refZ,title="asclogistic50D95(z,p,refZ)");

#--sel function: asclogistic5095(z,p,refZ)----
params = c( 50,  #--size at which selectivity = 0.50 (logit-scale mean)
           125); #--size at which selectivity = 0.95
refZ  = 0;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(asclogistic5095,z,params,refZ,title="asclogistic5095(z,p,refZ)");

#--sel function: asclogistic(z,p,refZ)----
params = c( 50,  #--size at which selectivity = 0.50 (logit-scale mean)
           0.5); #--slope at selectivity = 0.95
refZ  = 0;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(asclogistic,z,params,refZ,title="asclogistic(z,p,refZ)");

#--sel function: const_sel(z,p,refZ)----
#----function should not really be differentiable
params = c(1); #--constant value
refZ  = 0;     #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(const_sel,z,params,refZ,title="const_sel(z,p,refZ)");

#--sel function: ascnormal(z,p,refZ)----
params = c(25,   #--width of ascending limb
           100);  #--size at which ascending limb reaches 1
refZ  = 50;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(ascnormal,z,params,refZ,title="ascnormal(z,p,refZ)");

#--sel function: ascnormal2(z,p,refZ)----
params = c(0.5,   #--selectivity at size = refZ
           100);  #--size at which ascending limb reaches 1
refZ  = 50;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(ascnormal2,z,params,refZ,title="ascnormal2(z,p,refZ)");

#--sel function: ascnormal2a(z,p,refS)----
params = c(50,   #--size at sel = refS
           100); #--size at which ascending limb reaches 1
refS  = 0.5;      #--reference *selectivity*
#--set `f` to function in global environment to be tested
compareSelFun(ascnormal2a,z,params,refS,title="ascnormal2a(z,p,refS)");

#--sel function: ascnormal2b(z,p,refS)----
params = c(150,  #--size at which ascending limb reaches 1
            75); #--delta from size at 1 to size at which selectivity=refS
refS  = 0.5;     #--reference *selectivity*
#--set `f` to function in global environment to be tested
compareSelFun(ascnormal2b,z,params,refS,title="ascnormal2b(z,p,refS)");

#--sel function: ascnormal3(z,p,refS)----
params = c(50,   #--delta from max possible size (refZ[1]) at which ascending limb could reach 1
           0.5); #--selectivity at size=refZ[2]
refZ  = c(185,   #--max possible size at which the curve could reach 1 (e.g., max(z))
           75);  #--reference *size* at which curve reaches params[2]
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,s=ascnormal3(z,params,refZ)),aes(x=z,y=s)) + geom_line();
compareSelFun(ascnormal3,z,params,refZ,title="ascnormal3(z,p,refZ)");

#--sel function: dblnormal4(z,p,refZ)----
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            40,   #--offset to size at which descending limb departs from 1
            35); #--width of descending limb
refZ  = 0;       #--reference *size*
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,s=dblnormal4(z,params,refZ)),aes(x=z,y=s)) + geom_line();
compareSelFun(dblnormal4,z,params,refZ,title="dblnormal4(z,p,refZ)");

#--sel function: dblnormal4a(z,p,refZ)----
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            0.5, #--scaled increment to params[1] at which descending limb departs from 1
            35); #--width of descending limb
refZ  = 185;       #--reference *size*: max possible size (e.g., max(z))
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,s=dblnormal4a(z,params,refZ)),aes(x=z,y=s)) + geom_line();
compareSelFun(dblnormal4a,z,params,refZ,title="dblnormal4a(z,p,refZ)");

#--sel function: dblnormal6(z,p,refZ)----
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
           125,  #--size at which descending limb departs from 1
            35,  #--width of descending limb
           0.1,  #--floor of ascending limb
           0.1); #--floor of descending limb
refZ  = 185;       #--reference *size*: max possible size (e.g., max(z))
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,s=dblnormal6(z,params,refZ)),aes(x=z,y=s)) + geom_line();
compareSelFun(dblnormal6,z,params,refZ,title="dblnormal6(z,p,refZ)");





