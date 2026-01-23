#--script to test selectivity functions meant for RTMB
require(ggplot2);
require(RTMB);
dirPrj = rstudioapi::getActiveProject();
source(file.path(dirPrj,"R/MiscFunctions.R"),local=TRUE);
source(file.path(dirPrj,"R/SelectivityFunctions.R"),local=TRUE);
source(file.path(dirPrj,"R/getDFR_sdreport.R"),local=TRUE);
source(file.path(dirPrj,"R/compareSelFun.R"),local=TRUE);

#--define size bins
z = seq(25,180,by=5);

#--sel function: const_sel(z,p,refZ)----
#----function should not really be differentiable(?)
params = c(1); #--constant value
refZ  = 100;     #--reference *size* (not used!)
#--set `f` to function in global environment to be tested
compareSelFun(const_sel,z,params,refZ,title="const_sel(z,p,refZ)");

#--sel function: asclogistic(z,p,refZ)----
params = c(100.00,  #--size at which selectivity = 0.50 (logit-scale mean)
             0.05); #--width (1/slope) at z50
refZ  = 125;        #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(asclogistic,z,params,refZ,title="asclogistic(z,p,refZ)",verbose=TRUE);

#--sel function: asclogistic1(z,p,refZ)----
params = c(100,   #--size at which selectivity = 0.50 (logit-scale mean)
            20); #--width (1/slope) at z50
refZ  = 125;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(asclogistic1,z,params,refZ,title="asclogistic1(z,p,refZ)");

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

#--sel function: ascnormal1(z,p,refZ)----
params = c(25,   #--width of ascending limb
           100);  #--size at which ascending limb reaches 1
refZ  = 50;       #--reference *size*
#--set `f` to function in global environment to be tested
compareSelFun(ascnormal1,z,params,refZ,title="ascnormal(z,p,refZ)");

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

#--sel function: stackedLogistic1(z,p,refZ)----
params = c(0.2,  #--omega: -weighting factor on the first curve
            75,  #--mnZ1: size at inflection point for 1st logistic curve
            10,  #--sdZ1: sd for 1st logistic curve
           145,  #--mnZ2: size at inflection point for 2nd logistic curve
            10); #--sdZ2: sd for 2nd logistic curve
refZ  = 185;     #--not used
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,s=stackedLogistic1(z,params,refZ)),aes(x=z,y=s)) + geom_line();
compareSelFun(stackedLogistic1,z,params,refZ,title="stackedLogistic1(z,p,refZ)");

#--sel function: selSpline(z,p,knots)----
#--set up a spline curve
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            0.5, #--scaled increment to params[1] at which descending limb departs from 1
            35); #--width of descending limb
refZ  = 185;       #--reference *size*: max possible size (e.g., max(z))
yz = dblnormal4(z,params,refZ);
yz = ifelse(yz<0.0001,0.0001,ifelse(yz>0.9999,0.999,yz));
nz = length(z);
nk = 6; ik = seq(from=1,to=nz,length.out=nk);
xk = z[ik];
lgtyk = log(yz[ik]/(1-yz[ik]));
params = lgtyk;
knots = xk;
sf = selSpline(z,params,knots);
refZ  = 185;     #--not used
#--set `f` to function in global environment to be tested
dfr1 = tibble::tibble(z=z,s=selSpline(z,params,knots));
dfr2 = tibble::tibble(z=z,s=yz);
dfr3 = tibble::tibble(z=xk,s=exp(lgtyk)/(1+exp(lgtyk)));
ggplot(dfr1,aes(x=z,y=s)) + geom_line() + geom_point() +
  geom_point(data=dfr2,colour="blue") +
  geom_point(data=dfr3,colour="red",fill=NA,shape=23,size=6)
compareSelFun(selSpline,z,params,knots,title="selSpline(z,params,knots)");

#--sel function: selSplineClmpd(z,p,knots)----
#--set up a spline curve "clamped" on both ends
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            0.5, #--scaled increment to params[1] at which descending limb departs from 1
            35); #--width of descending limb
yz = dblnormal4(z,params,refZ);
yz = ifelse(yz<0.0001,0.0001,ifelse(yz>0.9999,0.999,yz));
nz = length(z);
nk = 6; ik = seq(from=3,to=nz-3,length.out=nk);
xk = z[ik];
lgtyk = log(yz[ik]/(1-yz[ik]));
params = lgtyk;
knots = xk;
sf = selSplineClmpd(z,params,knots);
refZ  = 185;     #--not used
#--set `f` to function in global environment to be tested
dfr1 = tibble::tibble(z=z,s=selSpline(z,params,knots));
dfr2 = tibble::tibble(z=z,s=yz);
dfr3 = tibble::tibble(z=xk,s=exp(lgtyk)/(1+exp(lgtyk)));
ggplot(dfr1,aes(x=z,y=s)) + geom_line() + geom_point() +
  geom_point(data=dfr2,colour="blue") +
  geom_point(data=dfr3,colour="red",fill=NA,shape=23,size=6)
compareSelFun(selSplineClmpd,z,params,knots,title="selSplineClmpd(z,params,knots)");

#--sel function: selSplineClmpdLeft(z,p,knots)----
#--set up a spline curve "clamped" on the left end
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            0.5, #--scaled increment to params[1] at which descending limb departs from 1
            35); #--width of descending limb
yz = dblnormal4(z,params,refZ);
yz = ifelse(yz<0.0001,0.0001,ifelse(yz>0.9999,0.999,yz));
nz = length(z);
nk = 6; ik = seq(from=3,to=nz-3,length.out=nk);
xk = z[ik];
lgtyk = log(yz[ik]/(1-yz[ik]));
params = lgtyk;
knots = xk;
sf = selSplineClmpdLeft(z,params,knots);
refZ  = 185;     #--not used
#--set `f` to function in global environment to be tested
dfr1 = tibble::tibble(z=z,s=selSpline(z,params,knots));
dfr2 = tibble::tibble(z=z,s=yz);
dfr3 = tibble::tibble(z=xk,s=exp(lgtyk)/(1+exp(lgtyk)));
ggplot(dfr1,aes(x=z,y=s)) + geom_line() + geom_point() +
  geom_point(data=dfr2,colour="blue") +
  geom_point(data=dfr3,colour="red",fill=NA,shape=23,size=6)
compareSelFun(selSplineClmpdLeft,z,params,knots,title="selSplineClmpdLeft(z,params,knots)");

#--sel function: selSplineClmpdRight(z,p,knots)----
#--set up a spline curve "clamped" on the right end
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            0.5, #--scaled increment to params[1] at which descending limb departs from 1
            35); #--width of descending limb
yz = dblnormal4(z,params,refZ);
yz = ifelse(yz<0.0001,0.0001,ifelse(yz>0.9999,0.999,yz));
nz = length(z);
nk = 6; ik = seq(from=3,to=nz-3,length.out=nk);
xk = z[ik];
lgtyk = log(yz[ik]/(1-yz[ik]));
params = lgtyk;
knots = xk;
sf = selSplineClmpdRight(z,params,knots);
refZ  = 185;     #--not used
#--set `f` to function in global environment to be tested
dfr1 = tibble::tibble(z=z,s=selSpline(z,params,knots));
dfr2 = tibble::tibble(z=z,s=yz);
dfr3 = tibble::tibble(z=xk,s=exp(lgtyk)/(1+exp(lgtyk)));
ggplot(dfr1,aes(x=z,y=s)) + geom_line() + geom_point() +
  geom_point(data=dfr2,colour="blue") +
  geom_point(data=dfr3,colour="red",fill=NA,shape=23,size=6)
compareSelFun(selSplineClmpdRight,z,params,knots,title="selSplineClmpdRight(z,params,knots)");




