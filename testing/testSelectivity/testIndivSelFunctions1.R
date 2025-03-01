#--script to test selectivity functions meant for RTMB
require(ggplot2);
require(RTMB);
dirPrj = rstudioapi::getActiveProject();
source(file.path(dirPrj,"R/MiscFunctions.R"),local=TRUE);
source(file.path(dirPrj,"R/SelectivityFunctions1.R"),local=TRUE);
source(file.path(dirPrj,"R/getDFR_sdreport.R"),local=TRUE);
source(file.path(dirPrj,"testing/testSelectivity/compareSelFun1.R"),local=TRUE);

#--define size bins
z = seq(25,180,by=5);

#--sel function: const_sel(z,p,pRefZ)----
#----function should not really be differentiable(?)
pCnst = c(1); #--constant value
pRefZ = 100;  #--reference *size* (not used!)
#--set `f` to function in global environment to be tested
res = compareSelFun1("const_sel",z,
                     pCnst=pCnst,pRefZ=pRefZ,
                     map_=list(pRefZ=factor(NA)),
                     title="const_sel(z,pCnst,pRefZ)",
                     verbose=TRUE);

#--sel function: asclogistic(z,p,pRefZ)----
z50 = 100.00; #--size at which selectivity = 0.50 (logit-scale mean)
slp =   0.05; #--slope at z50
pRefZ  = 125;  #--reference *size*
#--set `f` to function in global environment to be tested
res = compareSelFun1("asclogistic",z,
                     pZ50=z50,pSlp=slp,pRefZ=pRefZ,
                     map_=list(pRefZ=factor(NA)),
                     title="asclogistic(z,z50,slp,pRefZ)",verbose=TRUE);

#--sel function: asclogistic1(z,p,pRefZ)----
z50 = 100;   #--size at which selectivity = 0.50 (logit-scale mean)
wdZ = 20; #--width (1/slope) at z50
pRefZ  = 125;       #--reference *size*
#--set `f` to function in global environment to be tested
res = compareSelFun1("asclogistic1",z,
                     pZ50=z50,pWdZ=wdZ,pRefZ=pRefZ,
                     map_=list(pRefZ=factor(NA)),
                     title="asclogistic1(z,p,pRefZ)",
                     verbose=TRUE);

#--sel function: asclogistic5095(z,p,pRefZ)----
z50 = 50;   #--size at which selectivity = 0.50 (logit-scale mean)
z95 = 125;  #--size at which selectivity = 0.95
pRefZ  = 0; #--reference *size*
#--set `f` to function in global environment to be tested
res = compareSelFun1("asclogistic5095",z,
                     pZ50=z50,pZ95=z95,pRefZ=pRefZ,
                     map_=list(pRefZ=factor(NA)),
                     title="asclogistic5095(z,p,pRefZ)",
                     verbose=TRUE);

#--sel function: asclogistic50D95(z,p,pRefZ)----
z50   = 50;  #--size at which selectivity = 0.5 (logit-scale mean)
z9550 = 75;  #--difference between sizes at 95\% and 50\%-selected
pRefZ  = 0;       #--reference *size*
#--set `f` to function in global environment to be tested
res = compareSelFun1("asclogistic50D95",z,
                     pZ50=z50,pZ9550=z9550,pRefZ=pRefZ,
                     map_=list(pRefZ=factor(NA)),
                     title="asclogistic50D95(z,p,pRefZ)",
                     verbose=TRUE);

#--sel function: ascnormal1(z,p)----
ascMnZ = 100;  #--size at which ascending limb reaches 1
ascWdZ =  25;  #--width of ascending limb
#--set `f` to function in global environment to be tested
res = compareSelFun1("ascnormal1",z,
                     pAscZ1=ascMnZ,pAscWdZ=ascWdZ,
                     title="ascnormal(z,p,pRefZ)",
                     verbose=TRUE);

#--sel function: ascnormal2(z,p,pRefZ)----
ascMnZ  = 100;   #--size at which ascending limb reaches 1
ascRefS = 0.5;   #--selectivity at size = pRefZ
pRefZ   = 50;    #--reference *size*
#--set `f` to function in global environment to be tested
res = compareSelFun1("ascnormal2",z,
                     pAscZ1=ascMnZ,pAscRefS=ascRefS,pRefZ=pRefZ,
                     map_=list(pRefZ=factor(NA)),
                     title="ascnormal2(z,p,pRefZ)",
                     verbose=TRUE);

#--sel function: ascnormal2a(z,p,refS)----
pAscZ1  = 100; #--size at which ascending limb reaches 1
pZatRefS = 50;  #--size at sel = refS
pRefS    = 0.5; #--reference *selectivity*
#--set `f` to function in global environment to be tested
res = compareSelFun1("ascnormal2a",z,
                     pAscZ1=pAscZ1,pZatRefS=pZatRefS,pRefS=pRefS,
                     map_=list(pRefS=factor(NA)),
                     title="ascnormal2a(z,p,refS)",
                     verbose=TRUE);

#--sel function: ascnormal2b(z,p,refS)----
pAscZ1   = 150;  #--size at which ascending limb reaches 1
pDZ2RefS =  75;  #--delta from size at 1 to size at which selectivity=refS
pRefS    = 0.5;  #--reference *selectivity*
#--set `f` to function in global environment to be tested
res = compareSelFun1("ascnormal2b",z,
                     pAscZ1=pAscZ1,pDZ2RefS=pDZ2RefS,pRefS=pRefS,
                     map_=list(pRefS=factor(NA)),
                     title="ascnormal2b(z,p,refS)",
                     verbose=TRUE);

#--sel function: ascnormal3(z,p,refS)----
pDZ1   = 50;   #--delta from max possible size (pRefZ[1]) at which ascending limb could reach 1
pSatZ2 =  0.5; #--selectivity at size = pRefZ2
pMxZ1  = 185;  #--max possible size at which the curve could reach 1 (e.g., max(z))
pRefZ2 =  75;   #--reference *size* at which curve reaches pSatZ2
#--set `f` to function in global environment to be tested
res = compareSelFun1("ascnormal3",z,
                     pDZ1=pDZ1,pSatZ2=pSatZ2,pMxZ1=pMxZ1,pRefZ2=pRefZ2,
                     map_=list(pMxZ1=factor(NA),
                               pRefZ2=factor(NA)),
                     title="ascnormal3(z,p)",
                     verbose=TRUE);

#--sel function: dbllogistic5095(z,p,pRefZ)----
params = c( 50,  #--ascending limb size at 50% selected
           100,  #--ascending limb size at 95% selected
           125,  #--descending limb size at 95% selected
           150); #--descending limb size at 50% selected
pRefZ  = 0;       #--reference *size*
#--set `f` to function in global environment to be tested
res = compareSelFun1(dbllogistic5095,z,
                     params,pRefZ,
                     title="dbllogistic5095(z,p,pRefZ)");

#--sel function: dbllogistic(z,p,pRefZ)----
params = c( 50,  #--ascending limb size at which selectivity = 0.5 (logit-scale mean)
           0.5,  #--ascending limb slope at 50% selected
           125,  #--descending limb size at which selectivity = 0.5 (logit-scale mean)
           0.5); #--descending limb size at 50% selected
pRefZ  = 0;       #--reference *size*
#--set `f` to function in global environment to be tested
res = compareSelFun1(dbllogistic,z,
                     params,pRefZ,
                     title="dbllogistic(z,p,pRefZ)");

#--sel function: dblnormal4(z,p,pRefZ)----
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            40,   #--offset to size at which descending limb departs from 1
            35); #--width of descending limb
pRefZ  = 0;       #--reference *size*
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,s=dblnormal4(z,params,pRefZ)),aes(x=z,y=s)) + geom_line();
res = compareSelFun1(dblnormal4,z,
                     params,pRefZ,
                     title="dblnormal4(z,p,pRefZ)");

#--sel function: dblnormal4a(z,p,pRefZ)----
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            0.5, #--scaled increment to params[1] at which descending limb departs from 1
            35); #--width of descending limb
pRefZ  = 185;       #--reference *size*: max possible size (e.g., max(z))
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,s=dblnormal4a(z,params,pRefZ)),aes(x=z,y=s)) + geom_line();
res = compareSelFun1(dblnormal4a,z,
                     params,pRefZ,
                     title="dblnormal4a(z,p,pRefZ)");

#--sel function: dblnormal6(z,p,pRefZ)----
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
           125,  #--size at which descending limb departs from 1
            35,  #--width of descending limb
           0.1,  #--floor of ascending limb
           0.1); #--floor of descending limb
pRefZ  = 185;       #--reference *size*: max possible size (e.g., max(z))
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,s=dblnormal6(z,params,pRefZ)),aes(x=z,y=s)) + geom_line();
res = compareSelFun1(dblnormal6,z,
                     params,pRefZ,
                     title="dblnormal6(z,p,pRefZ)");

#--sel function: stackedLogistic1(z,p,pRefZ)----
params = c(0.2,  #--omega: -weighting factor on the first curve
            75,  #--mnZ1: size at inflection point for 1st logistic curve
            10,  #--sdZ1: sd for 1st logistic curve
           145,  #--mnZ2: size at inflection point for 2nd logistic curve
            10); #--sdZ2: sd for 2nd logistic curve
pRefZ  = 185;     #--not used
#--set `f` to function in global environment to be tested
ggplot(tibble::tibble(z=z,s=stackedLogistic1(z,params,pRefZ)),aes(x=z,y=s)) + geom_line();
res = compareSelFun1(stackedLogistic1,z,
                     params,pRefZ,
                     title="stackedLogistic1(z,p,pRefZ)");

#--sel function: selSpline(z,p,knots)----
#--set up a spline curve
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            0.5, #--scaled increment to params[1] at which descending limb departs from 1
            35); #--width of descending limb
pRefZ  = 185;       #--reference *size*: max possible size (e.g., max(z))
yz = dblnormal4(z,params,pRefZ);
yz = ifelse(yz<0.0001,0.0001,ifelse(yz>0.9999,0.999,yz));
nz = length(z);
nk = 6; ik = seq(from=1,to=nz,length.out=nk);
xk = z[ik];
lgtyk = log(yz[ik]/(1-yz[ik]));
params = lgtyk;
knots = xk;
sf = selSpline(z,params,knots);
pRefZ  = 185;     #--not used
#--set `f` to function in global environment to be tested
dfr1 = tibble::tibble(z=z,s=selSpline(z,params,knots));
dfr2 = tibble::tibble(z=z,s=yz);
dfr3 = tibble::tibble(z=xk,s=exp(lgtyk)/(1+exp(lgtyk)));
ggplot(dfr1,aes(x=z,y=s)) + geom_line() + geom_point() +
  geom_point(data=dfr2,colour="blue") +
  geom_point(data=dfr3,colour="red",fill=NA,shape=23,size=6)
res = compareSelFun1(selSpline,z,
                     params,knots,
                     title="selSpline(z,params,knots)");

#--sel function: selSplineClmpd(z,p,knots)----
#--set up a spline curve "clamped" on both ends
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            0.5, #--scaled increment to params[1] at which descending limb departs from 1
            35); #--width of descending limb
yz = dblnormal4(z,params,pRefZ);
yz = ifelse(yz<0.0001,0.0001,ifelse(yz>0.9999,0.999,yz));
nz = length(z);
nk = 6; ik = seq(from=3,to=nz-3,length.out=nk);
xk = z[ik];
lgtyk = log(yz[ik]/(1-yz[ik]));
params = lgtyk;
knots = xk;
sf = selSplineClmpd(z,params,knots);
pRefZ  = 185;     #--not used
#--set `f` to function in global environment to be tested
dfr1 = tibble::tibble(z=z,s=selSpline(z,params,knots));
dfr2 = tibble::tibble(z=z,s=yz);
dfr3 = tibble::tibble(z=xk,s=exp(lgtyk)/(1+exp(lgtyk)));
ggplot(dfr1,aes(x=z,y=s)) + geom_line() + geom_point() +
  geom_point(data=dfr2,colour="blue") +
  geom_point(data=dfr3,colour="red",fill=NA,shape=23,size=6)
res = compareSelFun1(selSplineClmpd,z,
                     params,knots,
                     title="selSplineClmpd(z,params,knots)");

#--sel function: selSplineClmpdLeft(z,p,knots)----
#--set up a spline curve "clamped" on the left end
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            0.5, #--scaled increment to params[1] at which descending limb departs from 1
            35); #--width of descending limb
yz = dblnormal4(z,params,pRefZ);
yz = ifelse(yz<0.0001,0.0001,ifelse(yz>0.9999,0.999,yz));
nz = length(z);
nk = 6; ik = seq(from=3,to=nz-3,length.out=nk);
xk = z[ik];
lgtyk = log(yz[ik]/(1-yz[ik]));
params = lgtyk;
knots = xk;
sf = selSplineClmpdLeft(z,params,knots);
pRefZ  = 185;     #--not used
#--set `f` to function in global environment to be tested
dfr1 = tibble::tibble(z=z,s=selSpline(z,params,knots));
dfr2 = tibble::tibble(z=z,s=yz);
dfr3 = tibble::tibble(z=xk,s=exp(lgtyk)/(1+exp(lgtyk)));
ggplot(dfr1,aes(x=z,y=s)) + geom_line() + geom_point() +
  geom_point(data=dfr2,colour="blue") +
  geom_point(data=dfr3,colour="red",fill=NA,shape=23,size=6)
res = compareSelFun1(selSplineClmpdLeft,z,
                     params,knots,
                     title="selSplineClmpdLeft(z,params,knots)");

#--sel function: selSplineClmpdRight(z,p,knots)----
#--set up a spline curve "clamped" on the right end
params = c(100,  #--size at which ascending limb reaches 1
            50,  #--width of ascending limb
            0.5, #--scaled increment to params[1] at which descending limb departs from 1
            35); #--width of descending limb
yz = dblnormal4(z,params,pRefZ);
yz = ifelse(yz<0.0001,0.0001,ifelse(yz>0.9999,0.999,yz));
nz = length(z);
nk = 6; ik = seq(from=3,to=nz-3,length.out=nk);
xk = z[ik];
lgtyk = log(yz[ik]/(1-yz[ik]));
params = lgtyk;
knots = xk;
sf = selSplineClmpdRight(z,params,knots);
pRefZ  = 185;     #--not used
#--set `f` to function in global environment to be tested
dfr1 = tibble::tibble(z=z,s=selSpline(z,params,knots));
dfr2 = tibble::tibble(z=z,s=yz);
dfr3 = tibble::tibble(z=xk,s=exp(lgtyk)/(1+exp(lgtyk)));
ggplot(dfr1,aes(x=z,y=s)) + geom_line() + geom_point() +
  geom_point(data=dfr2,colour="blue") +
  geom_point(data=dfr3,colour="red",fill=NA,shape=23,size=6)
res = compareSelFun1(selSplineClmpdRight,z,
                     params,knots,
                     title="selSplineClmpdRight(z,params,knots)");




