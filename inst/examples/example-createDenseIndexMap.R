#--create some example dimensions
require(tibble);
#----years
vYs = 2001:2005;                  attr(vYs,"dmnms")<-"y";
#----seasons
vCs = c("fall","spring");         attr(vCs,"dmnms")<-"c";
#----regions
vRs = c("All EBS");               attr(vRs,"dmnms")<-"r";
#----sexes
vXs = c("MALE","FEMALE");         attr(vXs,"dmnms")<-"x";
#----maturity state nested in sex (working outward)
vMs = list(  MALE=c("IMMATURE","MATURE"),
           FEMALE=c("IMMATURE","MATURE"));
attr(vMs,"dmnms")<-c("x","m");
#----shell condition nested in maturity state/sex (working outward)
vSs = list(MALE=  list(IMMATURE=c("NEW SHELL"),
                        MATURE=c("NEW SHELL","OLD SHELL")),
          FEMALE=list(IMMATURE=c("NEW SHELL"),
                        MATURE=c("NEW SHELL","OLD SHELL")));
attr(vSs,"dmnms")<-c("x","m","s");
#----ages nested in shell condition/maturity state/sex (working outward)
vAs = list(MALE=  list(IMMATURE=list(`NEW SHELL`="NO AGEs"),
                         MATURE=list(`NEW SHELL`=c("PTM 0"),
                                     `OLD SHELL`=c("PTM 1","PTM 2"))),
          FEMALE=list(IMMATURE=list(`NEW SHELL`="NO AGEs"),
                        MATURE=list(`NEW SHELL`=c("PTM 0"),
                                    `OLD SHELL`=c("PTM 1","PTM 2"))));
attr(vAs,"dmnms")<-c("x","m","s","a");
#----sizes nested in age/shell condition/maturity state/sex (working outward)
vZs = list( MALE=  list(IMMATURE=list(`NEW SHELL`=list(`NO AGEs`=seq(40,60,5))),
                          MATURE=list(`NEW SHELL`=list(`PTM 0`=seq(80,120,10)),
                                      `OLD SHELL`=list(`PTM 1`=seq(80,120,10),
                                                       `PTM 2`=seq(80,120,10)))),
          FEMALE=list(IMMATURE=list(`NEW SHELL`=list(`NO AGEs`=seq(25,50,5))),
                        MATURE=list(`NEW SHELL`=list(`PTM 0`=seq(50,100,10)),
                                    `OLD SHELL`=list(`PTM 1`=seq(50,100,10),
                                                     `PTM 2`=seq(50,100,10)))));
attr(vZs,"dmnms")<-c("x","m","s","a","z");

#--create maps from the 1-d index to the dimension levels for various combinations
#----of the example dimensions.
require(rtmbGMACS);
dfrSparse = createSparseDimsMap(y=vYs);   #--defined combinations
dfrFull   = createDenseDimsMap(dfrSparse); #--all possible combinations
dfrSparse = createSparseDimsMap(y=vYs,c=vCs);
dfrFull   = createDenseDimsMap(dfrSparse); #--all possible combinations
dfrSparse = createSparseDimsMap(y=vYs,z=vZs);
dfrFull   = createDenseDimsMap(dfrSparse); #--all possible combinations
dfrSparse = createSparseDimsMap(a=vAs);
dfrFull   = createDenseDimsMap(dfrSparse); #--all possible combinations
dfrSparse = createSparseDimsMap(z=vZs);
dfrFull   = createDenseDimsMap(dfrSparse); #--all possible combinations
dfrSparse = createSparseDimsMap(a=vAs,z=vZs);
dfrFull   = createDenseDimsMap(dfrSparse); #--all possible combinations
dfrSparse = createSparseDimsMap(y=vYs,c=vCs,r=vRs,x=vXs,m=vMs,s=vSs,a=vAs,z=vZs);
dfrFull   = createDenseDimsMap(dfrSparse); #--all possible combinations
