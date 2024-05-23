#--Model Dimensions for Tanner crab
#--regional dimensions
dRs = "all EBS";          attr(dRs,"dmnms")<-"r";
#--temporal dimensions
dYs = 1982:2022;          attr(dYs,"dmnms")<-"y";
dCs = c(0.62,0.01,0.37);  attr(dCs,"dmnms")<-"c";

#--biological dimensions
#----sizes nested in age/shell condition/maturity state/sex (working outward)
dZs = list( MALE=  list(IMMATURE=list(`NEW SHELL`=seq(25,180,5),
                                      `OLD SHELL`=seq(25,180,5)),
                          MATURE=list(`NEW SHELL`=seq(25,180,5),
                                      `OLD SHELL`=seq(25,180,5))),
          FEMALE=list(IMMATURE=list(`NEW SHELL`=seq(25,180,5),
                                    `OLD SHELL`=seq(25,180,5)),
                        MATURE=list(`NEW SHELL`=seq(25,180,5),
                                    `OLD SHELL`=seq(25,180,5))));
attr(dZs,"dmnms")<-c("x","m","s","z");

#--fleets
dFs = c("TCF","SCF","RKF","GFA","GFF","GFT","NMFS"); attr(dFs,"dmnms")<-"f";

#--additional time blocks
lstTBs = list(`1`=list(years=1982:2004,label="pre-rationalization"),
              `2`=list(years=2005:2022,label="post-rationalization"));
