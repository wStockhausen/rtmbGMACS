#--test parameterization using lme4 approach to specifications
require(lme4);

dirPrj = file.path(getwd(),"../.."); #--rstudioapi::getActiveProject();
source(file.path(dirPrj,"R","DimensionsUtilities.R"),local=TRUE);
source(file.path(dirPrj,"R","DimensionsFunctions.R"),local=TRUE);
source(file.path(dirPrj,"R","ragged_array.R"),local=TRUE);

#--define some dimensions
dZs = list( MALE=  list(IMMATURE=list(`NEW SHELL`=seq(25,180,5),
                                      `OLD SHELL`=seq(25,180,5)),
                          MATURE=list(`NEW SHELL`=seq(25,180,5),
                                      `OLD SHELL`=seq(25,180,5))),
          FEMALE=list(IMMATURE=list(`NEW SHELL`=seq(25,130,5),
                                    `OLD SHELL`=seq(25,130,5)),
                        MATURE=list(`NEW SHELL`=seq(25,130,5),
                                    `OLD SHELL`=seq(25,130,5))));
attr(dZs,"dmnms")<-c("x","m","s","z");
dfrZs = createSparseDimsMap(z=dZs);
dfrYZs = createSparseDimsMap(y=c(`2020`=2020,`2021`=2021,`2022`=2022,`2023`=2023),
                             z=dZs);

f = sparse_idx~1+x+(1+m|y);
glF = lme4::glFormula(f,dfrYZs);

reX = as.matrix(glF$reTrms$Zt)
View(reX)
