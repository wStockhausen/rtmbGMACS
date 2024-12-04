#--script to test the RTMB GMACS framework and the objective function
#--
require(tibble);
#--
dirPrj = rstudioapi::getActiveProject();

#--create inputs list
#----this list includes the following elements:
#------model configuration information ("cfg")
#------model control information ("ctl")
#------data to be fit ("dat")
#-----model optimization information ("fit")

inputs<-list();

#----create the configuration information
cfg = list();
#----define the dimensions, create the dimension maps
if (!requireNamespace("rtmbGMACS")){
  source(file.path(dirPrj,"R","DimensionsUtilities.R"));
  source(file.path(dirPrj,"R","DimensionsFunctions.R"));
}
#----years
vYs = 2001:2005;  attr(vYs,"dmnms")<-"y";
TB1 = 2001:2003; #--time block 1
TB2 = 2004:2005; #--time block 2
#----seasons
vCs = 1:4;        attr(vCs,"dmnms")<-"c";
#----sizes nested in age/shell condition/maturity state/sex (working outward)
vZs = list( MALE=  list(IMMATURE=list(`NEW SHELL`=seq(40,60,5),
                                      `OLD SHELL`=seq(40,60,5)),
                          MATURE=list(`NEW SHELL`=seq(80,120,10),
                                      `OLD SHELL`=seq(80,120,10))),
          FEMALE=list(IMMATURE=list(`NEW SHELL`=seq(25,50,5),
                                    `OLD SHELL`=seq(25,50,5)),
                        MATURE=list(`NEW SHELL`=seq(50,100,10),
                                    `OLD SHELL`=seq(50,100,10))));
attr(vZs,"dmnms")<-c("x","m","s","z");

lst = createDimsMaps(y=vYs,c=vCs,z=vZs);
dmsAllS2D = lst$dfrS2D;#--sparse-to-dense map
dmsAllD2S = lst$dfrD2S;#--dense-to-sparse map
lst = createDimsMaps(z=vZs);
dmsPopS2D = lst$dfrS2D;#--sparse-to-dense map
dmsPopD2S = lst$dfrD2S;#--dense-to-sparse map
rm(lst);
cfg[["dms"]]=list(All=list(S2D=dmsAllS2D,D2S=dmsAllD2S),
                  Pop=list(S2D=dmsPopS2D,D2S=dmsPopD2S));

#--define the fleets
flts=c("TCF"=list(type="fishery",season=2,timing="end"),
       "NMFS"=list(type="survey",season=1,timing="start"));
cfg[["fleets"]] = flts;

#----create the parameter info
#--TODO: env covars, random effects on function or parameters?, additive or multiplicative?

#----define recruitment parameters
#------devs could be "random effects" on lnRbar, with variance lnSigR
#------what about time blocks?
#------define season to apply? --> at "fcn" level
#------define "group" and sie into which to recruit
#------how to define/apply mirroring? --> use negative id to refer to mirrored quantity
#----R_t(y) = exp(lnRbar_t + lnDevsR_t(y) + sum_i[lnRecECs[i]*env_cov[i]]);
#----lnDevsR_t(y)~dnorm(0,lnSigR_t) or lnDevsR_t(y)~AR1(rhoRec_t, lnSigR_t)
recInfo  = list();
recInfo[["fcns"]]    = list(`1`=list(id=1,type="noSR",season="last",dms=list(y=vYs),params=c("lnRbar","lnDevsR","lnSigR")));
recInfo[["lnRbar"]]  = list(`1`=list(fcn=1,id=1,type="scalar",iv=10,lb=5,ub=15,phz=1,jtr=TRUE,rdm=FALSE,pr="none",p1=NA,p2=NA));
recInfo[["lnSigR"]]  = list(`1`=list(fcn=1,id=1,type="scalar",iv= 0,lb=-5,ub=5,phz=1,jtr=TRUE,rdm=FALSE,pr="none",p1=NA,p2=NA));
recInfo[["lnDevsR"]] = list(`1`=list(fcn=1,id=1,type="vector",dfr=tibble::tibble(y=vYs,iv= 0,lb=-5,ub=5,phz=1,jtr=TRUE,rdm=FALSE,pr="normal",p1=0,p2=2)));
recInfo[["prZsRec"]] =

#--define mortality parameters (g indicates population "group": e.g., xms)
#----lnM_tg(y) = lnM_tg + lnMo_tg + lnDevsM_tg(y) + sum_i[lnECsM_tg[i]*env_cov[i]]);
#----lnDevsM_tg(y)~dnorm(0,lnSigM_tg) or lnDevsR_tg(y)~AR1(rhoM_tg, lnSigM_tg);
#----lnM_tgz(y) = Lorenzen function for size variation
nmInfo = list();
nmInfo[["fcns"]] = list(`1`=list(id=1,type="lnM",dms=list(y=TB1,g=)))









