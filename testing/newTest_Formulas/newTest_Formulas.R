#--testing R package "Formula" capabilities
dirPrj =  rstudioapi::getActiveProject();
#dirPrj= getwd();
source(file.path(dirPrj,"R","MiscFunctions_ReadParse.R"));
source(file.path(dirPrj,"R","MiscFunctions_Text.R"));
source(file.path(dirPrj,"R","DimensionsFunctions.R"));
source(file.path(dirPrj,"R","DimensionsUtilities.R"));
source(file.path(dirPrj,"R","defineYearBlocks.R"));

vRs="EBS";
miz = as.character(seq(25,75,5));
mmz = as.character(seq(50,75,5));
fiz = as.character(seq(25,65,5));
fmz = as.character(seq(30,65,5));
vXs=list(male=list(imm=list(`new`=miz,
                             `old`=miz),
                    mat=list(`new`=mmz,
                             `old`=mmz,
                             `very`=mmz)
                  ),
          female=list(imm=list(`new`=fiz,
                               `old`=fiz),
                      mat=list(`new`=fmz,
                               `old`=fmz)
                     )
         ); attr(vXs,"dmnms")<-c("x","m","p","z");
ys = as.character(2015:2024); names(ys) = ys;
ss = as.character(1);         names(ss) = ss;
dmsC    = createSparseDimsMap(r=vRs,x=vXs);
dmsYS   = createSparseDimsMap(y=ys,s=ss);
dmsYSC  = createSparseDimsMap(y=ys,s=ss,r=vRs,x=vXs);
dimsZBC = dmsC |> dplyr::mutate(lft=as.numeric(z)-2.5,mid=as.numeric(z),rgt=as.numeric(z)+2.5);

#--define time block for NMFS Q's
lstNMFS<-list(preCV =as.character(c(2015:2019)),       #--pre gear change
              preGC =as.character(c(2021:2022)),       #--pre gear change
              postGC=as.character(c(2023:2025)));      #--post gear change
tbNMFS = lstNMFS |> defineYearBlocks();
##--expand to dimensions map for NMFS Q's
dmsNMFS<-tbNMFS |>
         defineYearBlocks(dimsMdl=dmsYSC) |>
         dplyr::mutate(f="NMFS");

#--main effects----
##--crossed effects----
f = "~1+yblk*m";
dfrMF = tbNMFS |> dplyr::cross_join(tibble::tribble(~m,"imm","mat")); #--model frame
Frmla = Formula::Formula(eval(parse(text=f)));
MdFrm = Formula:::model.frame.Formula(Frmla,data=dfrMF);
###--treatment and sum contrasts
ctrs = list(yblk="contr.treatment",m="contr.sum");
dfr   = convertToDimsFactors(MdFrm,NULL,contrasts=ctrs,verbose=TRUE); str(dfr);
MdMtx = Formula:::model.matrix.Formula(Frmla,data=dfr,rhs=NULL);
colnames(MdMtx)
head(dplyr::bind_cols(as.data.frame(MdMtx),MdFrm));
###--helmert and poly contrasts----
ctrs = list(yblk="contr.helmert",m="contr.poly");
dfr   = convertToDimsFactors(MdFrm,NULL,contrasts=ctrs,verbose=TRUE); str(dfr);
MdMtx = Formula:::model.matrix.Formula(Frmla,data=dfr,rhs=NULL);
colnames(MdMtx)
head(dplyr::bind_cols(as.data.frame(MdMtx),MdFrm));

#--smooths: NEED TO USE mgcv functions to do this
##--single s() smooth, no by variable
require(mgcv);
f = "s(z,k=10)";
sso = eval(parse(text=f));                       #--smooth specification formula
sc  = mgcv::smoothCon(sso,data=dmsNMFS |> dplyr::mutate(z=as.numeric(z)),knots=NULL); #--smooth construct
Xf  = sc[[1]]$X;                                 #--model matrix for (first) smooth term
Xf;
##--single s() smooth, with by variable
f = "s(z,k=10,by=yblk)";
sso = eval(parse(text=f));                       #--smooth specification formula
sc  = mgcv::smoothCon(sso,data=dmsNMFS |> dplyr::mutate(z=as.numeric(z)),knots=NULL); #--smooth construct
Xf  = sc[[1]]$X;                                 #--model matrix for (first) smooth term
Xf;
##--single ti() smooth with intercept, no by variable
require(mgcv);
f = "ti(z,k=10)";
sso = eval(parse(text=f));                       #--smooth specification formula
sc  = mgcv::smoothCon(sso,data=dmsNMFS |> dplyr::mutate(z=as.numeric(z)),knots=NULL); #--smooth construct
Xf  = sc[[1]]$X;                                 #--model matrix for (first) smooth term
Xf;
##--single ti() smooth, with by variable
f = "ti(z,k=10,by=yblk)";
sso = eval(parse(text=f));                       #--smooth specification formula
sc  = mgcv::smoothCon(sso,data=dmsNMFS |> dplyr::mutate(z=as.numeric(z)),knots=NULL); #--smooth construct
Xf  = sc[[1]]$X;                                 #--model matrix for (first) smooth term
Xf;
PredictMat(sc,data=(dmsNMFS |> dplyr::mutate(z=as.numeric(z))))

##--single ti() smooth, with factor smooth by-variable resulting in random wiggly curves
f = "ti(z,y,k=10,bs='fs')";
sso = eval(parse(text=f));                       #--smooth specification formula
sc  = mgcv::smoothCon(sso,data=dmsNMFS |> dplyr::mutate(z=as.numeric(z),y=factor(y)),knots=NULL); #--smooth construct
Xf  = sc[[1]]$X;                                 #--model matrix for (first) smooth term
Xf;

#--random effects
dfrMF = tbNMFS; #--model frame
f = "~0 + y|yblk";
ctrs = FALSE;
Frmla = Formula::Formula(eval(parse(text=f)));
MdFrm = Formula:::model.frame.Formula(Frmla,data=dfrMF);
MdMtx = Formula:::model.matrix.Formula(Frmla,data=convertToDimsFactors(MdFrm,NULL,contrasts=ctrs),rhs=NULL);
