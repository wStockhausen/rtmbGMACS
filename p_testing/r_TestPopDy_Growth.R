#--test Growth portion of population dynamics----
dirPrj = rstudioapi::getActiveProject();
if (FALSE){
  require(rtmbGMACS);
} else {
  source(file.path(dirPrj,"R","DimensionsFunctions.R"))
  source(file.path(dirPrj,"R","DimensionsUtilities.R"))
  source(file.path(dirPrj,"R","MiscFunctions.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Dataframe.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Transforms.R"))
  source(file.path(dirPrj,"R","readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R","readParamInfo_Growth_PrGr.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType2.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Growth_PrGr.R"))
  source(file.path(dirPrj,"R","calcGrowth_PrGr.R"))
  source(file.path(dirPrj,"R","calcPopDynamics.R"));
}

vRs="EBS";attr(vRs,"dmnms")<-"r";
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
ys = as.character(2020:2024); names(ys) = ys;
ss = as.character(1);         names(ss) = ss;
dmsC    = createSparseDimsMap(r=vRs,x=vXs);
dmsYS   = createSparseDimsMap(y=ys,s=ss);
dmsYSC  = createSparseDimsMap(y=ys,s=ss,r=vRs,x=vXs);

#--test maturation: imm -> mat----
mMPZ_MPZ = createMatrixNextMaturityState(dimMPZ<-dmsC |> dplyr::filter(r=="EBS",x=="female") |> dplyr::select(m,p,z));
mMPZ_MPZ = createMatrixNextMaturityState(dimMPZ<-dmsC |> dplyr::filter(r=="EBS",x=="male")   |> dplyr::select(m,p,z));

#--test molting: p -> p0-----
mMPZ_MPZ = createMatrixPto0(dimMPZ<-dmsC |> dplyr::filter(r=="EBS",x=="female") |> dplyr::select(m,p,z));
mMPZ_MPZ = createMatrixPto0(dimMPZ<-dmsC |> dplyr::filter(r=="EBS",x=="male")   |> dplyr::select(m,p,z));

#--test non-molting: p -> p+1----
mMPZ_MPZ = createMatrixPtoPplus1(dimMPZ<-dmsC |> dplyr::filter(r=="EBS",x=="female") |> dplyr::select(m,p,z));
mMPZ_MPZ = createMatrixPtoPplus1(dimMPZ<-dmsC |> dplyr::filter(r=="EBS",x=="male")   |> dplyr::select(m,p,z));

#--test doGrowth for one time step for one r, x----
##--read parameter infor file and create params list
conn  = file.path(dirPrj,"testing/testGrowth_PrGr/inputSpecs_Growth_PrGr.function.txt");
res   = readParamInfo_Growth_PrGr(conn,FALSE);
infoPrGr = extractParamInfo_Growth_PrGr(res,dims=list(dmsYSC=dmsYSC),FALSE);
params = list(pPrGr_MPs=lstPrGr$MPs$params);
lstMatsPrGr = calcGrowth_PrGr(dmsYSC,infoPrGr,params);

dimMPZ<-dmsC |> dplyr::filter(r=="EBS",x=="female") |> dplyr::select(m,p,z);
mPtoP0  = createMatrixPto0(dimMPZ);
mPtoPp1 = createMatrixPtoPplus1(dimMPZ);
mNextMatState=createMatrixNextMaturityState(dimMPZ);
n_mpz = 1:nrow(dimMPZ);
lst1 = list(vPrMolt=dimMPZ$m=="imm",
            mGrw,
            vPrMat=(dimMPZ$m=="imm")*1/(1+exp(-(as.numeric(dimMPZ$z)-max(as.numeric((dimMPZ |> dplyr::filter(m=="imm"))$z)))/1)),
            mPtoP0=mPtoP0,
            mPtoPp1=mPtoPp1,
            mNextMatState=mNextMatState);
np_mpz = doGrowth<(n_mpz,lst);
