#--testing mgcv smoothers in RTMB
#--smooth specification
require(mgcv);
require(ggplot2);
dirPrj = here::here();
source(file.path(dirPrj,"R/SelectivityFunctions.R"));

#--create model dimensions----
source(file.path(dirPrj,"R","MiscFunctions_Text.R"));
source(file.path(dirPrj,"R","DimensionsFunctions.R"));
source(file.path(dirPrj,"R","DimensionsUtilities.R"));
## The following string would be in the "dims"" file:
str=paste(
'MODEL_DIMS
  y <- 1948:2024;                           #--years
  s <- 1;                                   #--seasons
  r <- "EBS";                               #--regions
  x <- c("male","female");                  #--sex classes
  m <- c("immature","mature");              #--maturity state classes
  p <- c("new_shell","old_shell");          #--post-molt ages
  zc <- seq(24.5,184.5,5);                  #--size bin cutpoints
  z  <- 0.5*(zc[2:(length(zc))]+zc[1:(length(zc)-1)]); #--size bin midpoints
END')
strv = stringr::str_split_1(str,"\\n") |> extractLines("MODEL_DIMS","END");
lstDims = parseStrAsList(strv,verbose=TRUE);
dimsMdl = createSparseDimsMap(!!!(lstDims[names(lstDims)!="zc"]));#--drop zc (size bin cutpoints)
dimsZBC = tibble::tibble(lft=lstDims$zc[1:(length(lstDims$zc)-1)],#--size bin dims (left, center,right)
                         mid=lstDims$z,
                         rgt=lstDims$zc[2:(length(lstDims$zc))]);

#--create simple smooth specification, smooth construct, and model matrix
dfrZ = dimsMdl |> dplyr::select(z) |> dplyr::distinct(); #--dataframe for smooth construct
sso = eval(parse(text="s(z,k=10)"));             #--smooth specification
sc  = mgcv::smoothCon(sso,data=dfrZ,knots=NULL); #--smooth construct
Xf  = sc[[1]]$X;                                 #--model matrix for (first) smooth term

tmp = dfrZ |> dplyr::mutate(v=asclogistic(dfrZ$z,c(100,0.05)));
bta = coef(lm(v~Xf-1,data=tmp));
tmp$p = (Xf %*%bta)[,1];
ggplot(tmp |> tidyr::pivot_longer(c(v,p)),aes(x=z,y=value,colour=name)) +
  geom_line() + geom_point();

#--use "by" variable
dfrYZ = dimsMdl |> dplyr::select(y) |> dplyr::distinct() |>
          dplyr::mutate(y=factor(as.character(y)),
                        z50=rnorm(dplyr::n(),100,10)) |>
          dplyr::cross_join(dimsMdl |> dplyr::select(z) |> dplyr::distinct()) |>
          dplyr::mutate(rp=z50,
                        v=asclogistic(z,pZ50=z50,pSlp=0.05)); #--dataframe for smooth construct (as well as data to fit)
sso = eval(parse(text="s(z,k=10,by=y)"));         #--smooth specification
sc  = mgcv::smoothCon(sso,data=dfrYZ,knots=NULL); #--smooth construct
Xf  = Matrix::Matrix(sc[[1]]$X);                  #--sparse model matrix for first smooth term
i   = 2;
while(i<length(sc)){
  Xf = cbind(Xf,Matrix::Matrix(sc[[i]]$X));        #--create model matrix by appending "by" model matrices
  i  = i+1;
}

#--fit using lm----
Xf = as.matrix(Xf);
bta = coef(lm(v~Xf-1,data=dfrYZ));
tmp = dfrYZ;
tmp$p = (Xf %*%bta)[,1];
ggplot(tmp |> tidyr::pivot_longer(c(v,p)) |> dplyr::filter(as.character(y)!="2024"),
  aes(x=z,y=value,colour=name,shape=name,group=paste(as.character(y),name))) +
  geom_path(aes(linetype=name)) + geom_point() +
  facet_wrap(~y);

View(tmp |> tidyr::pivot_longer(c(v,p)) |> dplyr::filter(as.character(y)=="2024"));

#--fit using RTMB----
require(Matrix);
require(TMB);
require(RTMB);
Xf  = Matrix::Matrix(sc[[1]]$X);                  #--sparse model matrix for first smooth term
i   = 2;
while(i<=length(sc)){
  Xf = cbind(Xf,Matrix::Matrix(sc[[i]]$X));       #--create sparse model matrix by appending "by" model matrices
  i  = i+1;
}
data = list(dfrYZ=dfrYZ,
            Xf=Xf);
#--create model/objective function----
objfn<-function(pars){
  p = pars$p;
  cat("p = \n",p[1:10],"\n\n");
  cat("pars$p = \n",pars$p[1:10],"\n\n");
  cat("str(p) = \n");cat(str(p),"\n\n");

  obs<-OBS(as.vector(data$dfrYZ$v));         #--observations
  cat("obs = \n",head(obs),"\n\n");
  cat("obs = \n"); cat(str(obs),"\n\n");

  Xf = AD(data$Xf);
  #cat("Xf = \n",head(as.matrix(Xf),n=c(10,10)),"\n\n");
  #cat("str(Xf) = \n"); cat(str(Xf),"\n\n")
  #Xf = setAs(data$Xf,"adsparse");
  #cat("Xf = \n",head(as.matrix(Xf),n=c(10,10)),"\n\n");
  cat("str(Xf) = \n"); cat(str(Xf),"\n\n")

  prd = Xf %*% p; #--predicted values
  # prd = t(p) %*% t(Xf) ; #--predicted values
  cat("str(prd) = \n"); cat(str(prd),"\n\n");

  cat("prd = \n",head(as.vector(prd)),"\n\n")
  ADREPORT(pars);

  REPORT(prd);

  # obs %~% dnorm(as.vector(prd),AD(1),log=TRUE); #--adds -log-likelihood to hidden variable `.nll`
  nll = AD(0);
  nll = nll-sum(dnorm(obs,as.vector(prd),1,log=TRUE));
  cat("nll = \n",nll,"\n\n");
  nll;
}
pars = list(p=as.double(rep(1,ncol(Xf))));
objfn(pars);
#--MakeADFun the model----
#cat("MakeADFun'ing\n");
mdl = RTMB::MakeADFun(objfn,parameters=pars);

