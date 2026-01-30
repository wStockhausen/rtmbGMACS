#--functions for covariance/precision matrix calculations

#' Build a template for a diagonal precision matrix
#' @param np - number of (RE) parameters per group
#' @param ng - #--number of groups (levels) in grouping factor
#' @return a sparse diagonal matrix `T` of type Matrix::dgCMatrix
#' @details Parameter ordering is considered to be within group, then by group.
#' @export
#' @examplesIf FALSE
#' np = 2; #--number of (RE) parameters/group
#' ng = 5;  #--number of groups (levels) in grouping factor
#' isigma = rep(1:np,ng);
#' T = Matrix::.sparseDiagonal(q,x=isigma,shape="g"); #--template for Q
#' Tp = buildTemplateQ.diag(np,ng);
#' T-Tp
#'
buildTemplateQ.diag<-function(np,ng){
  isigma = rep(1:np,ng);
  Qt = Matrix::.sparseDiagonal(np*ng,x=isigma,shape="g"); #--template for Q
  return(Qt);
}

#' Calculate the values to fill a diagonal precision matrix (`Q`)
#' @param covpars - vector of input covariance parameters
#' @param idx - index vector mapping `covpars` values to `Q`
#' @return vector of values (inverse variances) for the "x" slot in `Q`
#' @details The values in `covpars` are standard deviations (\code{sigma}s) on the ln-scale.
#' These are converted to variances before being mapped to `Q`. Covariance parameters are ordered
#' by RE parameters within each group.
#'
#' @examplesIf FALSE
#' # example code
#' ng = 5;  #--number of groups (levels) in grouping factor
#' np = 2; #--number of (RE) parameters/group
#' #--number of covariance parameters/group = 2 (1 variance parameter for each RE parameter)
#' Qt = buildTemplateQ.diag(np,ng); #--template for actual Q
#' idx = Qt@x;                      #--"map" from covpars to Q
#' ##--R context
#' covpars = log(c(1,0.75))/2; #--std dev parameters on ln-scale (1 per RE parameter)
#' Q = Qt;
#' Q@x = fillQ.diag(covpars,idx);
#' Q
#' ##--RTMB context
#' covpars = RTMB::AD(log(c(1,0.5))/2,force=TRUE); #--advector
#' Q = RTMB::AD(Qt,force=TRUE);   #--adsparse matrix
#' Q@x = fillQ.diag(covpars,idx);
#' as.matrix(Q);
#' @export
#'
fillQ.diag<-function(covpars,idx){exp(-2*covpars)[idx];}

#' Build a template for a precision matrix with AR1 structure
#' across parameters within a group.
#' @param np - number of (RE) parameters per group
#' @param ng - number of groups (levels) in grouping factor
#' @return a sparse matrix `T` of type Matrix::dgCMatrix
#' @details In `T`, parameter ordering is considered to be within group, then by group.
#' If `sequential`=TRUE, the related covariance parameter vector should be ordered as
#' c(sigma's, rho's) (i.e., sequentially). If `sequential`=FALSE, sigma, rho values should be
#' interleaved as c(sigma[1],rho[1],sigma[2],rho[2],...).
#'
#' @export
#' @examplesIf FALSE
#' ng = 2; #--number of groups (levels) in grouping factor
#' np = 5; #--number of (RE) parameters
#' #--number of covariance parameters is 2/group (i.e., sigma, rho)
#' T = buildTemplateQ.ar1(np,ng);
#' T
#' T@x; #--covariance parameter vector indices
#'
#' T = buildTemplateQ.ar1(np,ng,sequential=FALSE);
#' T
#' T@x; #--covariance parameter vector indices
#'
#' @export
#'
buildTemplateQ.ar1<-function(np,ng,sequential=TRUE){
  lstQ = list();
  for (g in 1:ng){
    #--g = 1;
    if (sequential){
      isigma = rep(g,   each=np);
      irho   = rep(ng+g,each=(np-1));
    } else {
      isigma = rep(2*g-1,each=np);
      irho   = rep(2*g,  each=np);
    }
    # cat("ar1 gth npars=",np,"\n",
    #     "isigma:",isigma,"\n",
    #     "irho  :",irho,"\n")
    Qp = Matrix::bandSparse(np,k=c(0,1),diagonals=list(isigma,irho),sym=TRUE);#--looks ok for RTMB (via RTMB::AD(Qp));
    lstQ[[g]] = Qp;
  }
  Q = buildBlockDiagonalMat(lstQ);
  return(as(Q,"generalMatrix"));
}

#' Calculate the values to fill a precision matrix (Q) with an AR1 structure
#' across parameters within a group.
#' @param covpars - vector of input covariance parameters
#' @param ng - number of groups (levels) in grouping factor
#' @param idx - index vector mapping `covpars` values to Q
#' @return vector of values for the "x" slot in Q
#' @details The values in `covpars` are standard deviations (\code{sigma}s) on the ln-scale and
#' correlation parameters on the symmetric logit scale (-Inf,Inf) corresponding to (-1,1) on
#' the arithmetic scale. These are converted to variances and correlation s before being mapped to Q.
#'
#' The ordering of values within `covapar` is s1,c1,s2,c2,..., where s1 is sigma for the first group,
#' c is the correlation within the first group.
#'
#' @examplesIf FALSE
#' # example code
#' ng = 2; #--number of groups (levels) in grouping factor
#' np = 5; #--number of (RE) parameters
#' #--number of covariance parameters is 2/group (i.e., sigma, rho)
#' Qt = buildTemplateQ.ar1(np,ng,sequential=TRUE); #--template for actual Q
#' idx = Qt@x; #--"map" from covpars to Q
#' ##--R context
#' ###--because `sequential`=TRUE when building `Qt`, the covariance parameters are
#' ###--ordered as c(sigma's, rho's).
#' covpars = c(log(c(1,0.75))/2,symlogit(c(-0.5,0.5)));
#' Q = Qt;
#' Q@x = fillQ.ar1(covpars,ng,idx);
#' Q
#' ##--RTMB context
#' covpars = RTMB::AD(c(log(c(1,0.75))/2,symlogit(c(-0.5,0.5))),force=TRUE); #--advector
#' Q = RTMB::AD(Qt,force=TRUE);   #--adsparse matrix
#' Q@x = fillQ.ar1(covpars,ng,idx);
#' as.matrix(Q);
#' @export
#'
fillQ.ar1 <-function(covpars,ng,idx,sequential=TRUE){
  if (sequential) {
    v = c(exp(-2*covpars[1:ng]),symlogistic(covpars[(ng+1):(2*ng)]))[idx];
  } else {
    #--fill in
  }
  return(v)
}

#' @title Convert a precision matrix to a covariance matrix
#' @param Q - precision matrix
#' @return covariance matrix
#' @details Uses [solve] to convert the precision matrix to a covariance matrix. If
#' Q is a sparse matrix class, the resulting covariance matrix may not be.
#' @examplesIf FALSE
#' # example code
#' ##--R context
#' Q = mkPrecMat(5,c(0,0.5),"ar1",sparse=TRUE);
#' C = convertPrecMatToCovMat(Q);
#' t = getThetaFromCovMat(C);
#' Cp = getCorrMatFromTheta(t);
#' Qp = convertCovMatToPrecMat(Cp);
#'
#' ##--RTMB context
#' cov_pars = RTMB::AD(c(0,0.5),force=TRUE);   #--force cov_pars to advector
#' Q = mkPrecMat(5,cov_pars,"ar1",sparse=TRUE);#--Q is adsparse
#' as.matrix(Q);
#' C = convertPrecMatToCovMat(Q);              #--C is adsparse
#' as.matrix(C);
#'
#' @export
#'
convertPrecMatToCovMat<-function(Q){solve(Q)};

#' @title Convert a covariance matrix to a precision matrix
#' @param C - covariance matrix
#' @return precision matrix
#' @details Uses [solve] to convert the covariance matrix to a precision matrix. If
#' C is a sparse matrix class, the resulting precision matrix may not be.
#' @examplesIf FALSE
#' # example code
#' ##--R context
#' Q = mkPrecMat(5,c(0,0.5),"ar1",sparse=TRUE);
#' C = convertPrecMatToCovMat(Q);
#' t = getThetaFromCovMat(C);
#' Cp = getCorrMatFromTheta(t);
#' Qp = convertCovMatToPrecMat(Cp);
#'
#' ##--RTMB context
#' cov_pars = RTMB::AD(c(0,0.5),force=TRUE);   #--forcing to advector
#' Q = mkPrecMat(5,cov_pars,"ar1",sparse=TRUE);#--Q is adsparse
#' as.matrix(Q);
#' C = convertPrecMatToCovMat(Q);              #--C is adsparse
#' as.matrix(C);
#'
#' @export
#'
convertCovMatToPrecMat<-function(C){solve(C)};





#' @title Get the `theta` (Cholesky factor) parameterization of a covariance matrix
#' @description Function to get `theta` (Cholesky factor) parameterization of
#' a covariance matrix; from Balint Tamasi, TMB user's group list.
#' @param S a covariance matrix
#' @param corrs.only return only values corresponding to the correlation matrix parameters?
#' @return the corresponding \code{theta} parameter vector
#' @details If `corrs.only`=TRUE, \code{theta} corresponds to a strictly lower triangle matrix.
#' If `corrs.only`=FALSE, then `theta` corresponds to a lower triangle matrix. If `n` is the length of the main diagonal of `S`, then
#' the first `n` elements of the returned \code{theta} vector are given by
#' \code{log(diag(S))/2}.
#'
#' @examplesIf FALSE
#' # example code
#' ##--ordinary R
#' S = matrix(c(2,0.3,0.3,2),nrow=2);
#' getThetaFromCovMat(S);
#' getCorrMatFromTheta(getThetaFromCovMat(S),strictly.lower.triangle=FALSE)
#' getThetaFromCovMat(S,corrs.only=TRUE);
#'
#' ##--with RTMB
#' require(RTMB);
#'
#' ###--NOT the way  to test RTMB context(!!)
#' ####--error arises from RTMB approach to making R Cholesky factorization RTMB-compatible
#' S = matrix(c(2,0.3,0.3,2),nrow=2)
#' S = RTMB::AD(S,force=TRUE);
#' getThetaFromCovMat(S)                  #--generates C stack error!
#'
#' ###--the following run fine
#' S = matrix(c(2,0.3,0.3,2),nrow=2)
#' tp = MakeTape(getThetaFromCovMat,S);
#' tp(S)
#' tp = MakeTape(function(x){getThetaFromCovMat(x,corrs.only=TRUE)},S);
#' tp(S)
#'
#' @importFrom stats cov2cor
#' @export
#'
getThetaFromCovMat <- function(S, corrs.only=FALSE) {
  logsd <- log(diag(S))/2;         #--vector of log(sd)'s (log of sqrt of diagonals)
  cr <- cov2cor(S);                #--convert covariance to correlation matrix
  Lu <- chol(cr);                  #--Cholesky factorization (upper triangle)
  Lu <- Lu %*% diag(1 / diag(Lu)); #--scale
  corrs <- Lu[upper.tri(Lu)];      #--correlation parameters as upper triangular elements
  if (corrs.only) return(corrs);
  ret <- c(logsd,corrs);
  return(ret)
}

#' @title Get a correlation (or covariance) matrix from a Cholesky factor (lower triangle) `theta` parameterization
#' @description Function to get a correlation (or covariance) matrix from a Cholesky factor (lower triangle) `theta` parameterization.
#' @param t - a "\code{theta} vector for a Cholesky matrix, by row
#' @param strictly.lower.triangle - if TRUE, `theta` corresponds to a strictly lower triangle Cholesky factor; if FALSE, it includes values along the main diagonal
#' @param return.cov- flag to return covariance matrix if strictly.lower.triangle=FALSE (default=FALSE)
#' @return if `strictly.lower.triangle` is TRUE (default), the corresponding correlation matrix; otherwise, a covariance matrix
#' @details If `strictly.lower.triangle` is FALSE, the first `n` elements of the input `theta` are interpreted as the elements of a `n` x `n`
#' diagonal scaling matrix, with values equal to `log(sd)/2`, where `sd^2` is the vector of variances scaling the correlation
#' matrix to a covariance matrix.
#'
#' @examplesIf FALSE
#' # example code
#' ##--ordinary R
#' tp = 0.3144855;
#' getCorrMatFromTheta(tp);
#' t = getThetaFromCovMat(getCorrMatFromTheta(tp),corrs.only=FALSE);
#' getCorrMatFromTheta(t,strictly.lower.triangle=FALSE);
#'
#' S = matrix(c(2,0.3,0.3,2),nrow=2);
#' t = getThetaFromCovMat(S,corrs.only=FALSE);
#' Sp = getCorrMatFromTheta(t,strictly.lower.triangle=FALSE,return.cov=TRUE);
#'
#' #--round trip
#' t - as.theta.vcov(getCorrMatFromTheta(t),TRUE);
#'
#' ##--with RTMB
#' require(RTMB);
#' tp = RTMB::MakeTape(getCorrMatFromTheta,numeric(1));
#' tp(t);
#' getCorrMatFromTheta(RTMB::AD(t));
#' #--can't do roundtrip as example because of `chol` call in `getThetaFromCovMat`.
#' @export
#'
getCorrMatFromTheta<-function(t,strictly.lower.triangle=TRUE,return.cov=FALSE){
  d = getMatDimFromLwrTriLen(length(t),strictly.lower.triangle=strictly.lower.triangle);
  Ll = diag(x=1,nrow=d,names=FALSE);
  if (strictly.lower.triangle) {lt = t;} else {lt = t[(d+1):length(t)];}
  Ll[lower.tri(Ll)] <- lt;
  #diag(Ll) = 1.0;
  LtL = Ll %*% t(Ll);
  s = diag(1/sqrt(diag(LtL)));
  C = s %*% (LtL) %*% t(s);    #--don't really need to transpose `s`
  if ((!strictly.lower.triangle)&&return.cov)
    C = diag(exp(t[1:d])) %*% C %*% diag(exp(t[1:d]));
  return(C);
}

#' @title Make a precision matrix for (possibly multiple) AR1 or iid processes
#' @description Function for creating a precision matrix for
#' (possibly multiple) AR1 or iid (if `rho`=0) processes
#' @usage mkPrecMat.ar1(n, sigma=1, rho=0, sdX=NULL, sparse=TRUE)
#' @param n integer vector > 0, number of elements in each AR1 process (no default)
#' @param sigma - vector with the standard deviation for each white noise process
#' @param rho - vector with the correlation coefficient between successive values for each AR1 process (default=0)
#' @param sdX - vector with standard deviations of the observed process (default=1)
#' @param sparse - (T/F) Should the matrix be of class 'dsCMatrix' (default=0)
#'
#' @return a precision matrix with block-diagonal AR1 structure and attributes "params",
#' a vector with the values for `n`, `sigma` and `rho`, and "qType", which is "diag"
#' if `rho` is 0 and "ar1" otherwise.
#'
#' @details A precision matrix for multiple independent AR1 and/or iid processes can be
#' created by making `n`, `sigma` (or `sdX`), and `rho` vectors. `sigma` (or `sdX`), and `rho`
#' are expanded (or truncated) to the length of `n`.
#'
#' Note that the variance associated with individual observations is
#'
#' $$ var(X_i) = E(X^2_i)-\mu^2 = \frac{\sigma^2}{1-\rho^2} $$
#'
#' and thus
#'
#' $$ \sigma = \sqrt{{1-\rho^2}} \cdot sd(X) $$
#'
#' and consequently the precision matrix can (optionally) be specified in terms of the standard
#' deviation of the observed process (`sdX`).
#'
#' @examplesIf FALSE
#' # example code
#' ##--R context
#' Q <- mkPrecMat.ar1(n=5,sigma=2,rho=0.5,sparse=TRUE);
#' Q <- mkPrecMat.ar1(n=c(5,5),sigma=2,rho=0.5,sparse=TRUE);
#' C <- solve(Q); #--covariance matrix
#'
#' Qp <- mkPrecMat.ar1(n=5,sigma=0,rho=0.5, sdX=2/sqrt(1-0.5^2), sparse=TRUE);
#'
#'
#' ##--RTMB context
#' require(RTMB);
#' Q <- mkPrecMat.ar1(n=5,sigma=2, rho=0.5, sdX=0, sparse=TRUE, force=TRUE);
#' C <- solve(Q); #--covariance matrix
#' Qp <- mkPrecMat.ar1(n=5, sigma=0, rho=0.5, sdX=2/sqrt(1-0.5^2), sparse=TRUE, force=TRUE);
#' Q-Qp
#' Qp <- mkPrecMat.ar1(n=5, sigma=0, rho=RTMB::AD(0.5,force=TRUE), sdX=RTMB::AD(2/sqrt(1-0.5^2),force=TRUE), sparse=TRUE, force=FALSE);
#' @importFrom Matrix Diagonal Matrix
#' @importFrom RTMB AD
#' @md
#' @export
#'
mkPrecMat.ar1 <- function(n, sigma=1, rho=0, sdX=0, sparse=TRUE, force=FALSE){
  #--testing: n=5;sigma=2;rho=0.5;sparse=TRUE;
  #if (abs(rho)>=1) warning("rho >= 1. not a valid AR1 process.")
  if (any(sapply(list(sigma,rho,sdX),class)=="advector")) force=TRUE;
  nQ = length(n);
  vn = n;
  vrho    = RTMB::AD(rep_len(rho,length.out=nQ),force=force);
  vsigma  = RTMB::AD(rep_len(sigma,length.out=nQ),force=force);
  vsdX    = RTMB::AD(rep_len(sdX,length.out=nQ),force=force);
  vsigma  = vsigma + sqrt((1-vrho^2))*vsdX;#--either sigma or sdX should be non-zero, not both
  if (!force) cat("vsigma:",vsigma,"\n");
  #if(sigma <= 0) stop("sigma parameter must be greater than 0.");
#  d <- c(RTMB::AD(1.0,force=force),rep(1.0+rho^2,(n-2)),RTMB::AD(1.0,force=force));
  lstQ = list();
  for (i in 1:nQ){
    n     = vn[i];
    rho   = vrho[i];
    sigma = vsigma[i];
    d <- RTMB::AD(rep(1,n),force=force);
    d[2:(n-1)] <- rep(1.0+rho^2,(n-2));
    u <- rep(-1*rho,(n-1));
    Q <- RTMB::AD(Matrix::Matrix(0,n,n),force=force);
    Q[getRowColIndices_band(n,c(-1,0,1))] = c(u,d,u);
    D <- RTMB::AD(Matrix::Matrix(0,n,n),force=force);
  #  D[getRowColIndices_diag(n,0)] = RTMB::AD(1.0,force=force)/sigma;
    D[getRowColIndices_diag(n,0)] = 1.0/sigma;
    Q <- D %*% Q %*% D;
    if(!sparse) Q <- RTMB::AD(Matrix::Matrix(Q, sparse=FALSE),force=force);
    lstQ[[i]] = Q;
  }
  Qt = buildBlockDiagonalMat(lstQ,force=force);
  attr(Qt,"qType")  = "ar1";
  attr(Qt,"params") = c(n=vn,sigma=vsigma,rho=vrho);
  return(Qt);
}

#'
#' @title Make a precision matrix
#' @description Function to make a precision matrix.
#' @param n - dimension of the matrix
#' @param pars - vector of parameter values determining the precision matrix
#' @param type - type of precision matrix
#' @param sparse - flag to create a sparse or dense matrix
#' @param force - force an RTMB adsparse representation
#' @return a precision matrix (a Matrix type or RTMB "adsparse" type)
#' @details
#' If `force` is FALSE, the returned matrix will be an "adsparse" matrix only if
#' `pars` is an "advector"; otherwise, it will be a [Matrix] type.
#'
#' The `type`s of precision matrices that can be created using this function are currently:
#' \itemize{
#'  \item{"diag" - diagonal iid matrix with variance `V` along the diagonal}
#'  \item{"iid" - synonym for "diag"}
#'  \item{"ar1" - AR1 precision matrix with variance `V` and correlation coefficient rho (latter parameterized as rho=symlogistic(par))}
#' }
#'
#' Variances `V` are generally parameterized on the ln scale such that \eqn{V = exp(2*par)}.
#'
#' @export
#'
mkPrecMat<-function(n,pars,type,sparse=TRUE,force=FALSE){
  if (!RTMB:::ad_context()) print(pars);
  Q = NULL;
  if (type %in% c("diag","iid")) {
    Q = mkPrecMat.ar1(n,sigma=exp(2*pars[1]),sparse=sparse,force=force);
  } else
  if (type %in% c("ar1")) {
    Q = mkPrecMat.ar1(n,sigma=exp(2*pars[1]),rho=symlogistic(pars[2]),sparse=sparse,force=force);
  }
  return(Q);
}

