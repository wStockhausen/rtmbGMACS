#--misc. covariance/correlation matrix-related functions

#--RE variance-covariance related code
##--original code adapted from reformulas, glmmTMB (enum.R, glmmTMB.R)

#--from glmmTMB/R/glmmTMB.R: get_mat_dim
#'@title From length of (strictly) lower triangle, compute dimension of matrix
#'@param ntri - length of vector determining the (strictly) lower triangle
#'@param strictly.lower.triangle - flag (TRUE/FALSE)
#'@return the row or column dimension of the related square matrix
#'@details If strictly.lower.triangle=TRUE, the triangular matrix dimension
#'is strictly lower (i.e., below the main diagonal), if FALSE, it includes the
#'main diagonal.
#'@examplesIf FALSE
#'getMatDimFromLwrTriLen(getLwrTriLenFromMatDim(8));
#'getMatDimFromLwrTriLen(getLwrTriLenFromMatDim(8,strictly.lower.triangle=FALSE),strictly.lower.triangle=FALSE);
#'@export
getMatDimFromLwrTriLen <- function(ntri,strictly.lower.triangle=TRUE) {
    n = as.integer(round(0.5 * (1 + sqrt(1 + 8 * ntri))));
    if (!strictly.lower.triangle) n = n-1;
    return(n);
}

#'@title From dimension of square matrix, compute length of strictly lower triangle
#'@param matdim - dimension (number of rows or columns) of square matrix
#'@param strictly.lower.triangle - flag (TRUE/FALSE; default=TRUE)
#'@return the length of a vector defining the (strictly, if `strictly.lower.triangle`=FALSE) lower triangle.
#'@examplesIf FALSE
#'t1 = getLwrTriLenFromMatDim(8);
#'t2 = getLwrTriLenFromMatDim(8,strictly.lower.triangle=FALSE);
#'@export
getLwrTriLenFromMatDim <- function(matdim,strictly.lower.triangle=TRUE) {
  n = sum(1:(matdim-1)) + (!strictly.lower.triangle)*matdim;
  return(n);
}

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
#' S = matrix(c(2,0.3,0.3,2),nrow=2)
#' S = RTMB::AD(S,force=TRUE);
#' getThetaFromCovMat(S)          #--generates C stack error!
#' ###--the following runs fine
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

#' @title Get a matrix of row/column indices along the diagonal of an n x n matrix.
#' @param n -size of square matrix
#' @param diag - integer offset from main diagonal
#' @return integer matrix of row/column indices
#' @export
#'
getRowColIndices_diag<-function(n,diag){
  if (diag>=0) return(cbind(1:(n-diag), (diag+1):n));
  cbind((1-diag):n, 1:(n+diag));
}

#' @title Get a matrix of row/column indices along the diagonals of an n x n matrix.
#' @param n -size of square matrix
#' @param diag - vector of integer offsets from the main diagonal
#' @return integer matrix of row/column indices
#' @export
#'
getRowColIndices_band<-function(n,diags){
  m = getRowColIndices_diag(n,diags[1]);
  i = 2;
  while (i<=length(diags)){
    m = rbind(m,getRowColIndices_diag(n,diags[i]));
    i = i+1;
  }
  return(m);
}

#' @title Convert a precision matrix to a covariance matrix
#' @param Q - precision matrix
#' @return covariance matrix
#' @details Uses [solve] to convert the precision matrix to a covariance matrix. If
#' Q is a sparse matrix class, the resulting covariance matrix may not be.
#' @examplesIf FALSE
#' # example code
#' ##--R context
#' Q = mkPrecMat(5,c(0,0.5),"ar1",sparse=TRUE,force=FALSE);
#' M = convertPrecMatToCovMat(Q);
#' t = getThetaFromCovMat(M);
#' Mp = getCorrMatFromTheta(t);
#' Qp = convertCovMatToPrecMat(Mp);
#'
#' ##--RTMB context
#' Q = mkPrecMat(5,RTMB::AD(c(0,0.5),force=TRUE),"ar1",sparse=TRUE,force=FALSE);
#' as.matrix(Q);
#' C = convertPrecMatToCovMat(Q);
#' as.matrix(C);
#'
#' @export
#'
convertPrecMatToCovMat<-function(Q){solve(Q)};

convertCovMatToPrecMat<-function(M){solve(M)};

#' @title Get the starting rows/columns of blocks for a list of matrices to be sub-matrices in a "super" block-diagonal matrix
#' @description Function to get the starting rows/columns of blocks for a list of matrices to be sub-matrices in a "super" block-diagonal matrix.
#' @param lst - a list of block-diagonal matrices
#' @return a vector of starting rows/columns of blocks in the "super" block-diagonal matrix
#' @export
#'
getBlockDiagStarts<-function(lst){
  nrw = sapply(lst,nrow);
  bdi = lapply(lst,getBlockDiagInfo);
  if (sum(bdi[[1]])>0) {rws = bdi[[1]];} else {rws = 0;}
  j = length(rws); trw = nrw[1]; i=2;
  while (i <= length(lst)){
    rws = c(rws,trw+bdi[[i]]);
    trw = trw+nrw[i-1];
    i = i+1;
  }
  rws
}

#' @title Get the starting rows/columns of blocks in a block-diagonal matrix
#' @description Function to get the starting rows/columns of blocks in a block-diagonal matrix.
#' @param lst - a matrix
#' @return a vector of starting rows/columns of blocks in the block-diagonal matrix
#' @details If the matrix is not block-diagonal (i.e., it does not yet have a "block_start" attribute),
#' then its "block" starts at the first row (using a 0-based index). If the matrix is block-diagonal, its
#' "block_start" attribute (a 0-based vector of rows at which blocks start) is returned.
#' @export
#'
getBlockDiagInfo<-function(x){
  bdi = attr(x,"block_start");
  if (is.null(bdi)) bdi = 0;
  return(bdi);
}

#' @title Build a block-diagonal matrix from building block matrices
#' @description Build a block-diagonal matrix from building block matrices.
#' @param ... - individual matrices or a list of matrices
#' @return block-diagonal matrix, either an "adsparse" matrix or inheriting from [Matrix] class [Matrix::CsparseMatrix].
#' @importFrom Matrix .bdiag .M2C Matrix
#' @importFrom RTMB AD
#' @export
#'
buildBlockDiagonalMat<-function(...,force=FALSE){
  #--start processing
  if (RTMB:::ad_context()||force){
    #--RTMB context: can't use Matrix::.M2C or Matrix::.bdiag
      if ((n <- ...length()) == 0L) {
      #--empty input
      M = RTMB::AD(new("dgCMatrix"),force=force);
    } else {
      if (n > 1L) {
        lst = list(...); #--list of block diagonal matrices
      } else {
        lst <- ..1;
        if (!is.list(lst)) lst = list(lst); #--`lst` is a single matrix, wrap it in a list
      }
      #--`lst` should be a list of matrices
      trw = sum(sapply(lst,nrow));
      rws = getBlockDiagStarts(lst); #--zero-based
      M = RTMB::AD(Matrix::Matrix(0,trw,trw),force=force);
      Mx = M@x;#--0-length (numeric or ad)vector
      Mi = M@i;#--0-length integer vector
      rowstart=0;
      diffp = integer();#--0-length integer vector from which to compute M@p
      for (l in 1:length(lst)){
        #--l=1;
        Mp = RTMB::AD(lst[[l]],force=force);#--lth (sub) matrix
        Mpx = RTMB::AD(Mp@x,force=force);
        Mpi = Mp@i; #--row indices
        Mpp = Mp@p;
        Mx = c(Mx,Mpx);
        Mi = c(Mi,rowstart+Mpi);
        rowstart = rowstart + nrow(Mp);
        diffp=c(diffp,diff(Mpp));
      }
      M@x = Mx;
      M@i = as.integer(Mi);
      M@p = as.integer(c(0,cumsum(diffp)));
      attr(M,"block_start") = rws;
    }
  } else {
    #--Matrix context
    M = NULL;
    if ((n <- ...length()) == 0L) {
      #--empty input
      M = new("dgCMatrix");
    } else if (n > 1L) {
      lst = list(...);
      rws = getBlockDiagStarts(lst);
      M = Matrix::.M2C(Matrix::.bdiag(lst));
      attr(M,"block_start") = rws;
    } else if (!is.list(x <- ..1)) {  # `x` is now the input object
      #--`x` is a single matrix
      bdi = getBlockDiagInfo(x);
      M = as(x, "CsparseMatrix");
      attr(M,"block_start") = bdi;
    } else if (length(x) == 1L) {
      #--`x` is a single matrix wrapped in a list
      bdi = getBlockDiagInfo(x[[1L]]);
      M = as(x[[1L]], "CsparseMatrix");
      attr(M,"block_start") = bdi;
    } else {
      #--`x` is a list of matrices
      rws = getBlockDiagStarts(x);
      M = Matrix::.M2C(Matrix::.bdiag(x));
      attr(M,"block_start") = rws;
    }
  }
  return(M);
}
