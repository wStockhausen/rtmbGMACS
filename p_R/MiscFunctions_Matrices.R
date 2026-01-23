#--misc. covariance/correlation matrix-related functions

#--RE variance-covariance related code
##--original code adapted from reformulas, glmmTMB (enum.R, glmmTMB.R)

#--from glmmTMB/R/glmmTMB.R
#'@title From length of (strictly) lower triangle, compute dimension of matrix
#'
get_matdim <- function(ntri) {
    as.integer(round(0.5 * (1 + sqrt(1 + 8 * ntri))))
}

#--from glmmTMB/R/glmmTMB.R
##' transform correlation parameters to and from glmmTMB parameterization
##' @param theta vector of internal correlation parameters (elements of scaled Cholesky factor, in \emph{row-major} order)
##' @param return_val return a vector of correlation values from the lower triangle ("vec"), or the full correlation matrix ("mat")?
##' @return a vector of correlation values (\code{get_cor}) or glmmTMB scaled-correlation parameters (\code{put_cor})
##' @details
##' \code{\link{get_cor}} transforms from the glmmTMB parameterization (components of a \code{theta} parameter vector) to correlations;
##' \code{\link{put_cor}} does the inverse transformations, from correlations to \code{theta} values.
##'
##' These functions follow the definition at \url{http://kaskr.github.io/adcomp/classdensity_1_1UNSTRUCTURED__CORR__t.html}:
##' if \eqn{L} is the lower-triangular matrix with 1 on the diagonal and the correlation parameters in the lower triangle,
##' then the correlation matrix is defined as \eqn{\Sigma = D^{-1/2} L L^\top D^{-1/2}}{Sigma = sqrt(D) L L' sqrt(D)},
##' where \eqn{D = \textrm{diag}(L L^\top)}{D = diag(L L')}. For a single correlation parameter \eqn{\theta_0}{theta0}
##' (i.e. the correlation in a 2x2 correlation matrix), this works out to \eqn{\rho = \theta_0/\sqrt{1+\theta_0^2}}{rho = theta0/sqrt(1+theta0^2)}.
##' The \code{get_cor} function returns the elements of the lower triangle of the correlation matrix, in column-major order.
##'
##' These functions also work for AR1 correlation parameters.
##' @examplesIf FALSE
##' th0 <- 0.5
##' stopifnot(all.equal(get_cor(th0), th0/sqrt(1+th0^2)))
##' set.seed(101)
##' ## pick 6 values for a random 4x4 correlation matrix
##' print(C <- get_cor(rnorm(6), return_val = "mat"), digits = 3)
##' ## transform a correlation matrix to a theta vector
##' cor_mat <- matrix(c(1,0.3,0.1,
##'                     0.3,1,0.2,
##'                     0.1,0.2,1), ncol = 3)
##' put_cor(cor_mat, "mat")
##' put_cor(cor_mat[lower.tri(cor_mat)], "vec")
##' ## test: round-trip
##' stopifnot(all.equal(get_cor(put_cor(C), return_val = "mat"), C))
##' @export
get_cor <- function(theta, return_val = c("vec", "mat")) {
    return_val <- match.arg(return_val)
    n <-  get_matdim(length(theta))
    R <- diag(n)
    R[upper.tri(R)] <- theta
    R[] <- crossprod(R) # R <- t(R) %*% R
    scale <- 1 / sqrt(diag(R))
    R[] <- scale * R * rep(scale, each = n) # R <- cov2cor(R)
    if (return_val == "mat") return(R)
    return(R[lower.tri(R)])
}

#--from glmmTMB/R/glmmTMB.R
##' @rdname get_cor
##' @param C a correlation matrix
##' @param input_val input a vector of correlation values from the lower triangle ("vec"), or the full correlation matrix ("mat")?
#' @return the corresponding \code{theta} parameter vector without logsd scale parameter
##' @details This function is similar to [as.theta.vcov()], except this function
##' @export
put_cor <- function(C, input_val = c("mat", "vec")) {
    input_val <- match.arg(input_val);
    if (input_val == "vec") {
        ## construct matrix from input vector
        M <- diag(get_matdim(length(C)));
        M[lower.tri(M)] <- C;
        M[upper.tri(M)] <- t(M)[upper.tri(M)];
        C <- M;
    }
    cc2 <- chol(C);
    scale <- diag(cc2);
    cc2 <- cc2 %*% diag(1/scale);
    cc2[upper.tri(cc2)];
}

##' Get theta parameterisation of a covariance structure
## from Balint Tamasi, TMB user's group list
## FIXME: names based on dimnames of Sigma?
#' @param Sigma a covariance matrix
#' @param corrs.only return only values corresponding to the correlation matrix parameters?
#' @return the corresponding \code{theta} parameter vector
#' @importFrom stats cov2cor
as.theta.vcov <- function(Sigma, corrs.only=FALSE) {
  logsd <- log(diag(Sigma))/2; #--log(sd)
  cr <- cov2cor(Sigma);        #--convert covariance to correlation matrix
  cc <- chol(cr)               #--Cholesky factorization
  cc <- cc %*% diag(1 / diag(cc)) #--scale
  corrs <- cc[upper.tri(cc)]      #--correlation parameters as upper triangular elements
  if (corrs.only) return(corrs)
  ret <- c(logsd,corrs)
  return(ret)
}

#' @title Make a precision matrix for an AR1 process
#' @description Function for creating a precision matrix for
#' an AR1 process
#' @usage mkPrecMat.ar1(M, sigma=NULL, rho=0, sdX=NULL, sparse=TRUE)
#' @param n int > 0, number of elements in the AR1 process (default=NULL)
#' @param sigma - the standard deviation of the white noise process
#' @param rho - the correlation coefficient between successive values (default=0)
#' @param sdX - standard deviation of the observed process (default=NULL)
#' @param sparse - (T/F) Should the matrix be of class 'dsCMatrix' (default=TRUE)
#'
#' @return a precision matrix with an AR1 structure and attributes "params",
#' a vector with the values for `sigma` and `rho`, and "qType", which is "diag"
#' if `rho` is 0 and "ar1" otherwise
#'
#' @details Note that the variance associated with individual observations is
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
#' Q <- mkPrecMat.ar1(n=5,sigma=2,rho=0.5,sparse=TRUE);
#' C <- solve(Q); #--covariance matrix
#' Qp <- mkPrecMat.ar1(n=5,rho=0.5, sdX=2/sqrt(1-0.5^2), sparse=TRUE);
#' Q-Qp
#'
#' @importFrom Matrix bandSparse Diagonal Matrix
#' @md
#' @export
#'
mkPrecMat.ar1 <- function(n, sigma=1, rho=0, sdX=0, sparse=TRUE, force=FALSE){
  #--testing: n=5;sigma=2;rho=0.5;sparse=TRUE;
  #if (abs(rho)>=1) warning("rho >= 1. not a valid AR1 process.")
  rho = RTMB::AD(rho,force=force);
  sigma = RTMB::AD(sigma,force=force) + sqrt((1-rho^2))*RTMB::AD(sdX,force=force);#--either sigma or sdX should be non-zero, not both
  #if(sigma <= 0) stop("sigma parameter must be greater than 0.");
  d <- c(RTMB::AD(1.0,force=force),rep(1.0+rho^2,(n-2)),RTMB::AD(1.0,force=force));
  u <- rep(-1*rho,(n-1));
  # Q = Matrix::bandSparse(n, k=c(-1,0,1),diagonals=list(u,d,u),repr="C",symmetric=FALSE);
  Q <- RTMB::AD(Matrix::Matrix(0,n,n),force=force);
  Q[getRowColIndices_band(n,c(-1,0,1))] = c(u,d,u);
  D <- RTMB::AD(Matrix::Matrix(0,n,n),force=force);
  D[getRowColIndices_diag(n,0)] = RTMB::AD(1.0,force=force)/sigma;
  Q <- D %*% Q %*% D;
  if(!sparse) Q <- RTMB::AD(Matrix::Matrix(Q, sparse=FALSE),force=force);
  attr(Q,"qType")  = "ar1";
  attr(Q,"params") = c(sigma=sigma,rho=rho);
  return(Q);
}

mkPrecMat<-function(n,pars,type,sparse=TRUE,force=FALSE){
  if (!RTMB:::ad_context()) print(pars);
  Q = NULL;
  if (type %in% c("diag")) {
    Q = mkPrecMat.ar1(n,sigma=exp(pars[1]),sparse=sparse,force=force);
  } else
  if (type %in% c("ar1")) {
    Q = mkPrecMat.ar1(n,sigma=exp(pars[1]),rho=syminvlogit(pars[2]),sparse=sparse,force=force);
  }
  return(Q);
}

#' @title Get a matrix of row/column indices along the diagonal of an n x n matrix.
#' @param n -size of square matrix
#' @param diag - integer offset from main diagonal
#' @return integer matrix of row/column indices
#' @export
getRowColIndices_diag<-function(n,diag){
  if (diag>=0) return(cbind(1:(n-diag), (diag+1):n));
  cbind((1-diag):n, 1:(n+diag));
}

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
#' @export
convertPrecMatToCovMat<-function(Q){solve(Q)};

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
    # if ((n <- ...length()) == 0L) {
    #   #--empty input
    #   M = RTMB::AD(new("dgCMatrix"),force=force);
    # } else {
    #   if (n > 1L) {
    #     lst = list(...); #--list of block diagonal matrices
    #   } else {
    #     lst <- ..1;
    #     if (!is.list(lst)) lst = list(lst); #--`lst` is a single matrix, wrap it in a list
    #   }
    #   #--`lst` should be a list of matrices
    #   trw = sum(sapply(lst,nrow));
    #   rws = getBlockDiagStarts(lst); #--zero-based
    #   M = RTMB::AD(Matrix::Matrix(0,trw,trw),force=force);
    #   start=0;
    #   for (l in 1:length(lst)){
    #     Mp = RTMB::AD(lst[[l]],force=force);#--lth (sub) matrix
    #     n = nrow(Mp);
    #     i = (1:n);
    #     cs = as.matrix(tidyr::expand_grid(i,i)); #--2-column matrix of locations in M corresponding to locations in Mp
    #     M[cs+start] = Mp[cs];
    #     start = start + n;
    #   }
    #   attr(M,"block_start") = rws;
    # }
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
