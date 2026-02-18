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
#' @return block-diagonal matrix, either an "adsparse" matrix or inheriting from [Matrix] class `CsparseMatrix`.
#' @importFrom Matrix .bdiag .M2C Matrix
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
