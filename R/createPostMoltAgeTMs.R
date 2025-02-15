#--create post-molt age transition matrices
createPostMoltAgeTMs<-function(dmsC){
  nCs = nrow(dmsC);
  tmPMAI = matrix(0,nCs,nCs);#--matrix to increment post-molt age by one age category
  tmPMAD = matrix(0,nCs,nCs);#--matrix to decrement post-molt age to first age category
  ps = attr(dmsC,"dmlvs")$p
  np = length(ps);

  #  given row index irC pointing to rowC and p=rp,
  #  want to find column index icC such that colC has p=pc=min(rp+1,max(ps)) and otherwise rowC = colC
  ## increment rp to desired pc in rowC, then do join to colC such that rowC r,x,m,p,z = colC r,x,p,z
  ## resulting sparse_idx's are row and column indices for 1 in transition matrix.
  pC = ps[min(as.integer(dmsC$p)
  for (ir in 1:nCs){ #--loop over rows
    #--for testing: ir = 1;
    rowC = dmsC[ir,] |> dplyr::rename(ir=sparse_idx);
    cp = ps[min(as.integer(rowC$p)+1,np];
    colC =
    for (ic in 1:nCs) {#--loop over columns
      #--for testing: ic = 1;
      cp = dmsC[ic,]$p;
    }#--ic loop
  }#--ir loop

  return(list(tmPMAI=tmPMAI,tmPMAD=tmPMAD));
}
