#--create post-molt age transition matrices (TODO: evaluate speed of assignment methods)
#' @title Create transition matrices between post-molt age categories
#' @description Function to create two transition matrices between post-molt age categories
#' @param dmsC - the model dimensions dataframe for population classes (i.e., dims$dmsC when setting up model)
#' @return list with two transition matrices (`tmPMAD` and `tmPMAI`) between post-molt age categories (see Details below)
#' @details This function creates transition matrices for crab classes that will (`tmPMAD`) or will not
#' (`tmPMAI`) undergo molting during a growth season. These matrices contain only constants (0's or 1's)
#' and thus can be created outside the objective function and should be included as part of the `inputs` list.
#'
#' `tmPMAI` increments the post-molt age category for non-molting crab by one category (e.g., "new-shell" to "old_shell").
#' If nN is a vector of the numbers of crabs, by model population class, at the start of a growth season that will not molt,
#' then
#'
#' $$nN^+ = tmPMAI %*% nN$$
#'
#' "moves" crab from post-molt age class `p` to `p+1` (or leaves them in the last class, which acts as an aggregator).
#'
#' `tmPMAD` decrements the post-molt age category for molting crab to the first category (i.e., post-molt age 0-e.g., "new_shell").
#' If mN is a vector of the numbers of crabs, by model population class, at the start of a growth season that will molt,
#' then
#'
#' $$mN^+ = tmPMAD %*% mN$$
#'
#' "moves" crab from post-molt age class `p` to the first class (post-molt age 0--e.g., "new_shell").
#'
#' @import dplyr
#'
#' @md
#' @export
#'
createPostMoltAgeTMs<-function(dmsC){
  nCs = nrow(dmsC);
  ps = attr(dmsC,"dmlvs")$p
  np = length(ps);

  ##--column is where you start, row is where you end up
  colCs = dmsC |> dplyr::rename(col_=sparse_idx,p_col=p);
  rowCs = dmsC |> dplyr::rename(row_=sparse_idx,p_row=p);

  #--create matrix that moves crab that WILL NOT MOLT to next post-molt age----
  ##  given col index icC pointing to colC (reflecting initial categories) and p=cp,
  ##  want to find row index irC such that rowC (final categories) has p=rp=min(cp+1,max(ps)) and otherwise colC = rowC
  ### increment cp to desired pc in colC, then do join to rowC such that colC r,x,m,p,z = rowC r,x,p,z
  ### resulting sparse_idx's are column and row indices for 1 in transition matrix.
  ##--match up cp -> rp
  crsC = colCs |>
           dplyr::inner_join(rowCs,by=dplyr::join_by(r, x, m, z),relationship="many-to-many") |>
           dplyr::filter(p_row==ps[min(as.integer(p_col)+1,np)]);
  ##--TODO: choose faster of the following 2 methods
  ##--assign by indexing using for loop
  tmPMAI = matrix(0,nCs,nCs);#--matrix to increment post-molt age by one age category
  for (i in 1:nCs)
    tmPMAI[crsC$row_[i],crsC$col_[i]] = 1;
  View(tmPMAI);
  ##--assign using an index array (probably faster than above?)
  tmPMAI = matrix(0,nCs,nCs);#--matrix to increment post-molt age by one age category
  tmPMAI[array(c(crsC$row_,crsC$col_),dim=c(nCs,2))] = 1;
  View(tmPMAI);

  #--create matrix that moves crab that WILL MOLT to first post-molt stage (i.e., post-molt age 0)----
  ##  given col index icC pointing to colC (reflecting initial categories) and p=cp,
  ##  want to find row index irC such that rowC (final categories) has p=rp=ps[1] and otherwise colC = rowC
  crsC = colCs |>
           dplyr::inner_join(rowCs,by=dplyr::join_by(r, x, m, z),relationship="many-to-many") |>
           dplyr::filter(p_row==ps[1]);
  ##--TODO: choose faster of the following 2 methods
  ##--assign by indexing using for loop
  tmPMAD = matrix(0,nCs,nCs);#--matrix to decrement post-molt age to first age category
  for (i in 1:nCs)
    tmPMAD[crsC$row_[i],crsC$col_[i]] = 1;
  View(tmPMAD);
  ##--assign using an index array (probably faster than above?)
  tmPMAD = matrix(0,nCs,nCs);#--matrix to increment post-molt age by one age category
  tmPMAD[array(c(crsC$row_,crsC$col_),dim=c(nCs,2))] = 1;
  View(tmPMAD);

  return(list(tmPMAI=tmPMAI,tmPMAD=tmPMAD));
}
