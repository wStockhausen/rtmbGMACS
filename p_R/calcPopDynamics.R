#'
#' @title Run population dynamics
#' @description Run population dynamics for model time period.
#'
#' @param dims - model dimensions list
#' @param options - model options list
#' @param params - model parameters list
#' @param verbose - flag to print diagnostic info
#'
#' @return list with results from running the population dynamics model
#'
#' @details At each season `s`, the population in year `y`is advanced using
#'
#'   np = T %*% exp(-Z*dt) %*% G %*% n + R if `dt` > 0 otherwise
#'
#'   np = T %*% exp(-F)    %*% G %*% n + R if `dt` = 0
#'
#' where `dt` is season length (as a fraction of a year),
#' `n` is the population abundance vector at the start of the season,
#' `np` is the population abundance vector at the end of the season,
#' `T` is a matrix representing movement among regions,
#' `Z` = `M` + `F` is a diagonal matrix representing total mortality *rates*,
#' where `M` represents natural mortality rates and `F` represents
#' total fishing mortality *rates* in the case where the season is *not* instantaneous,
#' `G` is a matrix representing growth, which encompasses the probability of molting,
#' the probability of maturing, and actual growth (conditioned on having molted),
#' and R is recruitment. Movement, growth, and recruitment are treated as instantaneous
#' processes occurring at the beginning of a season, after which mortality due to
#' natural processes and/or fishing occurs (depending on season length), and then
#' recruitment occurs at the end of the season.
#'
#' Fishery capture processes are modeled by individual fleets as
#'
#' a_f = c_f/(M+F) %*% (1-exp(-Z*dt)) %*% n if `dt` >0
#'
#' a_f = c_f/F %*% (1-exp(-F)) %*% n if `dt` = 0
#'
#' r_f = rho_f * a_f;
#'
#' d_f = (1-rho_f) * af
#'
#' m_f = r_f + h_f * d_f
#'
#' where `a_f` is the abundance captured by fleet f,
#' `c_f` is the capture rate of fleet f,
#' `r_f` is the abundance retained by fleet f,
#' `rho_f` is the fraction retained by flet f,
#' `d_f` is the discard abundance (bycatch) by fleet f,
#' `m_f` is the abundance killed (fishing mortality) by fleet f,
#' and `h_f` is the handling mortality on discards in fleet f
#'
#' "Survey" capture processes are modeled as
#'
#' c_v = s_v %*% n or
#'
#' c_v = s_v %*% np
#'
#' depending on whether the "survey" occurs at the beginning or end of the season,
#' {TODO: how to handle CPUE indices from a fishery fleet when season is not instantaneous?}
#' where `c_v` is the abundance captured in the survey (scaled to the population)
#' and `s_v` is the catchability of the survey (which includes fully-selected catchability,
#' size-specific selectivity, and availability).
#'
#'
#' @export
#'
# calcPopDynamics<-function(dims,options,params,verbose){
#   nAtZ[1];
#   #--loop over years
#   for (y in years){
#     #--loop over seasons w/in a year
#     for (s in lstSsnInfo[[y]]){
#       dt = s$dt;                #--time step
#       #--all dynamics except movement occur within a region
#       ##--loop over regions for non-movement dynamics
#       for (r in lstRegInfo[[y]][[s]]){
#         #--all dynamics within a region is sex-specific
#         ##--loop over sexes
#         for (x in lstSexInfo[[y]][[s]][[r]]){
#           n_ysrx_mpz = n_ysrxmpz[y,s,r,x,,,]; #--initial abundance
#           lstSrvInfo = getSurveyInfo(y,s,r,x);
#           if (!is.null(lstSrvInfo)&&lstSrvInfo$start){
#             ##--determine survey abundance at start of time step
#             # for (v in lstSrvInfo){
#             #   if (v$start){
#             #     vQ_vysrx_mpz = getSurveyCatchability(v,y,r,x); #--vector (diagonal matrix) of size-specific catchability
#             #     n_vysrx_mpz[v,y,s] = mQ_vysrx_mpz * n_ysrx_mpz;#--survey abundance
#             #   }
#             # }#--survey loop
#           }
#
#           #--determine growth and fishery dynamics
#           lstGrwInfo  = getGrowthInfo(y,s,r,x);
#           lstFshInfo  = getFshInfo(y,s,r,x);
#           ##--instantaneous season (dt==0): no natural mortality
#           if (dt == 0){
#             if (!is.null(lstGrwInfo)){
#               ###--growth
#               n_mpz = doGrowth(n_mpz,lstGrwInfo);
#             }
#             if (!is.null(lstFshInfo)&&lstFshInfo$instantaneous){
#               ###--instantaneous fisheries catch and mortality
#               lstFshCat = doInstFsh(n_mpz,lstFshInfo)
#             }
#           } else {
#             lstNatMortInfo = getNatMortInfo(y,s,r,x);
#             if (!is.null(lstFshInfo)&&(!lstFshInfo$instantaneous)){
#               ###--instantaneous fisheries catch and mortality
#               lstFshCat = doFsh(n_mpz,lstFshInfo,lstNatMortInfo);
#               n_mpz = lstFshCat$n_mpz;#--need to do something else
#             } else {
#               n_mpz = doNatMort(n_mpz,lstNatMortInfo);
#             }
#           }
#
#           ##--determine survey abundance at end of time step
#           lstSrvInfo = getSurveysInfo(y,s,r,x,end=TRUE);
#           lstSrvCat = doSurveys
#
#         } #--sex loop
#
#       #--do movement among regions
#
#       #--do recruitment
#
#     }#--season loop
#   }#--year loop
#
# }

#' @title Apply growth to population represented by n_mpz for y,s,r,x
#' @param n_mpz - vector of abundance classified by m (maturity state), p (post-molt age), z (size)
#' @param lst - growth info for crab in y,s,r,x (i.e., lstGrw[[y]][[s]][[r]][[x]])
#' @return list with element vector n_mpz after growth (molt/growth/maturation)
#' @details `lst` is a list in which the following elements may (or have to be) present:
#'
#' \itemize{
#'   \item{maturity - string indicating whether molt to maturity depends on "pre-molt size" or "post-molt size" [required]}
#'   \item{mPtoP1 - dgCMatrix object to increment post-molt age by 1 for non-molting crab [required]}
#'   \item{mPto0 - dgCMatrix object to decrement post-molt age to 0 for molting crab [required]}
#'   \item{mImmToMat - dgCMatrix object to increment maturity state from immature to mature [required]}
#'   \item{mT - dgCMatrix object representing combined molting/growth transitions [optional]}
#'   \item{mPrGr - dgCMatrix object representing size transition matrix for molting crab when maturity is determined by post-molt size}
#'   \item{mPrGrImm - dgCMatrix object representing size transition matrix for crab not undergoing the molt to maturity}
#'   \item{mPrGrMat - dgCMatrix object representing size transition matrix for crab undergoing the molt to maturity}
#'   \item{vPrMat - vector representing probability of molting to maturity}
#'   \item{mPto0 - dgCMatrix object representing combined molting/growth transitions [optional]}
#' }
#'
#'  The value of `maturity` determines the order in which the probability of maturing is applied, either before growth occurs
#'  ("pre-molt size") or after growth occurs ("post-molt size"). If `mT` is present, then it is used to determine molting and growth under the assumption that the diagonal elements
#'  represent the probability of not molting (and thus not growing or maturing) while the remaining upper triangle elements
#'  represent the combined probability of molting and growing. If `mT` is not present, then `mPrGr` (if `maturity`=="post-molt size')
#'  or both `mPrGrImm `and `mPrGrMat` (if `maturity`=="pre-molt size") must be present (the latter allows growth to differ between
#'  crab that grow without maturing and those that undergo the molt to maturity).
#'
#' @export
#'
doGrowth<-function(n_mpz,lst){
  ret = list();
  if (!is.null(lst)){
    #--do growth
    if (!is.null(lst$mT)) {
      ##--growth based on transition matrix lst$mT (molting, growth, and aging combined)
      ###--mT is an upper triangular matrix (no negative growth)
      ####--diag(mT) is a vector reflecting the fraction of crab that did not grow and presumably didn't molt
      ####--triu(mT,k=1) is an upper triangular matrix with 0's along the diagonal reflecting the fractions of crab that molted and grew
      ###--crab that molted and grew also return to post-molt age 0
      ###--crab that did not molt and grow have post-molt age incremented by 1 (until max age tracked)
      np_ng_mpz =      diag(lst$mT)  * n_mpz;#--crab that did not grow       : need to increment post-molt age by 1
      np_gr_mpz = (1.0-diag(lst$mT)) * n_mpz;#--crab that will molt and grow : will need to decrement post-molt age to 0
      if (lst$maturity=="pre-molt size"){
        #--molt to maturity depends on premolt size
        ##--apply probability that molt is to maturity based on premolt size
        np_mat_gr_mpz =      lst$vPrMat  * np_gr_mpz; #--crab that mature         : need to update maturity state
        np_imm_gr_mpz = (1.0-lst$vPrMat) * np_gr_mpz; #--crab that stay immature  : no need to update maturity state
        ##--apply growth (only) transitions
        np_mat_gr_mpz = triu(lst$mT,k=1) %*% np_mat_gr_mpz;
        np_imm_gr_mpz = triu(lst$mT,k=1) %*% np_imm_gr_mpz;
      } else {
        #--molt to maturity depends on postmolt size
        ##--apply growth (only) transitions
        np_gr_mpz = triu(lst$mT,k=1) %*% np_gr_mpz;
        ##--apply probability that molt is to maturity based on postmolt size
        np_mat_gr_mpz =      lst$vPrMat  * np_gr_mpz; #--crab that mature         : need to update maturity state
        np_imm_gr_mpz = (1.0-lst$vPrMat) * np_gr_mpz; #--crab that stay immature  : no need to update maturity state
      }
      #--update post-molt ages
      np_ng_mpz     = lst$mPtoP1 %*% np_ng_mpz;     #--p incremented to p+1
      np_mat_gr_mpz = lst$mPto0  %*% np_mat_gr_mpz; #--p decremented to 0
      np_imm_gr_mpz = lst$mPto0  %*% np_imm_gr_mpz; #--p decremented to 0
      #--update maturity state
      np_mat_gr_mpz = lst$mImmToMat  %*% np_mat_gr_mpz;
      ##--re-combine vectors
      ret[[1]] = np_ng_mpz + np_imm_gr_mpz + np_mat_gr_mpz;
    } else {
      #--growth based on separate matrices for probabilities of molting and growth
      ##--apply probability of molting (doesn't update size or post-molt age)
      np_nmlt_mpz = (1.0-lst$vPrMolt) * np_mpz; #--crab that will not molt: need to increment postmolt age
      np_mlt_mpz  =      lst$vPrMolt  * np_mpz; #--crab that will molt    : need to decrement postmolt age to 0
      if (lst$maturity=="pre-molt size"){
        #--molt to maturity depends on premolt size
        ##--apply probability that molt is to maturity based on premolt size
        np_mat_mlt_mpz = lst$vPrMat       * np_mlt_mpz; #--need to update maturity state
        np_imm_mlt_mpz = (1.0-lst$vPrMat) * np_mlt_mpz; #--no need to update maturity state
        ##--apply growth transitions
        np_mat_mlt_mpz = lst$mPrGrMat %*% np_mat_mlt_mpz;
        np_imm_mlt_mpz = lst$mPrGrImm %*% np_imm_mlt_mpz;
        ##--update post-molt ages
        np_nmlt_mpz    = lst$mPtoP1 %*% np_nmlt_mpz;  # postmolt age p -> p+1
        np_mat_mlt_mpz = lst$mPto0 %*% np_mat_mlt_mpz;# postmolt age p -> 0
        np_imm_mlt_mpz = lst$mPto0 %*% np_imm_mlt_mpz;# postmolt age p -> 0
        ##--update maturity state
        np_mat_mlt_mpz = lst$mImmToMat %*% np_mat_mlt_mpz;# maturity state updated to mature
        ##--recombine
        ret[[1]] = np_nmlt_mpz + np_mat_mlt_mpz + np_imm_mlt_mpz;
      } else
      if (lst$maturity=="post-molt size") {
        #--molt to maturity depends on postmolt size
        ##--apply growth transitions
        np_mlt_mpz = lst$mPrGr %*% np_mlt_mpz;
        ##--apply probability that molt is to maturity based on postmolt size
        np_mat_mlt_mpz = lst$vPrMat       * np_mlt_mpz; #--need to update maturity state
        np_imm_mlt_mpz = (1.0-lst$vPrMat) * np_mlt_mpz; #--no need to update maturity state
        ##--update post-molt ages
        np_nmlt_mpz    = lst$mPtoP1 %*% np_nmlt_mpz;  # postmolt age p -> p+1
        np_mat_mlt_mpz = lst$mPto0 %*% np_mat_mlt_mpz;# postmolt age p -> 0
        np_imm_mlt_mpz = lst$mPto0 %*% np_imm_mlt_mpz;# postmolt age p -> 0
        ##--update maturity state
        np_mat_mlt_mpz = lst$mImmToMat %*% np_mat_mlt_mpz;# maturity state updated to mature
        ##--recombine
        ret[[1]] = np_nmlt_mpz + np_mat_mlt_mpz + np_imm_mlt_mpz;
      } else {
        stop("Error: `maturity` should be 'pre-molt size' or 'post-molt size' but was '",lst$maturity,"'.");
      }
    }
  }
  return(ret);
}

#' @title Create a matrix to advance maturity state to the next level
#' @description Function to create a matrix to advance maturity state from `m` to `m`+1 for crab.
#' @param dimMPZ - subset of a SparseDimsMap reflecting maturity states (m), post-molt age (p), and sizes (z) for a single sex
#' @return a permutation and aggregation matrix of 1's and 0's to update a vector with maturity from `m` to `m+1`
#' @details The resulting matrix can be used to update maturity state for crab that molted.
#'
#' vp = m %*% v;
#'
#' where `m` is the matrix, `v` is an abundance vector characterized by m,p,z, and `vp` is the resulting
#' abundance vector with the maturity state updated.
#'
#' @importFrom dplyr cross_join
#' @importFrom dplyr mutate
#' @export
#'
createMatrixNextMaturityState<-function(dimMPZ){
  nr = nrow(dimMPZ);
  levMs  = levels(dimMPZ$m);  #--levels of m
  unqMs  = unique(dimMPZ$m);  #--indices of m levels actually in dimMPZ
  nUMs = length(unqMs);
  nxtM = unqMs[c(2:nUMs,nUMs)];
  names(nxtM) = as.character(unqMs);
  dimMPZ_MPZ = dplyr::cross_join(dimMPZ |> dplyr::mutate(idx=dplyr::row_number()),
                                 dimMPZ |> dplyr::mutate(idx=dplyr::row_number())) |>
                 dplyr::mutate(val=as.numeric(((m.x==nxtM[as.character(m.y)])&(p.y==p.x)&(z.y==z.x))));
  mMPZ_MPZ = matrix(data=dimMPZ_MPZ$val,nrow=nr,ncol=nr,byrow=TRUE);
  rownames(mMPZ_MPZ) = paste0(dimMPZ$m," ",dimMPZ$p," ",dimMPZ$z);
  colnames(mMPZ_MPZ) = rownames(mMPZ_MPZ);
  return(mMPZ_MPZ);
}

#' Create a matrix to return post-molt category from `p` to 0 for crab that molted
#' @description Function to create a matrix to sreturn post-molt category from `p` to `0`for crab that molted.
#' @param dimMPZ - subset of a SparseDimsMap reflecting maturity states (m), post-molt age (p), and sizes (z) for a single sex
#' @return a permutation and aggregation matrix of 1's and 0's to update a vector with post-molt age states from `p` to `0`
#' @details The resulting matrix can be used to update post-molt age states for crab that molted.
#'
#' vp = m %*% v;
#'
#' where `m` is the matrix, `v` is an abundance vector characterized by m,p,z, and `vp` is the resulting
#' abundance vector with the post-molt states updated.
#'
#' @importFrom dplyr cross_join
#' @importFrom dplyr mutate
#' @export
#'
createMatrixPto0<-function(dimMPZ){
  nr = nrow(dimMPZ);
  ps = levels(dimMPZ$p);#--need to make sure this is min level w/in dimMPZ, not min of all levels (e.g., p=0,1 for imm; 0,1,2,3,4 for mat)
  dimMPZ_MPZ = dplyr::cross_join(dimMPZ |> dplyr::mutate(idx=dplyr::row_number()),
                                 dimMPZ |> dplyr::mutate(idx=dplyr::row_number())) |>
                 dplyr::mutate(val=as.numeric(((p.y==ps[1])&(m.y==m.x)&(z.y==z.x))));
  mMPZ_MPZ = matrix(data=dimMPZ_MPZ$val,nrow=nr,ncol=nr);
  rownames(mMPZ_MPZ) = paste0(dimMPZ$m," ",dimMPZ$p," ",dimMPZ$z);
  colnames(mMPZ_MPZ) = rownames(mMPZ_MPZ);
  return(mMPZ_MPZ);
}
#' @title Create a matrix to shift post-molt category from `p` to `p+1` (or max(`p`) as an accumulator state) for crab that didn't molt
#' @description Function to create a matrix to shift post-molt category from `p` to `p+1` (or max(`p`) as an accumulator state) for crab that didn't molt.
#' @param dimMPZ - subset of a SparseDimsMap reflecting maturity states (m), post-molt age (p), and sizes (z) for a single sex
#' @return a permutation and aggregation matrix of 1's and 0's to update a vector with post-molt age states from `p` to `p`+1 for crab that didn't molt.
#' @details The resulting matrix can be used to update post-molt age states for crab that did not undergo a molt.
#'
#' vp = m %*% v;
#'
#' where `m` is the matrix, `v` is an abundance vector characterized by m,p,z, and `vp` is the resulting
#' abundance vector with the post-molt states updated.
#'
#' @importFrom dplyr cross_join
#' @importFrom dplyr mutate
#' @export
#'
createMatrixPtoPplus1<-function(dimMPZ){
  nr     = nrow(dimMPZ);
  levPs  = levels(dimMPZ$p);  #--levels of p
  unqPs  = unique(dimMPZ$p);  #--indices of p levels actually in dimMPZ
  nUPs   = length(unqPs);
  nxtP   = unqPs[c(2:nUPs,nUPs)];
  names(nxtP) = as.character(unqPs);
  dimMPZ_MPZ  = dplyr::cross_join(dimMPZ |> dplyr::mutate(idx=dplyr::row_number()),
                                  dimMPZ |> dplyr::mutate(idx=dplyr::row_number())) |>
                dplyr::mutate(val=as.numeric(((p.x==nxtP[as.character(p.y)])&(m.y==m.x)&(z.y==z.x))));
  mMPZ_MPZ = matrix(data=dimMPZ_MPZ$val,nrow=nr,ncol=nr,byrow=TRUE);
  rownames(mMPZ_MPZ) = paste0(dimMPZ$m," ",dimMPZ$p," ",dimMPZ$z);
  colnames(mMPZ_MPZ) = rownames(mMPZ_MPZ);
  return(mMPZ_MPZ);
}
