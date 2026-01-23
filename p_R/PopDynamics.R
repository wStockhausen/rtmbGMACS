##--NOTE: these are functions from rsimTCSAM and are included simply as references
#'
#' @title Step the population numbers one season
#' @description Function to step the population numbers one season
#'
#' @param y - year index
#' @param s - season index
#'
#' @details
#' The population is propagated forward one season (i.e.) time step. Potential processes involved
#' in the process are:
#'  - molting/change in shell condition
#'  - maturation
#'  - growth
#'  - natural mortality
#'  - fishing mortality
#'  -
#'
#' @export
#'
stepPopNs<-function(y,s,N_z){
  #--calculate natural, fishing, and total mortality
  Z_z = M_z*dt;
  for (f in 1:nFlts){
    Fc_fz[f,] = Sel[f,]*F[f];                     #--fishery capture rate in f
    Fr_fz[f,] = R_fz[f,] * Fc_fz[f,];             #--retained catch mortality rate in f
    Fd_fz[f,] = (1-R_fz[f,]) * dmr[f] * Fc_fz[f,];#--discard catch mortality rate in f
    Ft_fz[f,] = Fr_fz[f,]+Fd_fz[f,];              #--discard catch mortality rate in f
    Ft_z = Ft_z + Fm_fz[f,];#--total fishing mortality
  }
  Z_z = M_z + Ft_z;

  #--calculate survival probabilities
  S_z = exp(-Z_z*dt);
  if (mort_first[y,i]){
    #--mortality first, then growth
    Np_z = N_z * S_z;    #--survivors: matrix multiplication not necessary
    D_z  = N_z * (1-S_z);#--total dead
    C_z  = 0;
    for (f in 1:nFlts){
      Ct_fz[f,] = (Fc_fz[f,]/Z_z) * D_z;#--total number captured in fishery f
      Rm_fz[f,] = (Fr_fz[f,]/Z_z) * D_z;#--retained catch mortality in fishery f
      Dm_fz[f,] = (Fd_fz[f,]/Z_z) * D_z;#--discard catch mortality in fishery f
      Cm_fz[f,] = (Ft_fz[f,]/Z_z) * D_z;#--total catch mortality in fishery f
      Cm_z = Cm_z + Cm_fz[f,];#--total catch mortality across all fisheries
    }
    #--determine catch mortality in fisheries
    if (doGrowth[y,s]){
      #--apply stage changes due to molting/shell condition/maturity/growth
      Np_z  = Np_z %*% StgTrns_zz;#--note matrix multiplication
    }
  } else {
    #--growth first, then mortality
    if (doGrowth[y,s]){
      #--apply stage changes due to molting/shell condition/maturity/growth
      Np_z = N_z %*% StgTrns_zz;#--note matrix multiplication
    }
    Np_z = Np_z * S_z;    #--survivors: matrix multiplication not necessary
    D_z  = Np_z * (1-S_z);#--total dead
    C_z  = 0;
    for (f in 1:nFlts){
      Ct_fz[f,] = (Fc_fz[f,]/Z_z) * D_z;#--total number captured in fishery f
      Rm_fz[f,] = (Fr_fz[f,]/Z_z) * D_z;#--retained catch mortality in fishery f
      Dm_fz[f,] = (Fd_fz[f,]/Z_z) * D_z;#--discard catch mortality in fishery f
      Cm_fz[f,] = (Ft_fz[f,]/Z_z) * D_z;#--total catch mortality in fishery f
      Cm_z = Cm_z + Cm_fz[f,];#--total catch mortality across all fisheries
    }
  }
  if (doRec[y,s]){
    Np_z = Np_z + Recr_z;
  }
  return(list(N_z=Np_z,Ct_fz=Ct_fz,Rm_fz=Rm_fz,Dm_fz=Dm_fz,Cm_fz=Cm_fz,Cm_z=Cm_z));
}

runAllYears<-function(Ni_z){
  N_z = calcInitN();
  N_ysz = array(0,dim=c(nY,nS,nZ),dimnames=c("y","s","z"));
  N_ysz[1,1,] = N_z;
  for (y in stYr:endYr){
      N_ysz[1,,] = 1; #--TODO: FIX ME
  }

}

runOneYear<-function(y,N_z){
  N_sz[1,]    = 0;#--need to dimension these and set to 0
  Ct_sfz[1,,] = 0;
  Rm_sfz[1,,] = 0;
  Dm_sfz[1,,] = 0;
  Cm_sfz[1,,] = 0;
  Cm_sz[1,]   = 0;
  N_sz[1,]    = N_z;#--need to dimension these
  for (s in 1:nS){
    lst = stepPopNs(y,s,N_z);
    N_sz[s+1,]    = lst$N_z;
    Ct_sfz[s+1,,] = lst$Ct_fz;
    Rm_sfz[s+1,,] = lst$Rm_fz;
    Dm_sfz[s+1,,] = lst$Dm_fz;
    Cm_sfz[s+1,,] = lst$Cm_fz;
    Cm_sz[s+1,]   = lst$Cm_z;
  }
  return(list(N_sz=N_sz,Ct_sfz=Ct_sfz,
              Rm_sfz=Rm_sfz,Dm_sfz=Dm_sfz,
              Cm_sfz=Cm_sfz,Cm_sz=Cm_sz));
}

