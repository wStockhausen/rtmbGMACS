##--DEFUNCT!!
# mapFunctionsToData_Allometry<-function(dfrZW1,dfrMapPIsToFcns){
#   ##--map functions/parameters to data----
#   dfrZWp = dfrZW1  |>
#              dplyr::left_join(dfrFMapPIsToFcns) |>
#              dplyr::arrange(obs_id);
#   lstZWp = list();
#   for (nm in names(lstModMtx)){
#     #--testing: nm = names(lstModMtx)[1];
#     lstModMtxPar = lstModMtx[[nm]];
#     if (lstModMtxPar$FEs$npars>0) {
#       dfrZWppFEs = dfrZWp |>
#                      dplyr::mutate(par_nm=nm,
#                                    fcn_idx=lstModMtxPar$fcn_idx,fcn_id=lstModMtxPar$fcn_id,fcn_nm=lstModMtxPar$fcn_nm,
#                                    par_idx=lstModMtxPar$par_idx,par_id=lstModMtxPar$par_id,type="FEs") |>
#                      dplyr::inner_join(lstModMtxPar$FEs$dfrMF);
#     } else {dfrZWppFEs = NULL;}
#     if (lstModMtxPar$ECs$npars>0) {
#       dfrZWppECs = dfrZWp |>
#                      dplyr::mutate(par_nm=nm,
#                                    fcn_idx=lstModMtxPar$fcn_idx,fcn_id=lstModMtxPar$fcn_id,fcn_nm=lstModMtxPar$fcn_nm,
#                                    par_idx=lstModMtxPar$par_idx,par_id=lstModMtxPar$par_id,type="ECs") |>
#                      dplyr::inner_join(lstModMtxPar$ECs$dfrMF);
#     } else {dfrZWppECs = NULL;}
#     if (lstModMtxPar$REs$npars>0) {
#       dfrZWppREs = dfrZWp |>
#                      dplyr::mutate(par_nm=nm,
#                                    fcn_idx=lstModMtxPar$fcn_idx,fcn_id=lstModMtxPar$fcn_id,fcn_nm=lstModMtxPar$fcn_nm,
#                                    par_idx=lstModMtxPar$par_idx,par_id=lstModMtxPar$par_id,type="REs") |>
#                      dplyr::inner_join(lstModMtxPar$REs$dfrMF);
#     } else {dfrZWppREs = NULL;}
#     lstZWp[[nm]] = dplyr::bind_rows(dfrZWppFEs,dfrZWppECs,dfrZWppREs);# |> tidyr::pivot_wider(names_from="type",values_from="row_idx");
#     rm(lstModMtxPar,dfrZWppFEs,dfrZWppECs,dfrZWppREs);
#   }
#   dfrZWpp = dplyr::bind_rows(lstZWp) |>
#               tidyr::pivot_wider(names_from="type",values_from="row_idx") |>
#               dplyr::arrange(obs_id);
#   if (!("FEs" %in% names(dfrZWpp))) dfrZWpp$FEs = NA;
#   if (!("ECs" %in% names(dfrZWpp))) dfrZWpp$ECs = NA;
#   if (!("REs" %in% names(dfrZWpp))) dfrZWpp$REs = NA;
#   #dfrZWpp  = dfrZWpp |> dplyr::rowwise() |> dplyr::mutate(row_idxs = list(c(FEs,ECs,REs)));#--NOTE: might want to drop row_idxs?
#   udfrZWpp = dfrZWpp |> dplyr::distinct(fcn_idx,par_idx,FEs,ECs,REs);#--unique function/model parameter combinations
#   rm(nm);
#
#   dfrZWppp = dfrZWpp |>
#                dplyr::rowwise() |>
#                dplyr::mutate(par_info=list(list(par_id=par_id,par_idx=par_idx,FEs=FEs,ECs=ECs,REs=REs))) |>
#                dplyr::ungroup() |>
#                dplyr::select(!c(par_nm,par_idx,FEs,ECs,REs)) |>
#                tidyr::pivot_wider(names_from="par_id",values_from=par_info) |>
#                dplyr::mutate(pars=list(list()),
#                              prd=NA_real_,
#                              nll=NA_real_,
#                              z=as.numeric(z),      #--convert character to numeric
#                              obs=as.numeric(obs)); #--convert character to numeric
#  return(dfrZWppp);
# }
#
