
mapParameterIndicesToFunctions<-function(lstCTL,lstModMtx){
  dfr1 = lstCTL$dfrFcnsExp |> dplyr::select(!par_keys) |>
           dplyr::inner_join(lstCTL$dfrMapParsToFcns,by="fcn_idx",relationship="many-to-many");

  lstDFRs = list();
  for (par_nm_ in names(lstModMtx)){
    #--testing: par_nm_ = names(lstModMtx)[1];
    lstMM = lstModMtx[[par_nm_]];
    par_key_ = lstMM$par_key;
    par_id_  = lstMM$par_id;
    lstIDXs = list();
    for (idxType_ in c("FEs","ECs","REs")){
      #--testing: idxType_ = "FEs";
      if (!is.null(lstMM[[idxType_]]$dfrMF))
        lstIDXs[[idxType_]] = lstMM[[idxType_]]$dfrMF |>
                                dplyr::mutate(par_nm=par_nm_,
                                              par_id=par_id_,.before=1) |>
                                dplyr::mutate(idxType=idxType_);
    }
    dfrIDXs = dplyr::bind_rows(lstIDXs);
    lstDFRs[[par_nm_]] = dfr1 |> dplyr::filter(par_key==par_key_) |>
                           dplyr::inner_join(dfrIDXs);
    rm(dfrIDXs,lstIDXs);
  }
  dfrDFRs = dplyr::bind_rows(lstDFRs) |>
              tidyr::pivot_wider(names_from="idxType",values_from="row_idx");
  if (!("FEs" %in% names(dfrDFRs))) dfrDFRs$FEs = NA_integer_;
  if (!("ECs" %in% names(dfrDFRs))) dfrDFRs$ECs = NA_integer_;
  if (!("REs" %in% names(dfrDFRs))) dfrDFRs$REs = NA_integer_;

  #--
  udfr = dfrDFRs |> dplyr::select(!c(par_key,par_id,par_nm,FEs,ECs,REs)) |> dplyr::distinct();
  jnby = names(udfr);
  lst = list();
  for (r in 1:nrow(udfr)){
    #--testing: r = 1;
    rw = udfr[r,];
    dfrp = rw |> dplyr::inner_join(dfrDFRs,by=jnby) |>
            dplyr::select(c(par_id,par_key,par_nm,FEs,ECs,REs));
    lstp = list();
    for (rp in 1:nrow(dfrp)){
      #--testing: rp = 1;
      rwp = dfrp[rp,];
      lstp[[rwp$par_id]] = as.list(rwp);
    }
    rw = rw |> dplyr::mutate(par_idxs=list(lstp));
    lst[[r]] = rw;
  }
  dfrMap = dplyr::bind_rows(lst);
  return(dfrMap);
}
