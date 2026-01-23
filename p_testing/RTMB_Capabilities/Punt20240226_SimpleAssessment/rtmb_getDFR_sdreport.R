getDFR_sdreport <- function(obj,report=TRUE,model=""){
  sdr    = sdreport(obj);
  #--determine max number of dimensions
  nDims = 0;
  if (report) {
    for (nm in names(sdr$env$ADreportDims)){
      nD = length(sdr$env$ADreportDims[[nm]]);
      if (nD>nDims) nDims = nD;
    }
  } else {
    for (nm in names(sdr$env$parameters)){
      dms = dim(sdr$env$parameters[[nm]]);
      nD  = ifelse(is.null(dms),1,length(dms));
      if (nD>nDims) nDims = nD;
    }
  }
  lstSDR = list();
  types = c("Est","Std")
  for (type in types){
    #--type=types[1];
    lst = as.list(sdr,type,report=report);
    for (nm in names(lst)){
      #--nm = names(lst)[1];
      nmp = paste0(nm,"_",type);
      var = lst[[nm]];
      if (is.null(dim(var))) dim(var)<-length(var);#--assign a dimension
      dms = dim(var);
      nDv = length(dms);
      dfrIdxs<-function(var){
        idxs = paste0("i",1:nDv);
        dfr = tidyr::expand_grid("{idxs[1]}":=(1:dms[1]));
        if (nDv>1){
          for (d in 2:length(dms)){
            vn = as.character(idxs[d]);
            dfr = dfr |> dplyr::cross_join(tidyr::expand_grid("{vn}":=(1:dms[d])));
          }
        }
        return(dfr);
      }
      lstSDR[[nmp]] = dplyr::bind_cols(
                        tibble::tibble(name=nm,nDim=nDv,type=type,value=as.vector(var)),
                        dfrIdxs(var)
                      );
    }
  }
  dfrSDR = dplyr::bind_rows(lstSDR) |> dplyr::mutate(model=model) |> 
             tidyr::pivot_wider(id_cols=c("model","name","nDim",paste0("i",1:nDims)),
                                names_from="type",values_from="value");
  return(dfrSDR);
}
