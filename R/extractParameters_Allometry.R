#--extract alllometry parameters
#'
#' @title Extract allometry parameters from parameter info list
#' @description Function to extract allometry parameters from parameter info list.
#' @param lst - parameter info list from [readParamInfo_Allometry()]
#' @param dims - dimensions list
#' @return a list
#' @details
#' @export
#'
extractParameters_Allometry<-function(lst,
                                      dims){
  if (FALSE){
    #TODO: delete this section--for development only when not running function
    dims = inputs$dims;
    lst = params_info$allometry;
  }
  lstAlls = NULL;
  if (!is.null(dims)) lstAlls = alls_GetLevels(dims);

  #--create `map` list for allometry
  map = list();

  #--expand allometry parameter information
  if (lst$option=="function"){

  } else if (lst$option=="data"){
    #--inputs are fixed values, so phase = -1
    tform = stringr::str_starts(lst$transform,"tf_") |>
              ifelse(lst$transform,paste0("tf_",lst$transform));
    eval(parse(text=paste0("tf<-",tform)));#--evaluate transformation
    if (is.character(lst$dfr$IV[1])){
      #--need to evaluate `value`s to get (fixed) parameter values
      dfrp = lst$dfr |> dplyr::select(!c(z,IV));
      dfrv = lst$dfr |> dplyr::select(IV);
      lst1 = list();
      for (r in 1:nrow(lst$dfr)){
        lst1[[r]] =  dfrp[r,] |>
                   dplyr::cross_join(
                     tibble::as_tibble(eval(parse(text=dfrv[r,]))) |> tidyr::pivot_longer(tidyselect::everything(),names_to="z")
                   )
      }
      dfr = dplyr::bind_rows(lst1);
      rm(dfrp,dfrv,lst1);
    } else {
      dfr = lst$dfr;
    }
    #--transform values, add bounds, phase < 0, priors info, parameter index
    dfr = dfr |>
            dplyr::mutate(IV=tf(IV),
                          LB=-Inf,
                          UB= Inf,
                          phz=-1,
                          PriorType="none",
                          Pr1=0,
                          Pr2=999) |>
            dplyr::mutate(pidx=row_number(),.before=1); #--parameter index
    #--extract parameter values
    pWatZ = dfr$IV;#--weight-at-size in kg
    map[["pWatZ"]] = factor(NA+pWatZ);

      |>
              expandDataframe(lstAlls=lstAlls);
    }
  } else {
    stop("option not recognized.")
  }
}
