#'
#' @title Get the 1-d index associated with a vector of dimension indices for a fully-crossed expansion
#' @description Function to get the 1-d index associated with a vector of dimension indices
#' for a fully-crossed expansion of the dimensions.
#' @param v - vector of dimension indices
#' @param n - vector of the number of corresponding dimension levels
#' @param rev - flag to get the index for the reverse order
#' @return an integer, the 1-d index into the fully-crossed dimensions array
#' @details If \code{rev} is false (the default), the index is calculated with the
#' first dimension cycling the fastest; otherwise, the index is calculated with the
#' last dimension cycling the fastest.
#'
#' @example inst/examples/example-getDenseDimsIndex.R
#'
#' @export
getDenseDimsIndex<-function(v,n,rev=FALSE){
    idx = 0; mi = 1;
    if (!rev){
        for (i in seq(1,length(v))){
            idx = idx + (v[i]-1)*mi;
            mi  = mi*n[i]
        }
    } else {
        for (i in rev(seq(1,length(v)))){
            idx = idx + (v[i]-1)*mi;
            mi  = mi*n[i]
        }
    }
    return(unname(idx)+1);
}

#'
#' @title Create a dataframe representation of a nested dimensions list by
#' recursively traversing it
#'
#' @description Function to create a dataframe representation of a
#' nested dimensions list by recursively traversing it.
#'
#' @param lst0 - dimensions list (or vector of levels for single dimension) to traverse
#' @param name - character string to associate with the terminal nesting level
#' @param level - counter for outer nesting level (leave at 0)
#' @param debug - flag (T/F) to print debugging info (default is FALSE)
#'
#' @return a tibble (see [tibble::tibble()]) with extra attribute "dmnms" copied
#' from \code{lst0} (or the value of \code{name} if \code{lst0} does not have a
#' "dmnms" attribute).
#'
#' @details Returns a tibble representation of the nested list by recursively
#' "unraveling" nested lists until the terminal nested level is reached.
#' The dimensions list should have a "dmnms" attribute that is a character vector
#' with the names of the associated dimensions (a vector of levels does not require
#' a dmnms attribute; the dimension name/dmnms attribute for the resulting tibble
#' will be the value of \code{name}).
#'
#' @example inst/examples/example-traverseDimsList.R
#'
#' @importFrom dplyr bind_rows select
#' @importFrom tibble tibble
#'
#' @export
#'
traverseDimsList<-function(lst0,name,level=0,debug=FALSE){
    col = paste0("v",level);
    if (debug) cat(paste0(col,":"),name,"\n");
    if (inherits(lst0,"list")) {
        dfr = NULL;
        for (name0 in names(lst0)){
            res   = traverseDimsList(lst0[[name0]],name0,level=level+1);
            level = res$level;#--reset level
            dfr1  = res$dfr;
            dfr1[[col]] = name;
            dfr   = dplyr::bind_rows(dfr,dfr1); #--copy
        }
    } else {
        if (debug) cat("\t",lst0,"\n");
        dfr = tibble::tibble(val=lst0,name=name);
        names(dfr)[2] = col;
        if (level>0) return(list(level=level-1,dfr=dfr));
        #--level==0 drops out to last bit
    }
    if (level==0) {
        dfr = dfr |> dplyr::select(rev(seq(1,ncol(dfr)-1,1)));
        if (is.null(attr(lst0,"dmnms"))){
            names(dfr) =  name;
        } else {
            names(dfr) = attr(lst0,"dmnms");
        }
        return(dfr);
    }
    return(list(level=level-1,dfr=dfr));
}

#'
#' @title Create a sparse index map (a tibble) for a set of (possibly nested) dimension levels
#'
#' @description Function to create a sparse 1-d index (a tibble) for a set of
#' (possibly nested) dimension levels. The index value can then be used to identify
#' the associated level for each dimension.
#'
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> name/value dimension pairs
#' @param debug - flag (T/F) to print debugging info
#'
#' @return a `DimsMap` object, which is a tibble (see [tibble::tibble()]) with
#' class "DimsMap" and attributes
#' \itemize{
#' \item{"dmnms" - vector of dimension names}
#' \item{"dmlvs" - list with non-nested dimension levels, by dimension name}
#' \item{"dmlns" - vector of dimension lengths, by dimension name}
#' }
#'
#' @details Each name in
#' \code{...} should be the name of the terminal dimension in the associated value (a
#' possibly-nested dimension list or a vector of levels). For each name, the
#' associated value is transformed to a tibble using [traverseDimsList()], then all the
#' tibbles are crossed to create a tibble with only the desired combinations of
#' actual dimension levels. The first column ("i") contains the 1-d index associated
#' with each combination of actual dimension levels.
#'
#' @example inst/examples/example-createSparseDimsMap.R
#'
#' @importFrom dplyr cross_join filter mutate row_number select
#' @importFrom rlang list2
#' @importFrom stringr str_sub str_subset
#'
#' @export
#'
createSparseDimsMap<-function(...,debug=FALSE){
    dfr = NULL;
    dots = rlang::list2(...);
    for (nm in names(dots)){
        #--testing: nm = names(dots)[1];
        dfrp = traverseDimsList(dots[[nm]],nm,debug=debug);
        if (is.null(dfr)) {
            dfr = dfrp;
            #--convert all dimensions to factors
            for (nm in names(dfr)) dfr[[nm]] = factor(dfr[[nm]]);
        } else {
            dfr = dfr |> dplyr::cross_join(dfrp,suffix=c("LFT","RGT"));
            nmsL = stringr::str_sub(stringr::str_subset(names(dfr),"LFT$"),end=-4);
            if (length(nmsL)>0){
                #--duplicate column names (ignoring trailing LFT/RGT) -> over-expanded rows
                #----keep only rows where ?LFT==?RGT
                for (nm in nmsL){
                  dfr = dfr |> dplyr::filter(.data[[paste0(nm,"LFT")]]==.data[[paste0(nm,"RGT")]]);
                }
                # #--remove "LFT" from column names, drop ?RGT columns
                nmsL = stringr::str_sub(stringr::str_subset(names(dfr),"LFT$"),end=-4);
                nmsA = stringr::str_subset(names(dfr),"RGT$",negate=TRUE);
                dfr  = dfr |> dplyr::select(dplyr::all_of(nmsA));
                nmsA = ifelse(nmsA %in% paste0(nmsL,"LFT"),
                              stringr::str_sub(nmsA,end=-4),nmsA);
                names(dfr) = nmsA;
            }
            #--convert all dimensions to factors
            for (nm in names(dfr)) dfr[[nm]] = factor(dfr[[nm]]);
        }#--!is.null(dfr);
    }#--nm
    dmsn = names(dfr);
    dmsl = list();
    dmsi = vector("numeric",length(dmsn));
    names(dmsi) = dmsn;
    for (n in dmsn){
        dmsl[[n]] = levels(dfr[[n]]);
        dmsi[n]   = length(dmsl[[n]]);
    }
    if (is.null(dfr)) return(NULL);
    attr(dfr,"dmtyp") <-"sparse"; #--dim type
    attr(dfr,"dmnms") <-dmsn;#--dim names
    attr(dfr,"dmlvs") <-dmsl;#--dim levels
    attr(dfr,"dmlns") <-dmsi;#--dim lengths
    dfrp = dfr |> createDimsFactors() |>
                 dplyr::arrange(dplyr::pick((tidyselect::everything())));
    dfr = dfrp |> dplyr::mutate(sparse_idx=dplyr::row_number(),.before=1);
    class(dfr) = c("DimsMap",class(dfr)[!("DimsMap" %in% class(dfr))]);#--add class attribute
    return(dfr);
}

#'
#' @title Create an index map using all (non-nested, fully-crossed) dimension levels
#'
#' @description Function to create an index map using all
#' (non-nested, fully-crossed) dimension levels
#'
#' @param map - a DimsMap object created using [createSparseDimsMap()]
#' @param debug - flag (T/F) to print debugging info
#'
#' @return a DimsMap object, which is atibble (see [tibble::tibble()]) with class "DimsMap"
#' and attributes
#' \itemize{
#' \item{"dmnms" - vector of dimension names}
#' \item{"dmlvs" - list with non-nested dimension levels, by dimension name}
#' \item{"dmlns" - vector of dimension lengths, by dimension name}
#' }
#' These attributes have the same values as the corresponding ones in \code{map}.
#'
#' @details The sparse dims map \code{map} is used to identify all levels of
#' each dimension included in the tibble. These are then used to create a
#' fully-crossed (dense) dims map (i.e., ignoring any nesting of dimensions).
#' The first column ("dense_idx") contains the 1-d index associated
#' with each combination of the full dimension levels.
#'
#' @example inst/examples/example-createDenseDimsMap.R
#'
#' @importFrom dplyr cross_join mutate row_number
#' @importFrom stringr str_sub str_subset
#'
#' @export
#'
createDenseDimsMap<-function(map,debug=FALSE){
    dmlvs = attr(map,"dmlvs");
    dfr = NULL;
    for (nm in names(dmlvs)) {
        #--nm = names(lst)[3];
        dfrp = traverseDimsList(dmlvs[[nm]],nm,debug=debug);
        if (is.null(dfr)){
            dfr = dfrp;
        } else {
            dfr  = dfr |> dplyr::cross_join(dfrp);
        }
    }
    attr(dfr,"dmtyp") <-"full";           #--dim type
    attr(dfr,"dmnms") <-attr(map,"dmnms");#--dim names
    attr(dfr,"dmlvs") <-dmlvs;            #--dim levels
    attr(dfr,"dmlns") <-attr(map,"dmlns");#--dim lengths
    dfr = dfr |> createDimsFactors() |>
                 dplyr::arrange(dplyr::pick((tidyselect::everything())));
    dfr = dfr |> dplyr::mutate(dense_idx=dplyr::row_number(),.before=1)
    class(dfr) = c("DimsMap",class(dfr)[!("DimsMap" %in% class(dfr))]);#--add class attribute
    return(dfr);
}

#'
#' @title Expand a dimensions map by additional dimensions
#' @description Function to expand a dimensions map by additional dimensions.
#' @param m1 - dimensions map to expand
#' @param m2 - dimensions map to expand by
#' @return a dimensions map with expanded dimensions (a tibble with class "DimsMap")
#' @details This function uses [dplyr::cross_join()] to cross join the two dimension maps to create
#' a map with expanded dimensions. The two maps should not have any dimensions in common.
#'
#' This is *not* a S3 generic method, so you must use `expand.DimsMap(m1,m2)`.
#'
#' If both maps are "dense", the result is a "dense" map.
#' If either is "sparse", the result is a "sparse" map.
#' The other map attributes are appropriately expanded.
#'
#' IMPORTANT NOTE: the (sparse or dense) index of `m1` is replicated `n` times, where `n` is
#' the number of rows in `m2`. It is the USER'S RESPONSIBILITY to re-number the index column as
#' appropriate for further use.
#'
#' @export
#'
expand.DimsMap<-function(m1,m2){
    m1_dmtyp<-attr(m1,"dmtyp"); #--dim type
    m1_dmnms<-attr(m1,"dmnms"); #--dim names
    m1_dmlvs<-attr(m1,"dmlvs"); #--dim levels
    m1_dmlns<-attr(m1,"dmlns"); #--dim lengths

    m2_dmtyp<-attr(m2,"dmtyp"); #--dim type
    m2_dmnms<-attr(m2,"dmnms"); #--dim names
    m2_dmlvs<-attr(m2,"dmlvs"); #--dim levels
    m2_dmlns<-attr(m2,"dmlns"); #--dim lengths

    dfr = dplyr::cross_join(m1,
                            m2 |> dplyr::select(all_of(m2_dmnms)));
    attr(dfr,"dmtyp") <-ifelse((m1_dmtyp=="dense")&&(m2_dmtyp=="dense"),"dense","sparse"); #--dim type
    attr(dfr,"dmnms") <-c(m1_dmnms,m2_dmnms);#--dim names
    attr(dfr,"dmlvs") <-c(m1_dmlvs,m2_dmlvs);#--dim levels
    attr(dfr,"dmlns") <-c(m1_dmlns,m2_dmlns);#--dim lengths
    class(dfr) = c("DimsMap",class(dfr)[!(class(dfr) %in% "DimsMap")]);#--add class attribute
    return(dfr);
}

#'
#' @title Test if an object is a DimsMap
#'
#' @description Function to test if an object is a DimsMap.
#'
#' @param x - object to test
#'
#' @return logical
#'
#' @details None.
#'
#' @export
#'
is.DimsMap<-function(x){
  return(inherits(x,"DimsMap"))
}

#'
#' @title Test if an object is a sparse DimsMap
#'
#' @description Function to test if an object is a sparse DimsMap.
#'
#' @param x - object to test
#'
#' @return logical
#'
#' @details None.
#'
#' @export
#'
is.SparseDimsMap<-function(x){
  res = FALSE;
  if (inherits(x,"DimsMap")){
    if (attr(dfr,"dmtyp")=="sparse") res = TRUE;
  }
  return(res);
}

#'
#' @title Test if an object is a dense DimsMap
#'
#' @description Function to test if an object is a dense DimsMap.
#'
#' @param x - object to test
#'
#' @return logical
#'
#' @details None.
#'
#' @export
#'
is.DenseDimsMap<-function(x){
  res = FALSE;
  if (inherits(x,"DimsMap")){
    if (attr(dfr,"dmtyp")=="dense") res = TRUE;
  }
  return(res);
}

#'
#' @title Create maps between sparse and dense dimension indices
#'
#' @description Function to create maps between sparse and dense dimension indices.
#'
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> name/value dimension pairs
#' @param debug - flag (T/F) to print debugging info
#'
#' @return a list with elements
#' \itemize{
#' \item{dfrS2D - tibble with mapping from sparse to dense (full) index values}
#' \item{dfrD2S - tibble with mapping from dense (full) to sparse index values}
#' }
#'
#' @details Each name in
#' \code{...} should be the name of the terminal dimension in the associated value (a
#' possibly-nested dimension list or a vector of levels), as in [createSparseDimsMap()].
#' In \code{dfrS2D}, the row index corresponds to the 1's-based sparse index value.
#' In \code{dfrF2S}, the row index corresponds to the 1's-based full index value.
#'
#' @example inst/examples/example-createDimsMaps.R
#'
#' @importFrom dplyr inner_join full_join mutate select
#'
#' @export
#'
createDimsMaps<-function(...,debug=FALSE){
    dfrSprs = createSparseDimsMap(...,debug=debug);
    dfrDens = createDenseDimsMap(dfrSprs,debug=debug);
    dmnms   = attr(dfrSprs,"dmnms");
    dfrS2D  = dfrSprs |> dplyr::inner_join(dfrDens,by=dmnms) |>
                dplyr::select(-1) |>
                dplyr::select(dense_idx,!dplyr::last_col());# |>
                #createDimsFactors();
    dfrD2S  = dfrDens |> dplyr::full_join(dfrSprs,by=dmnms) |>
                dplyr::mutate(sparse_idx=ifelse(is.na(sparse_idx),-1,sparse_idx)) |>
                dplyr::select(-1) |>
                dplyr::select(sparse_idx,!dplyr::last_col());# |>
                #createDimsFactors();
    return(list(dfrS2D=dfrS2D,dfrD2S=dfrD2S));
}

#'
#' @title Create the intrinsic model dimensions from the user model dimensions
#'
#' @description Function to create the intrinsic model dimensions from the user model dimension.
#'
#' @param udfr - dataframe representing user dimensions map (sparse or full)
#'
#' @return a dataframe representing the intrinsic model dimensions.
#'
#' @details The user dimensions defined by `udfr` are expanded to the full set of
#' intrinsic model dimensions ("y","s","r","x","m","a","p","z"), with default values
#' ("all") added for intrinsic dimensions undefined in the user dimensions.
#'    * `y` - year index
#'    * `s` - season index
#'    * `r` - region index
#'    * `x` - sex index
#'    * `m` - maturity state index
#'    * `a` - age index
#'    * `p` - post-molt age (or shell condition) index
#'    * `z` - size index
#' @example inst/examples/example-createIntrinsicDimensions.R
#'
#' @importFrom dplyr mutate select
#'
#' @export
#'
createIntrinsicDims<-function(udfr){
  #--intrinsic dimensions info
  idmnms = c("y","s","r","x","m","a","p","z");
  idefs  = c(y="all",
             s="all",
             r="all",
             x="all",
             m="all",
             a="all",
             p="all",
             z="all");

  #--user dimensions info
  utype  = attr(udfr,"dmtyp");  #--user dimensions type (sparse or dense)
  #----determine missing intrinsic dimensions
  udmnms = attr(udfr,"dmnms"); #--vector of user dim names
  mdmnms = idmnms[!(idmnms %in% udmnms)];#--missing intrinsic dimensions
  #----other dimensions info
  idmlvs = attr(udfr,"dmlvs"); #--list of user dim levels
  idmlns = attr(udfr,"dmlns"); #--vector of user dim lengths
  #--add columns, info for missing dimensions
  idfr = udfr;
  if (length(mdmnms)>0){
    for (dmnm in mdmnms){
      idfr = idfr |> dplyr::mutate("{dmnm}":=factor(idefs[dmnm]));
      idmlvs[dmnm] = idefs[dmnm];
      idmlns[dmnm] = 1;
    }
  }
  #--rearrange into canonical format
  idfr = idfr |> dplyr::select(tidyselect::all_of(idmnms));
  idmlvs = idmlvs[idmnms];
  idmlns = idmlns[idmnms];
  attr(idfr,"dmtyp") = paste("intrinsic-",utype);#--dimensions type
  attr(idfr,"dmnms") = idmnms; #--vector of dimension names
  attr(idfr,"dmlvs") = idmlvs; #--list of dimension levels
  attr(idfr,"dmlns") = idmlns; #--vector of dimension lengths
  return(idfr);
}
