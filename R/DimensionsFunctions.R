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
#' @example inst/examples/example-getFullDimsIndex.R
#'
#' @export
getFullDimsIndex<-function(v,n,rev=FALSE){
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
#'
#' @return a tibble (see [tibble::tibble()]) with attributes
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
createSparseDimsMap<-function(...){
    dfr = NULL;
    dots = rlang::list2(...);
    for (nm in names(dots)){
        #--testing: nm = names(dots)[1];
        dfrp = traverseDimsList(dots[[nm]],nm,debug=TRUE);
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
    attr(dfr,"dmnms") <-dmsn;#--dim names
    attr(dfr,"dmlvs") <-dmsl;#--dim levels
    attr(dfr,"dmlns") <-dmsi;#--dim lengths
    dfrp = dfr |> createDimsFactors() |>
                 dplyr::arrange(dplyr::pick((tidyselect::everything())));
    dfr = dfrp |> dplyr::mutate(sparse_idx=dplyr::row_number(),.before=1);
    return(dfr);
}

#'
#' @title Create an index map using all (non-nested, fully-crossed) dimension levels
#'
#' @description Function to create an index map using all
#' (non-nested, fully-crossed) dimension levels
#'
#' @param map - tibble created using [createSparseDimsMap()]
#'
#' @return a tibble (see [tibble::tibble()]) with attributes
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
#' The first column ("j") contains the 1-d index associated
#' with each combination of the full dimension levels.
#'
#' @example inst/examples/example-createFullIndexMap.R
#'
#' @importFrom dplyr cross_join mutate row_number
#' @importFrom stringr str_sub str_subset
#'
#' @export
#'
createFullDimsMap<-function(map){
    dmlvs = attr(map,"dmlvs");
    dfr = NULL;
    for (nm in names(dmlvs)) {
        #--nm = names(lst)[3];
        dfrp = traverseDimsList(dmlvs[[nm]],nm);
        if (is.null(dfr)){
            dfr = dfrp;
        } else {
            dfr  = dfr |> dplyr::cross_join(dfrp);
        }
    }
    attr(dfr,"dmnms") <-attr(map,"dmnms");#--dim names
    attr(dfr,"dmlvs") <-dmlvs;            #--dim levels
    attr(dfr,"dmlns") <-attr(map,"dmlns");#--dim lengths
    dfr = dfr |> createDimsFactors() |>
                 dplyr::arrange(dplyr::pick((tidyselect::everything())));
    dfr = dfr |> dplyr::mutate(dense_idx=dplyr::row_number(),.before=1)
    return(dfr);
}

#'
#' @title Create maps between sparse and dense dimension indices
#'
#' @description Function to create maps between sparse and dense dimension indices.
#'
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> name/value dimension pairs
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
createDimsMaps<-function(...){
    dfrSprs = createSparseDimsMap(...);
    dfrDens = createFullDimsMap(dfrSprs);
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
