#'
#' @title Create a template for the information necessary to specify the FE components of a model parameter
#' @description Function to create a template for the information necessary to specify the FE components of a model parameter.
#' @param formula - formula, possibly containing fixed effects
#' @param model_frame - model frame for parameter
#' @param contrasts - a list of contrasts (if not NULL)
#' @return a character vector
#'
#' @importFrom dplyr mutate row_number
#' @importFrom stringr fixed str_split_fixed
#' @importFrom tibble as_tibble
#'
#' @export
#'
makeFileTemplate_ParamsFE<-function(formula,
                                    model_frame,
                                    contrasts=NULL,
                                    drop.unused.levels=TRUE,
                                    debug=TRUE){
  lst = createParamsInfo_FEs(formula,
                             model_frame,
                             contrasts=contrasts,
                             drop.unused.levels=drop.unused.levels);
  str = "";
  if (is.null(lst)){
    str = c(str,"##--NO FIXED EFFECTS TERMS----");
    return(str);
  }
  for (trm in lst$feTerms){
    #--testing: trm = lst$reTerms[1];
    str = c(str,paste(trm,"   #--FE term"))
    strm = stringr::str_split_fixed(lst$re_pvs_lbls[[trm]],pattern=stringr::fixed("@"),n=3)
    colnames(strm) = c("term","param","group_level");
    dfr = tibble::as_tibble(strm) |> dplyr::select(-1) |>
            dplyr::mutate(value=0.0) |> dplyr::mutate(idx=dplyr::row_number(),.before=1);
    str = c(str,paste(nrow(dfr),"#--number of rows",collapse="\t\t"));
    str = c(str,paste(names(dfr),collapse="\t"));
    for (i in 1:nrow(dfr)) str = c(str,paste(dfr[i,],collapse="\t"));
  }
  if (debug) cat(str,sep="\n")
  return(str);
}
#'
#' @title Create a template for the information necessary to specify the RE components of a model parameter
#' @description Function to create a template for the information necessary to specify the RE components of a model parameter.
#' @param formula - formula, possibly containing random effects
#' @param model_frame - model frame for parameter
#' @param contrasts - a list of contrasts (if not NULL)
#' @return a character vector
#'
#' @importFrom dplyr mutate row_number
#' @importFrom stringr fixed str_split_fixed
#' @importFrom tibble as_tibble
#'
#' @export
#'
makeFileTemplate_ParamsRE<-function(formula,
                                    model_frame,
                                    contrasts=NULL,
                                    drop.unused.levels=TRUE,
                                    debug=TRUE){
  lst = createParamsInfo_REs(formula,
                             model_frame,
                             contrasts=contrasts,
                             drop.unused.levels=drop.unused.levels);
  str = "##--specifications for random effects-------------";
  if (is.null(lst)){
    str = c(str,"0 \t\t #--number of terms defining random effects\n##--NO RANDOM EFFECTS TERMS--------");
    return(str);
  }
  str = c(str,paste(length(lst$reTerms)," \t#--number of terms defining random effects"));
  ctr = 1;
  for (trm in lst$reTerms){
    #--testing: trm = lst$reTerms[1];
    #--RE parameter values for term
    str = c(str,paste("##--RE term",ctr,"----"))
    str = c(str,paste(trm,"\t#--RE term"))
    str = c(str,"##---parameter values----");
    strm = stringr::str_split_fixed(lst$re_pvs_lbls[[trm]],pattern=stringr::fixed("@"),n=3)
    colnames(strm) = c("term","param","group");
    dfr = tibble::as_tibble(strm) |> dplyr::select(-1)|>
            dplyr::mutate(idx=dplyr::row_number(),.before=1)|>
            dplyr::mutate(value=0.0);
    str = c(str,paste(0,"#--phase",collapse="\t\t"));
    str = c(str,paste(nrow(dfr),"#--number of rows",collapse="\t\t"));
    str = c(str,paste(names(dfr),collapse="\t"));
    for (i in 1:nrow(dfr)) str = c(str,paste(dfr[i,],collapse="\t"));

    #--covariance parameter values for term
    ##--testing: trm = lst$reTerms[1];
    str = c(str,"##---covariance parameter values----");
    str = c(str,paste(lst$covtype[[trm]],"\t\t#--covariance type"));
    ncps = attr(lst$lstQtp[[trm]],"num_covpars");
    str = c(str,paste(ncps,"\t\t#--number of covariance parameters"));
    str = c(str,paste(c("idx","mirror","value","phase"),collapse="\t"));
    for (i in 1:ncps){
      str = c(str,paste(i,0,0,0,collapse="\t\t"));
    }
  }
  if (debug) cat(str,sep="\n")
  return(str);
}
