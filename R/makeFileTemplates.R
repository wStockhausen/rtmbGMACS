#'
#' @title Create a template for the information necessary to specify all components of a model parameter
#' @description Function to create a template for the information necessary to specify all components of a model parameter.
#' @param parameter - parameter name
#' @param x - formula or list from [createParamInfo()]
#' @param model_frame - model frame for parameter (required if `x` is a formula)
#' @param contrasts - a list of contrasts (if not NULL)
#' @param drop.unused.levels - flag to drop unused levels
#' @return a character vector
#'
#' @importFrom rlang is_formula
#'
#' @export
#'
makeFileTemplate_Param<-function(parameter,
                                 x,
                                 model_frame=NULL,
                                 contrasts=NULL,
                                 drop.unused.levels=TRUE,
                                 offset=NULL,
                                 debug=TRUE){
  if (rlang::is_formula(x)){
    lst = createParamInfo(parameter,
                          x,
                          model_frame,
                          contrasts=contrasts,
                          drop.unused.levels=drop.unused.levels,
                          offset=offset,
                          sparse=FALSE);
  } else {lst=x;}

  str = paste0("#-----Parameter values specification for ",parameter,deparse(x$formula));
  str = c(str,paste(parameter,"\t\t#--parameter name"));
  #--fixed effects
  txt = makeFileTemplate_ParamFEs(x$piFEs,debug=FALSE);
  str = c(str,txt);
  #--random effects
  txt = makeFileTemplate_ParamREs(x$piREs,debug=FALSE);
  str = c(str,txt);
  #--smooths
  txt = makeFileTemplate_ParamSMs(x$piSMs);  #--TODO: finish this function
  str = c(str,txt);
  if (debug) cat(str,sep="\n");
  return(str)
}

#'
#' @title Create a template for the information necessary to specify the FE components of a model parameter
#' @description Function to create a template for the information necessary to specify the FE components of a model parameter.
#' @param x - formula, possibly containing fixed effects, or list from [createParamInfo_FeEs()]
#' @param model_frame - model frame for parameter (required if `x` is a formula)
#' @param contrasts - a list of contrasts (if not NULL)
#' @param drop.unused.levels - flag to drop unused levels
#' @return a character vector
#'
#' @importFrom dplyr mutate row_number
#' @importFrom rlang is_formula
#' @importFrom tibble tibble
#'
#' @export
#'
makeFileTemplate_ParamFEs<-function(x,
                                    model_frame=NULL,
                                    contrasts=NULL,
                                    drop.unused.levels=TRUE,
                                    debug=TRUE){
  if (rlang::is_formula(x)){
    lst = createParamInfo_FEs(x,
                               model_frame,
                               contrasts=contrasts,
                               drop.unused.levels=drop.unused.levels,
                               sparse=FALSE);
  } else {lst=x;}

  str = "##--specifications for fixed effects-------------";
  if (is.null(lst)){
    str = c(str,"##--NO FIXED EFFECTS TERMS----");
    return(str);
  }
  trm = deparse(lst$fe_form)
  str = c(str,paste(trm,"   #--FE terms"))
    dfr = tibble::tibble(param=lst$fe_pvs_lbls) |>
            dplyr::mutate(value=0.0,
                          phase=0,
                          lower=-Inf,
                          upper=Inf,
                          prior="none",
                          p1=NA,
                          p2=NA) |>
            dplyr::mutate(idx=dplyr::row_number(),
                          mirror=0,
                          .before=1);
    str = c(str,paste(nrow(dfr),"#--number of rows",collapse="\t\t"));
    str = c(str,paste(names(dfr),collapse="\t"));
    for (i in 1:nrow(dfr)) str = c(str,paste(dfr[i,],collapse="\t"));

  if (debug) cat(str,sep="\n")
  return(str);
}

#'
#' @title Create a template for the information necessary to specify the RE values for a model parameter
#' @description Function to create a template for the information necessary to specify the RE values for a model parameter.
#' @param x - formula, possibly containing fixed effects, or list from [createParamInfo_REs()]
#' @param model_frame - model frame for parameter
#' @param contrasts - a list of contrasts (if not NULL)
#' @return a character vector
#'
#' @importFrom dplyr mutate row_number
#' @importFrom rlang is_formula
#' @importFrom stringr fixed str_split_fixed
#' @importFrom tibble as_tibble
#'
#' @export
#'
makeFileTemplate_ParamREs<-function(x,
                                       model_frame,
                                       contrasts=NULL,
                                       drop.unused.levels=TRUE,
                                       debug=TRUE){
  if (rlang::is_formula(x)){
    lst = createParamInfo_REs(x,
                              model_frame,
                              contrasts=contrasts,
                              drop.unused.levels=drop.unused.levels);
  } else {lst=x;}

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
    strm = stringr::str_split_fixed(lst$re_pvs_lbls[[trm]],pattern=stringr::fixed("}{"),n=3)
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
    str = c(str,paste("idx","mirror","value","phase","lower","upper","prior","p1","p2",sep="\t"));
    for (i in 1:ncps){
      str = c(str,paste(i,0,0,0,-Inf,Inf,"none",NA,NA,sep="\t"));
    }
    ctr = ctr+1;
  }#--trm loop
  if (debug) cat(str,sep="\n")
  return(str);
}

#'
#' @title Create a template for the information necessary to specify the smooths for a model parameter TODO: finish this function!
#' @description Function to create a template for the information necessary to specify the smooths for a model parameter.
#' @param x - formula, possibly containing fixed effects, or list from [createParamInfo_SMs()]
#' @param model_frame - model frame for parameter
#' @param contrasts - a list of contrasts (if not NULL)
#' @return a character vector
#'
#' @importFrom dplyr mutate row_number
#' @importFrom rlang is_formula
#' @importFrom stringr fixed str_split_fixed
#' @importFrom tibble as_tibble
#'
#' @export
#'
makeFileTemplate_ParamSMs<-function(x,
                                    model_frame,
                                    contrasts=NULL,
                                    drop.unused.levels=TRUE,
                                    debug=TRUE){
  if (rlang::is_formula(x)){
    lst = createParamInfo_SMs(x,
                              model_frame,
                              contrasts=contrasts,
                              drop.unused.levels=drop.unused.levels);
  } else {lst=x;}

  str = "##--specifications for smooths-------------";
  if (is.null(lst)){
    str = c(str,"0 \t\t #--number of terms defining smooths\n##--NO SMOOTH TERMS--------");
    return(str);
  }
}#--makeFileTemplate_ParamSMs


