#--extract survey catchability functions parameters
#'
#' @title Extract survey catchability functions parameters from parameter info list
#' @description Function to extract survey catchability functions parameters from parameter info list.
#' @param lst - parameter info list from [readParamInfo_SurveysCatchability()]
#' @param dims - `dims` list from `setupModelDims`
#' @return a list (see details)
#' @details The format of the returned list depends on the `option` specified in `lst`.
#' If `lst$option` is "pre-specified" (i.e., the input parameter values are fixed values of weight-at-size
#' specified at least for every population category), the output list has elements
#' \itemize{
#' \item{option - format option}
#' \item{params - vector of weight-at-size values}
#' \item{map - vector for the RTMB `map` list indicating the parameters are fixed; i.e. not to be estimated}
#' \item{dfrIdx2Pars - dataframe with columns "pidx", dimension names (y,s,r,x,m,p,z), and parameter specs IV,LB,UB,phz,PriorType,Pr1,Pr2
#' indicating the un-expanded mapping from `pidx` to the parameter input values ("IV"). `pidx` is the index into the parameter vector.
#' The dimension columns are incidental and reflect the specification from the input text file.}
#' \item{dfrDims2Pars - dataframe with columns "pidx" and the dimension names (y,s,r,x,m,p,z).
#' The dimension levels are fully expanded such that each row corresponds to a single "multidimensional" level.
#' The value of `pidx` for each row indicates the index of the associated weight-at-size value in the parameter vector.}
#' }
#'
#' If `lst$option` is "function", the output list has elements
#' \itemize{
#' }
#'
#' @import dplyr
#'
#' @export
#'
extractParamInfo_SurveysCatchability<-function(lst,
                                               dims=NULL,
                                               verbose=TRUE){
  if (verbose) message("starting extractParameters_SurveysCatchability.")
  if (FALSE){
    #--NOTE: run this section if you are just stepping through the code for development purposes
    ##--assumes you've created `dims` via `dims = setupModelDims();`
    verbose = TRUE;
    lst = res;#--assumed output from `res  = readParamInfo_SurveysCatchability(conn,verbose=FALSE);`
  }

  #--expand Surveys parameter information----
  if (tolower(lst$option)=="function"){
    ##--option == "function"----
    ##--inputs are functions and parameters definitions
    out = extractParamInfoFunctionType1(lst,dims$dmsYSC,"Surveys_Catchability",
                                        xtra_cols=c("flt","sel_idx","avl_idx"),
                                        verbose=verbose);
    out$flts = lst$flts;
  } else if (tolower(lst$option)=="pre-specified"){
    ##--option == "pre-specified"----
    ###--inputs are pre-specified (fixed values)
    tform = stringr::str_starts(lst$transform,"tf_") |>
              ifelse(lst$transform,paste0("tf_",lst$transform));
    eval(parse(text=paste0("tf<-",tform)));#--evaluate transformation
    if (is.character(lst$dfr$value[1])){
      ###--values are character strings----
      ###--need to evaluate and transform `value`s to get (fixed) parameter values
      if (verbose){
        message("in extractParamInfo_Surveys: lst$dfr:")
        print(lst$dfr);
      }
      dfrp = lst$dfr |> dplyr::select(!c(z,value));
      dfrv = lst$dfr |> dplyr::select(value);
      lst1 = list();
      for (r in 1:nrow(lst$dfr)){
        lst1[[r]] =  dfrp[r,] |>
                   dplyr::cross_join(
                     tibble::as_tibble(eval(parse(text=dfrv[r,]))) |>
                       tidyr::pivot_longer(tidyselect::everything(),names_to="z",values_to="value")
                   )
      }
      dfr = dplyr::bind_rows(lst1);
      rm(dfrp,dfrv,lst1);
    } else {
      ####--values are numeric----
      #####--just copy lst4dfr
      dfr = lst$dfr;
    }
    ###--transform values, add bounds, phase < 0, priors info, parameter index----
    dfr = dfr |>
            dplyr::mutate(value=tf(value),
                          LB=-Inf,
                          UB= Inf,
                          phz=-1,
                          PriorType="none",
                          Pr1=0,
                          Pr2=999) |>
            dplyr::mutate(pidx=dplyr::row_number(),.before=1) |> #--parameter index
            dplyr::rename(IV=value);
    ###--add missing dimensions as "all"s
    dmnms=attr(dims$dmsYSC,"dmnms");
    for (dmnm in dmnms){
      if (!(dmnm %in% names(dfr))) dfr[[dmnm]] = "all";
    }
    dfr = dfr |> dplyr::select(pidx,flt,y,s,r,x,m,p,z,IV,LB,UB,phz,PriorType,Pr1,Pr2) |>
            dplyr::mutate(sel_idx=0,    #--no need for a sel function
                          avl_idx=0,    #--no need for a avl function
                          .after=flt);
    ###--extract parameter values----
    if (verbose) {
      message("in extractParameters_SurveysCatchability: resolved 'pre-specified' option values");
      print(dfr);
    }
    pSrvQ_FPs = dfr$IV;#--survey catchability value
    map = list(pSrvQ_FPs=factor(NA+pSrvQ_FPs));#--NAs indicate fixed values
    if (verbose) message("in extractParameters_Surveys: expanding dataframe.")
    ###--create list of all dimension levels in `dims$dmsYSC` to convert "all"'s to pop levels----
    ####--listAlls is a list with all individual dimension levels, by individual dimension y, s, r, x, m, a, p, z
    lstAlls = NULL;
    if (!is.null(dims)) lstAlls = alls_GetLevels(dims$dmsYSC,verbose=verbose);
    dfrp = dfr |>
              dplyr::select(!c(IV,LB,UB,phz,PriorType,Pr1,Pr2)) |>
              expandDataframe(lstAlls=lstAlls,verbose=verbose);#--expanded for "alls"
    out = list(option="pre-specified",
               flts=lst$flts,
               params=pSrvQ_FPs,
               map=map,
               dfrIdx2Pars=dfr,
               dfrDims2Pars=dfrp);
  } else {
    ##--option not recognized----
    stop("option not recognized.")
  }
  return(out);
}

if (FALSE){
  #--test "pre-specified" option with "vertical" format-----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Transforms.R"))
  source(file.path(dirPrj,"R/readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R/readParamInfo_SurveysCatchability.R"))
  source(file.path(dirPrj,"R/extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R/extractParamInfo_SurveysCatchability.R"))
  source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))
  dims = setupModelDims();
  conn = file.path(dirPrj,"testing/testSurveysCatchability/inputSpecs_SurveysCatchability.pre-specified-vertical.txt");
  res1v  = readParamInfo_SurveysCatchability(conn,verbose=FALSE);
  res2v = extractParamInfo_SurveysCatchability(res1v,dims,verbose=FALSE);
  View(res2v$dfrDims2Pars);
  View(res2v$dfrIdx2Pars);
}

if (FALSE){
  #--test "pre-specified" option with "horizontal" format-----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Transforms.R"))
  source(file.path(dirPrj,"R/readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R/readParamInfo_SurveysCatchability.R"))
  source(file.path(dirPrj,"R/extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R/extractParamInfo_SurveysCatchability.R"))
  source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))
  dims = setupModelDims();
  conn = file.path(dirPrj,"testing/testSurveysCatchability/inputSpecs_SurveysCatchability.pre-specified-horizontal.txt");
  res1h  = readParamInfo_SurveysCatchability(conn,verbose=FALSE);
  res2h = extractParamInfo_SurveysCatchability(res1h,dims,verbose=FALSE);
  View(res2h$dfrDims2Pars);
  View(res2h$dfrIdx2Pars);
}

if (FALSE){
  #--test "function" option----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Transforms.R"))
  source(file.path(dirPrj,"R/readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R/readParamInfo_SurveysCatchability.R"))
  source(file.path(dirPrj,"R/extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R/extractParamInfo_SurveysCatchability.R"))
  source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))
  dims = setupModelDims();
  conn = file.path(dirPrj,"testing/testSurveysCatchability/inputSpecs_SurveysCatchability.function.txt");
  res1f  = readParamInfo_SurveysCatchability(conn,verbose=FALSE);
  res2f = extractParamInfo_SurveysCatchability(res1f,dims,verbose=TRUE);
  View(res2f$Fcns);
  View(res2f$dfrUniqCmbs);
  View(res2f$dfrUHCs);
}
