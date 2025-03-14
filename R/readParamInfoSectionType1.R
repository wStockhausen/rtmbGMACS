#--read model parameter info section
#' @title Read a model parameter info section
#' @description Function to read a model parameter info section
#' @param conn - connection (or filename) to read from
#' @param verbose - flag to print extra info
#' @return nested list with elements (see **Details**)
#' @details If `option` = "function", the returned list has elements
#' \itemize{
#'   \item{option}
#'   \item{Fcns - list with elements}
#'   \itemize{
#'     \item{n - number of functions defined}
#'     \item{dfr - dataframe with function specs}
#'     \item{reflvls - list with number of reference levels definied (`n`) and a dataframe with the levels (`dfr`)}
#'   }
#'   \item{MPs - list with main parameters info; elements are}
#'   \itemize{
#'     \item{n - number of main parameters defined}
#'     \item{dfr - dataframe with main parameters specs}
#'   }
#'   \item{OPs - list with offset parameters info; elements are}
#'   \itemize{
#'     \item{n - number of offset parameters defined}
#'     \item{dfr - dataframe with offset parameters specs}
#'   }
#'   \item{DPs - list with devs parameters info; elements are}
#'   \itemize{
#'     \item{n - number of dev parameters defined}
#'     \item{dfr - dataframe with dev parameter specs}
#'     \item{reflvls - list with number of reference levels defined (`n`) and a dataframe with the levels (`dfr`)}
#'   }
#'   \item{REs - list with random effects parameters info; elements are}
#'   \itemize{
#'     \item{n - number of random effects parameters defined}
#'     \item{dfr - dataframe with random effects parameter specs}
#'     \item{reflvls - list with number of reference levels defined (`n`) and a dataframe with the levels (`dfr`)}
#'   }
#'   \item{ECs - list with environmental covariates info; elements are}
#'   \itemize{
#'     \item{n - number of environmental covariates defined}
#'     \item{dfr - dataframe with environmental covariates specs}
#'     \item{reflvls - list with number of reference levels defined (`n`) and a dataframe with the levels (`dfr`)}
#'   }
#'   \item{FPs - list with functional priors info; elements are}
#'   \itemize{
#'     \item{n - number of functional priors defined}
#'     \item{dfr - dataframe with functional priors specs}
#'   }
#' }
#'
#' If `option` = "pre-specified", a list with the following elements is returned:
#' \itemize{
#'   \item{option}
#'   \item{transform - name of function to transform input values to kg}
#'   \item{dfr - dataframe with value specs }
#' }
#'
#' @examples
#' # example code reading "vertical" pre-specified format
#' if (FALSE){
#'   conn="inputSpecs_Allometry.data_vertical.txt";
#'   lstResults = readParamInfo_Allometry(conn,verbose=TRUE);
#' }
#'
#' # example code reading "horizontal" pre-specified format
#' if (FALSE){
#'   conn="inputSpecs_Allometry.data_horizontal.txt";
#'   lstResults = readParamInfo_Allometry(conn,verbose=TRUE);
#' }
#'
#' # example code reading "function" format
#' if (FALSE){
#'   conn="inputSpecs_Allometry.function.txt";
#'   lstResults = readParamInfo_Allometry(conn,verbose=TRUE);
#' }
#'
#' @import dplyr
#' @import purrr
#' @import readr
#' @import stringr
#' @md
#' @export
#'
readParamInfoSectionType1<-function(lns,verbose=FALSE){
  iln = 1;
  #--parse input format option----
  if (verbose) cat("Starting readParamInfoSectionType1\n");
  lst     = extractTextSection(lns,1,1);
  option  = (stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--input format option
  #--parse text----
  if (tolower(option)=="function"){
    ##--parse "function" input option info----
    ###--parse function definitions----
    if (verbose) cat("reading function info\n");
    lst     = extractTextSection(lns,1,lst$end+1);
    nFcns   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of functions
    lst     = extractTextSection(lns,nFcns+1,lst$end+1);
    dfrFcns = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ function defs
    if (verbose) print(dfrFcns);
    lst     = extractTextSection(lns,1,lst$end+1);
    nLvls   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of reference levels
    if (nLvls>0){
      if (verbose) cat("reading function reference levels info\n");
      lst    = extractTextSection(lns,nLvls+1,lst$end+1);
      dfrRefLvls = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ function reference levels
      if (verbose) print(dfrRefLvls);
    } else {dfrRefLvls=NULL;}

    ##--parse parameter information (as necessary)----
    ###--Main parameters----
    if (verbose) cat("reading main parameters info\n");
    lst     = extractTextSection(lns,1,lst$end+1);
    nMPs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of main parameters
    if (nMPs>0){
      lst     = extractTextSection(lns,nMPs+1,lst$end+1);
      dfrMPs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ main parameter defs
      if (verbose) print(dfrMPs);
    } else {dfrMPs=NULL;}

    ###--parse offset parameters----
    if (verbose) cat("reading offset parameters info\n");
    lst     = extractTextSection(lns,1,lst$end+1);
    nOPs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of offset parameters
    if (nOPs>0){
      lst     = extractTextSection(lns,nOPs+1,lst$end+1);
      dfrOPs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ offset parameter defs
      if (verbose) print(dfrOPs);
    } else {dfrOPs=NULL;}

    ###--parse "devs" parameters----
    if (verbose) cat("reading devs parameters info\n");
    lst     = extractTextSection(lns,1,lst$end+1);
    nDPs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of devs parameters
    if (nDPs>0){
      lst     = extractTextSection(lns,nDPs+1,lst$end+1);
      dfrDPs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ devs parameter defs
      if (verbose) print(dfrDPs);
      ###--reference levels----
      lst     = extractTextSection(lns,1,lst$end+1);
      nLvls_DPs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of devs reference level defs
      if (nLvls_DPs>0){
        if (verbose) cat("reading devs parameter referece levels info\n");
        lst     = extractTextSection(lns,nLvls_DPs+1,lst$end+1);
        dfrRefLvls_DPs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ devs reference level def
        if (verbose) print(dfrRefLvls_DPs);
      } else {dfrRefLvls_DPs=NULL;}
    } else {dfrDPs=NULL; nLvls_DPs=0; dfrRefLvls_DPs=NULL;}

    ###--parse "RE" parameter specifications----
    if (verbose) cat("reading RE parameters info\n");
    lst     = extractTextSection(lns,1,lst$end+1);
    nREs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of REs parameter specifications
    if (nREs>0){
      lst     = extractTextSection(lns,nREs+1,lst$end+1);
      dfrREs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ RE parameter specifications
      if (verbose) print(dfrREs);
    } else {dfrREs=NULL;}

    ###--parse parameter-related environmental covariates----
    if (verbose) cat("reading parameter-related environmental covariates info\n");
    lst     = extractTextSection(lns,1,lst$end+1);
    nPECs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of environmental covariates
    if (verbose) cat("reading",nPECs,"parameter-related environmental covariate definitions.")
    if (nPECs>0){
      lst     = extractTextSection(lns,nPECs+1,lst$end+1);
      dfrPECs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ environmental covariate defs
      if (verbose) print(dfrPECs);
      ####--reference levels----
      lst     = extractTextSection(lns,1,lst$end+1);
      nLvls_PECs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of environmental covariate reference level defs
      if (verbose) cat("will read",nLvls_PECs,"parameter-related environmental covariate reference levels.")
      if (nLvls_PECs>0){
        if (verbose) cat("reading parameter-related environmental covariate reference levels info\n");
        lst     = extractTextSection(lns,nLvls_PECs+1,lst$end+1);
        dfrRefLvls_PECs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ environmental covariate reference level def
        if (verbose) print(dfrRefLvls_PECs);
      } else {dfrRefLvls_PECs=NULL;}
    } else {dfrPECs=NULL; nLvls_PECs=0; dfrRefLvls_PECs=NULL;}

    ###--parse function-related environmental covariates----
    if (verbose) cat("reading function-related environmental covariates info\n");
    lst     = extractTextSection(lns,1,lst$end+1);
    nFECs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of environmental covariates
    if (verbose) cat("reading",nFECs,"function-related environmental covariate definitions.")
    if (nFECs>0){
      lst     = extractTextSection(lns,nFECs+1,lst$end+1);
      dfrFECs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ environmental covariate defs
      if (verbose) print(dfrFECs);
      ####--reference levels----
      lst     = extractTextSection(lns,1,lst$end+1);
      nLvls_FECs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of environmental covariate reference level defs
      if (verbose) cat("will read",nLvls_FECs,"function-related environmental covariate reference levels.")
      if (nLvls_FECs>0){
        if (verbose) cat("reading function-related environmental covariate reference levels info\n");
        lst     = extractTextSection(lns,nLvls_FECs+1,lst$end+1);
        dfrRefLvls_FECs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ environmental covariate reference level def
        if (verbose) print(dfrRefLvls_FECs);
      } else {dfrRefLvls_FECs=NULL;}
    } else {dfrFECs=NULL; nLvls_FECs=0; dfrRefLvls_FECs=NULL;}

    ###--parse functional priors----
    if (verbose) cat("reading functional priors info\n");
    lst     = extractTextSection(lns,1,lst$end+1);
    nFPs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of functional priors
    if (verbose) cat("will read",nFPs,"functional priors definitions.")
    if (nFPs>0){
      if (verbose) cat("Reading functional priors info with",nFPs,"rows.");
      dfrFPs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ functional priors defs
      if (verbose) print(dfrFPs);
    } else {dfrFPs=NULL;}

    ###--create output list----
    if (verbose) cat("creating output list\n");
    out = list(option=option,
               Fcns=list(n=nFcns,dfr=dfrFcns,
                         reflvls=list(n=nLvls,dfr=dfrRefLvls)),        #--functions
               MPs=list(n=nMPs,dfr=dfrMPs),                            #--main parameters
               OPs=list(n=nOPs,dfr=dfrOPs),                            #--offset parameters
               DPs=list(n=nDPs,dfr=dfrDPs,
                        reflvls=list(n=nLvls_DPs,dfr=dfrRefLvls_DPs)), #--devs parameter vectors
               REs=list(n=nREs,dfr=dfrREs),                            #--RE parameter vectors
               PECs=list(n=nPECs,dfr=dfrPECs,
                        reflvls=list(n=nLvls_PECs,dfr=dfrRefLvls_PECs)), #--parameter-level env. covars
               FECs=list(n=nFECs,dfr=dfrFECs,
                        reflvls=list(n=nLvls_FECs,dfr=dfrRefLvls_FECs)), #--function-level env. covars
               FPs=list(n=nFPs,dfr=dfrFPs)                               #--functional priors
               );
  } else if (tolower(option)=="pre-specified") {
    ##--parse "pre-specified" input option info----
    if (verbose) cat("reading 'pre-specified' option info\n");
    lst  = extractTextSection(lns,1,lst$end+1);
    fmtHV = (stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--horizontal or vertical format?
    lst  = extractTextSection(lns,1,lst$end+1);
    nRws = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of rows to read
    lst  = extractTextSection(lns,1,lst$end+1);
    tfrm = (stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--transform function
    lst  = extractTextSection(lns,nRws+1,lst$end+1);
    dfr  = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ values
    ##--if tolower(fmtHV)=="vertical", no need to do anything(?)
    if (tolower(fmtHV)=="horizontal"){
      ###--convert horizontal format to vertical----
      nms = names(dfr);#--column names in dfr
      #--find index of name that starts with "values"
      idx = which(stringr::str_starts(nms,"values"));
      #--extract dimension identifier between the parentheses (e.g., "z" in "values(z):27")
      dim = extractTextBetweenParens(nms[idx]);
      #--extract the dimension value after ":" (e.g., "27" in "values(z):27")
      dmval = extractTextAfterString(nms[idx],":");
      #--extract all the dimension values from the column names
      ncols = length(nms);
      dmvals = c(dmval,nms[(idx+1):ncols]);
      #--convert to vertical format
      tblDM = tibble::as_tibble_col(dmvals,column_name=dim);
      tblLst = list();
      for (rw in 1:nRws){
        vals = as.vector(t(dfr[rw,idx:ncols]));
        tblLst[[rw]]  = dfr[rw,1:(idx-1)] |>
                          dplyr::cross_join(tblDM |> tibble::add_column(tibble::as_tibble_col(vals,column_name="value")));
      }
      dfr = dplyr::bind_rows(tblLst);
    }#--"horizontal" format
    ###--create output list----
    out = list(option=option,
               transform=tfrm,
               dfr=dfr);
  } else {
    stop(paste0("Option '",option,"' not recognized."));
  } #--option
  if (verbose) cat("Finished readParamInfoSectionType1\n");
  return(out);
}

