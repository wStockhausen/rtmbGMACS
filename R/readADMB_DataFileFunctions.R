#--functions to read input data file sections
#'
#' @title Read a slightly modified ADMB GMACS input data file
#' @description Function to read a slightly modified ADMB GMACS input data file.
#' @param fn - path to input data file
#' @return list
#' @details The input ADMB GMACS file should be modified in the following manner:
#' \itemize{
#' \item{The initial section (prior to the Catch Data section) should be removed or commented out.}
#' \item{The comment character "#" should be removed from the header line for input Catch Data, Relative Abundance Indices, Size Comps, Growth Data, and Environmental Indices dataframes}
#' }
#'
#' The output list has an attribute "type" = "input data list" with elements
#' \itemize{
#'   \item{dfrCDF - list of catch data dataframes}
#'   \item{dfrID  - list of index data dataframes}
#'   \item{dfrZCD - list of size comps dataframes}
#'   \item{dfrGrD - list of growth (molt increment) data dataframes}
#'   \item{dfrMOD - list of maturity ogive data dataframes}
#'   \item{dfrTD - list of tagging data dataframes}
#'   \item{dfrED - list of environmental data dataframes}
#' }
#' @export
#'
readADMB_DataFile<-function(fn){
  if (file.exists(fn)){
    lns = readLines(fn);

    dfrCD = readCatchData(lns);

    dfrID = readIndexData(lns);

    dfrZCD = readSizeComps(lns);

    dfrGrD = readGrowthData(lns);#--TODO: should this return a list of dataframes?

    #--maturity ogives
    #----TODO: add ability to read maturity ogives
    warning("Reading maturity ogive data has not yet been implemented!");
    dfrMOD = NULL;#--placeholder

    #--tagging data
    #----TODO: add ability to read tagging data
    dfrTD = NULL;#--placeholder
    warning("Reading tagging data has not yet been implemented!");

    dfrED = readEnvironmentalData(lns);
  } else {
    stop(paste0("File '",fn,"' does not exist!"));
  }
  lst = list(dfrCD=dfrCD,dfrID=dfrID,dfrZCD=dfrZCD,
             dfrGrD=dfrGrD,
             dfrMOD=dfrMOD,dfrTD=dfrTD,
             dfrED=dfrED);
  attr(lst,"type") = "input data list";
  return(lst);
}

#'
#' @title Read a catch data frame in GMACS input format 1
#' @importFrom dplyr mutate
#'
readCatchDataFrameFormat1<-function(lnsp){
  nlnsp = length(lnsp);
  ilnp = 1;
  ilnp = findNextLine(lnsp,ilnp+1);
  units_type = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  catch_type = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  fleet = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  x = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  m = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  s = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  nr = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  dfr = parseTextToDataframe(lnsp,ilnp,nr,col_names=TRUE) |>
          dplyr::mutate(f=fleet,x=x,m=m,s=s,
                        catch_type=tolower(catch_type),
                        units_type=tolower(units_type));
  return(list(dfr=dfr,iln=ilnp+nr-1));
}

#'
#' @title Read the Catch Data section of an input data file
#' @return a tibble (see [tibble::tibble]) with attribute "type" = "catch data dataframe"
#' @importFrom dplyr bind_rows
readCatchData<-function(lns){
  nlns = length(lns);
  iln = findKeyword(lns,"Catch Data");
  if (is.null(iln)) return(NULL);
  #--read input format
  iln = findNextLine(lns,iln+1);
  input_format = parseVal(lns[iln]);
  iln = findNextLine(lns,iln+1);
  nCDFs = parseVal(lns[iln]);
  if (nCDFs==0) return(NULL);
  lst = list();
  iln = iln+1;
  if (input_format==1){
    for (iCDF in 1:nCDFs) {
      res = readCatchDataFrameFormat1(lns[iln:nlns]);
      lst[[iCDF]] = res$dfr;
      iln = res$iln+1+iln;
    }
  } else {
    #TODO: read other formats?
  }
  dfr = dplyr::bind_rows(lst);
  attr(dfr,"type") = "catch data dataframe";
  return(dfr);
}

#'
#' @title Read an index data frame in GMACS input format 1
#' @importFrom dplyr mutate
#'
readIndexDataFrameFormat1<-function(lnsp){
  nlnsp = length(lnsp);
  ilnp = 1;
  ilnp = findNextLine(lnsp,ilnp+1);
  units_type = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  catch_type = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  fleet = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  x = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  m = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  s = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  nr = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  dfr = parseTextToDataframe(lnsp,ilnp,nr,col_names=TRUE) |>
          dplyr::mutate(f=fleet,x=x,m=m,s=s,
                        catch_type=tolower(catch_type),
                        units_type=tolower(units_type));
  return(list(dfr=dfr,iln=ilnp+nr-1));
}

#'
#' @title Read the Index Data section of an input data file
#' @return a tibble (see [tibble::tibble]) with attribute "type" = "relative abundance data dataframe"
#' @importFrom dplyr bind_rows
#'
readIndexData<-function(lns){
  nlns = length(lns);
  iln = findKeyword(lns,"Index Data");
  if (is.null(iln)) return(NULL);
  #--read input format
  iln = findNextLine(lns,iln+1);
  input_format = parseVal(lns[iln]);
  iln = findNextLine(lns,iln+1);
  nRADs = parseVal(lns[iln]);
  if (nRADs==0) return(NULL);
  lst = list();
  iln = iln+1;
  if (input_format==1){
    for (iRAD in 1:nRADs) {
      res = readIndexDataFrameFormat1(lns[iln:nlns]);
      lst[[iRAD]] = res$dfr;
      iln = res$iln+1+iln;
    }
  } else {
    #TODO: read other formats?
  }
  dfr = dplyr::bind_rows(lst);
  attr(dfr,"type") = "relative abundance dataframes";
  return(dfr);
}

#'
#' @title Read a size comps frame in GMACS input format 1
#' @importFrom dplyr mutate
#'
readSizeCompsDataFrameFormat1<-function(lnsp){
  nlnsp = length(lnsp);
  ilnp = 1;
  ilnp = findNextLine(lnsp,ilnp+1);
  catch_type = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  fleet = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  x = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  m = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  s = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  nr = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  nc = parseVal(lnsp[ilnp]);
  ilnp = findNextLine(lnsp,ilnp+1);
  dfr = parseTextToDataframe(lnsp,ilnp,nr,col_names=TRUE) |>
          dplyr::mutate(f=fleet,x=x,m=m,s=s,
                        catch_type=tolower(catch_type));
  return(list(dfr=dfr,iln=ilnp+nr-1));
}

#'
#' @title Read the Size Comps section of an input data file
#' @return a tibble (see [tibble::tibble]) with attribute "type" = "size comps dataframe"
#' @importFrom dplyr bind_rows
#'
readSizeComps<-function(lns){
  nlns = length(lns);
  iln = findKeyword(lns,"Size Comps");
  if (is.null(iln)) return(NULL);
  #--read input format
  iln = findNextLine(lns,iln+1);
  input_format = parseVal(lns[iln]);
  iln = findNextLine(lns,iln+1);
  nZCs = parseVal(lns[iln]);
  if (nZCs==0) return(NULL);
  lst = list();
  iln = iln+1;
  if (input_format==1){
    for (iZC in 1:nZCs) {
      res = readSizeCompsDataFrameFormat1(lns[iln:nlns]);
      lst[[iZC]] = res$dfr;
      iln = res$iln+1+iln;
    }
  } else {
    #--TODO: read other formats?
  }
  dfr = dplyr::bind_rows(lst);
  attr(dfr,"type") = "size comps dataframe";
  return(dfr);
}

#'
#' @title Read the Growth Data section of an input data file
#' @return a tibble (see [tibble::tibble]) with attribute "type" = "growth data dataframe"
#' @import dplyr
#'
readGrowthData<-function(lns){
  nlns = length(lns);
  iln = findKeyword(lns,"Growth Data");
  if (is.null(iln)) return(NULL);
  #--read input format
  iln = findNextLine(lns,iln+1);
  input_format = parseVal(lns[iln]);
  iln = findNextLine(lns,iln+1);
  nObs = parseVal(lns[iln]);
  if (nObs==0) return(NULL);
  iln = findNextLine(lns,iln+1);
  lst = list();
  xs = c("male","female");
  if (input_format==1){
  dfr = parseTextToDataframe(lns,iln,nObs,col_names=TRUE) |>
          dplyr::rowwise() |>
          dplyr::mutate(x=xs[sex]) |>
          dplyr::ungroup() |>
          dplyr::select(x,premolt_size,molt_inc,cv=CV);
  } else {
    #TODO: read other formats?
  }
  attr(dfr,"type") = "growth data dataframe";
  return(dfr);
}

#'
#' @title Read the Environmental Data section of an input data file
#' @return a tibble (see [tibble::tibble]) with attribute "type" = "environmental data dataframes"
#' @import dplyr
#'
readEnvironmentalData<-function(lns){
  nlns = length(lns);
  iln = findKeyword(lns,"Environmental Data");
  if (is.null(iln)) return(NULL);
  #--read input format
  iln = findNextLine(lns,iln+1);
  input_format = parseVal(lns[iln]);
  iln = findNextLine(lns,iln+1);
  nEIs = parseVal(lns[iln]);
  if (nEIs==0) return(NULL);
  iln = findNextLine(lns,iln+1);
  if (input_format==1){
    dfr1 = parseTextToDataframe(lns,iln,nEIs);
    iln = iln+nEIs;
    iln = findNextLine(lns,iln+1);
    #--determine number of rows in environmental indices dataframe
    nr = 0; for (rw in 1:nrow(dfr1)) nr = nr + (dfr1$`stop-year`-dfr1$`start-year`+1);
    dfr2 = parseTextToDataframe(lns,iln,nr);
    #--convert merged dataframe into list of individual dataframes, by env. index
    #--NOTE: possibly a dplyr function to do this in one step
    strEIs = dfr2 |> dplyr::distinct(`index/name`) |> dplyr::pull(`index/name`);
    lst = list();
    for (strEI in strEIs){
      lst[[strEI]] = dfr2 |> dplyr::filter(`index/name`==strEI);
    }
  } else {
    #TODO: read other formats?
  }
  dfr = dplyr::bind_rows(lst);
  attr(dfr,"type") = "environmental data dataframe";
  return(dfr);
}


