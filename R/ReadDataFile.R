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
#' @export
#'
readDataFile<-function(fn){
  if (file.exists(fn)){
    lns = readLines(fn);
    lstCDFs = readCatchDataframes(lns);
    lstRADs = readRelativeAbundanceData(lns);
    lstZCs = readSizeComps(lns);
    dfrGrw  = readGrowthData(lns);
    #--maturity ogives
    lstEIs = readEnvironmentalData(lns);
  } else {
    stop(paste0("File '",fn,"' does not exist!"));
  }
  return(list(lstCDFs=lstCDFs,lstRADs=lstRADs,dfrGrw=dfrGrw,lstEIs=lstEIs));
}

findKeyword<-function(lns,kw){
  nlns = length(lns);
  iln = 1;
  while((iln<nlns)&&
        stringr::str_starts(toupper(stringr::str_trim(lns[iln])),toupper(kw),negate=TRUE)){
    iln=iln+1;
  }
  if (iln==nlns) return(NULL);#--didn't find keyword in lns
  return(iln);
}

parseVal<-function(str,type=NULL){
  #--split by whitespace and drop anything after comment character
  strp1 = scan(text=str,what=character(),comment.char="#",quiet=TRUE,strip.white=TRUE);
  #--parse remaining string
  val = readr::parse_guess(strp1);
  return(val);
}

#'
#' @title Find next uncommented line
#' @description Function to find next uncommented line in a vector of character strings.
#' @param lns - character vector to search
#' @param iln - index at which to start search
#' @return index into `lns` at which the next uncommented line occurs
#' @export
#'
findNextLine<-function(lns,iln){
  nlns = length(lns);
  while((iln<nlns)&&
        (stringr::str_length(stringr::str_trim(lns[iln]))==0||
        stringr::str_starts(toupper(stringr::str_trim(lns[iln])),"#"))){
    iln=iln+1;
  }
  if (iln==nlns) return(NULL);#--didn't find an uncommented line
  return(iln);
}

#'
#' @title Read a catch data frame in GMACS input format 1
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
  dfr = readDataframe(lnsp,ilnp,nr,col_names=TRUE) |>
          dplyr::mutate(f=fleet,x=x,m=m,s=s,
                        catch_type=tolower(catch_type),
                        units_type=tolower(units_type));
  return(list(dfr=dfr,iln=ilnp+nr-1));
}

#'
#' @title Read Catch Data section of input data file
readCatchDataframes<-function(lns){
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

  }
  return(lst);
}

#'
#' @title Read a catch data frame in GMACS input format 1
#'
readRelativeAbundanceDataFrameFormat1<-function(lnsp){
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
  dfr = readDataframe(lnsp,ilnp,nr,col_names=TRUE) |>
          dplyr::mutate(f=fleet,x=x,m=m,s=s,
                        catch_type=tolower(catch_type),
                        units_type=tolower(units_type));
  return(list(dfr=dfr,iln=ilnp+nr-1));
}

#'
#' @title Read Relative Abundance Data section of input data file
#'
readRelativeAbundanceData<-function(lns){
  nlns = length(lns);
  iln = findKeyword(lns,"Relative Abundance Data");
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
      res = readRelativeAbundanceDataFrameFormat1(lns[iln:nlns]);
      lst[[iRAD]] = res$dfr;
      iln = res$iln+1+iln;
    }
  } else {

  }
  return(lst);
}

#'
#' @title Read a catch data frame in GMACS input format 1
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
  dfr = readDataframe(lnsp,ilnp,nr,col_names=TRUE) |>
          dplyr::mutate(f=fleet,x=x,m=m,s=s,
                        catch_type=tolower(catch_type));
  return(list(dfr=dfr,iln=ilnp+nr-1));
}

#'
#' @title Read Size Comps section of input data file
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

  }
  return(lst);
}

#'
#' @title Read Growth Data section of input data file
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
  dfr = readDataframe(lns,iln,nObs,col_names=TRUE) |>
          dplyr::rowwise() |>
          dplyr::mutate(x=xs[sex]) |>
          dplyr::ungroup() |>
          dplyr::select(x,premolt_size,molt_inc,cv=CV);
  } else {

  }
  return(dfr);
}

#'
#' @title Read Growth Data section of input data file
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
    dfr1 = readDataframe(lns,iln,nEIs);
    iln = iln+nEIs;
    iln = findNextLine(lns,iln+1);
    #--determine number of rows in environmental indices dataframe
    nr = 0; for (rw in 1:nrow(dfr1)) nr = nr + (dfr1$`stop-year`-dfr1$`start-year`+1);
    dfr2 = readDataframe(lns,iln,nr);
    #--convert merged dataframe into list of individual dataframes, by env. index
    #--NOTE: possibly a dplyr function to do this in one step
    strEIs = dfr2 |> dplyr::distinct(`index/name`) |> dplyr::pull(`index/name`);
    lst = list();
    for (strEI in strEIs){
      lst[[strEI]] = dfr2 |> dplyr::filter(`index/name`==strEI);
    }
  } else {

  }
  return(lst);
}

#'
#' @title Read dataframe
#'
readDataframe<-function(lns,iln,nlns,col_names=TRUE){
  if (col_names) nlns = nlns+1;
  lnsp = stringr::str_split_i(stringr::str_trim(lns[iln:(iln+nlns-1)]),"#",1);#--remove trailing comments
  lnsp = stringr::str_replace_all(lnsp,"\\s+"," ");                           #--replace all white space with single space
  dfr = readr::read_delim(I(stringr::str_trim(lnsp)),
                        delim=" ",
                        comment="#",
                        col_names=col_names);
  return(dfr);
}

