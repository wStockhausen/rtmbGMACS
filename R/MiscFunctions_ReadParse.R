#'
#' @title Remove all instances of a pattern from a character vector
#' @description Function to remove all instances of a pattern from a character vector.
#' @param str - character vector to search
#' @param pattern - regex of pattern to remove
#' @return character vector with pattern removed
#' @importFrom stringr str_remove_all
#'
rmChars<-function(str,pattern=""){
  strp = stringr::str_remove_all(str,pattern=pattern)
  return(strp);
}

#'
#' @title Advance (increment) the global counter `iln` by 1
#' @description Function to advance (increment) the global counter `iln` by 1.
#' @return NULL.
#' @details On return, `iln` in the global environment will have been advanced by 1
#'
advance<-function(){
  iln<<-iln+1;
  invisible(NULL)
}

#'
#' @title Check if a string matches a key word
#' @description Function to check if a string matches a key word.
#' @param str - string to check
#' @param kws - vector of key words to check against
#' @return TRUE or FALSE.
#' @details uses [stringr::str_trim()] to trim starting/ending whitespace.
#' @importFrom stringr str_trim
#'
isKeyWord<-function(str,kws = c("list","vector","matrix","dataframe")){
  return(stringr::str_trim(str) %in% kws)
}

#'
#' @title Find the line index corresponding to a keyword in a vector of text lines
#' @description Function to find the line index corresponding to a keyword in a vector of text lines.
#' @param lns - vector of text to check
#' @param kw - key word to check
#' @return index of `lns` at first occurrence of `kw`, or NULL if not found
#' @importFrom stringr str_starts
#' @importFrom stringr str_trim
#'
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

#'
#' @title Find the next uncommented line
#' @description Function to find the next uncommented line in a vector of character strings.
#' @param lns - character vector to search
#' @param iln - index at which to start search
#' @return index into `lns` at which the next uncommented line occurs
#' @import stringr
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
#' @title Parse a text vector to create a dataframe
#' @param lns - character vector of text lines to parse
#' @param iln - index to line at which to start parsing
#' @param nlns - number of lines to parse
#' @param col_names - flag (T/F) or column names (see [readr::read_delim()])
#' @param show_col_types - flag (T/F) to show column types (see [readr::read_delim()])
#' @return a tibble (see [tibble::tibble])
#'
#' @details Uses [readr::read_delim()] to parse the vector of text lines after
#' removing all trailing comments and replacing all white space with a
#' single space.
#'
#' @importFrom readr read_delim
#' @import stringr
#'
parseTextToDataframe<-function(lns,iln,nlns,col_names=TRUE,show_col_types=FALSE){
  if (col_names) nlns = nlns+1;
  lnsp = stringr::str_split_i(stringr::str_trim(lns[iln:(iln+nlns-1)]),"#",1);#--remove trailing comments
  lnsp = stringr::str_replace_all(lnsp,"\\s+"," ");                           #--replace all white space with single space
  dfr = readr::read_delim(I(stringr::str_trim(lnsp)),
                        delim=" ",
                        comment="#",
                        col_names=col_names,
                        show_col_types=show_col_types);
  return(dfr);
}

#'
#' @title Parse a string from an input data file to a value
#' @description Function to parse a string from an input data file to a value.
#' @param str - string to parse
#' @param type - expected value type (not used currently)
#' @return a value with mode as determined by [readr::parse_guess()]
#' @details Splits the string by whitespace and drops anything after the comment character,
#' using [scan()] with `comment.char`="#" and `strip.white`=TRUE followed by [readr::parse_guess()].
#' @importFrom readr parse_guess
#'
parseVal<-function(str,type=NULL){
  #--split by whitespace and drop anything after comment character
  strp1 = scan(text=str,what=character(),comment.char="#",quiet=TRUE,strip.white=TRUE);
  #--parse remaining string
  val = readr::parse_guess(strp1);
  return(val);
}

#'
#' @title Parse multiple text lines into a list with named elements
#' @description Function to parse multiple text lines into a list with named elements.
#' @param lns - the text lines parse
#' @param names - the names of the elements of the output list
#' @return named list of parsed values
#' @details The returned list has element names matching the `names` provided.
#' Each non-white space line that doesn't start with a comment character ('#') will
#' result in an element of the returned list. Vector-valued elements should have
#' all values on the same line. Matrices and other two-dimensional (or higher)
#' constructs cannot be paresd corretly using this function.
#'
#' IMPORTANT: A counter with the name `iln` must be defined in a parent environment
#' to this function. The parsing of lines in `lns` uses the value of iln
#' at the time the function is called as the starting index in `lns`.
#' When the function exits, the value of `iln` indicates the index of the last element
#' examined.
#'
#' @importFrom stringr str_split_1
#'
parseLines<-function(lns,names){
  lst = list();
  for (i in 1:length(names)) {
    iln<<-findNextLine(lns,iln);#--iln exists in a parent environment
    lst[[i]] = parseVal(stringr::str_split_1(lns[iln],"#")[1]); advance();
  }
  names(lst) = names;
  return(lst);
}

#'
#' @title Parse multiple text lines into a list with named elements
#' @description Function to parse multiple text lines into a list with named elements until reaching an ">EOL<" marker.
#' @param lns - the text lines parse
#' @param names - the names of the elements of the output list
#' @return named list of parsed values
#' @details Each element in the text to be parsed is identified as to type
#' using a line that starts with "name: type", where type is one of "vector", "matrix",
#' or "dataframe".
#'
#' IMPORTANT: A counter with the name `iln` should be defined in the global environment
#' while this function runs. The parsing of lines in `lns` uses the value of iln
#' at the time the function is called as the starting index in `lns`.
#' When the function exits, the value of `iln` indicates the index of the last element
#' examined.
#'
#' @importFrom stringr str_split_1
#'
parseList<-function(lns,verbose=FALSE){
  out = list();
  ln = stringr::str_trim(lns[iln]);#--iln must be in the global environment
  doBreak=FALSE;
  while(!stringr::str_starts(stringr::str_trim(ln),">EOL<")){
    while((nchar(ln)==0)||stringr::str_starts(ln,"#")){
      advance();
      ln = stringr::str_trim(lns[iln]);
      if (verbose) cat("in parseList:",iln,ln,"\n")
      if (stringr::str_starts(ln,">EOL<")) doBreak=TRUE;
    }
    if (doBreak) break;
    splt = stringr::str_trim(stringr::str_split_1(ln,":"));
    if (isKeyWord(splt[2])){
      objname = stringr::str_remove(splt[1],"--");
      advance();
      if (splt[2]=="list"){
        out[[objname]] = parseList(lns);
      } else if (splt[2]=="vector") {
        out[[objname]] = parseVector(lns);
      } else if (splt[2]=="matrix") {
        out[[objname]] = parseMatrix(lns);
      } else if (splt[2]=="dataframe") {
        out[[objname]] = parseDataframe(lns);
      }
    } else {
      out[[splt[1]]] = parseVal(splt[2]);
    }
    advance();
    ln = stringr::str_trim(lns[iln]);
    if (verbose) cat("in parseList: ",iln,ln,"\n")
  }
  return(out);
}

#'
#' @title Parse multiple text lines into a dataframe
#' @description Function to parse multiple text lines into a dataframe until reaching an ">EOD<" marker.
#' @param lns - the text lines parse
#' @param verbose - flag (T/F) to print debugging info
#' @return dataframe (a [tibble::tibble()]) parsed values
#' @details
#'
#' IMPORTANT: A counter with the name `iln` should be defined in the global environment
#' while this function runs. The parsing of lines in `lns` uses the value of iln
#' at the time the function is called as the starting index in `lns`.
#' When the function exits, the value of `iln` indicates the index of the last element
#' examined.
#'
#' @importFrom stringr str_split_1
#'
parseDataframe<-function(lns,verbose=FALSE){
  out = list();
  ln = stringr::str_trim(lns[iln]);
  #--set up tibble template
  cnms = stringr::str_trim(as.character(stringr::str_split_1(ln,"[:blank:]+")));
  if (any(cnms=="")){
    ids = which(cnms=="");
    for (i in ids) cnms[i] = paste0("V",i);
  }
  ncs = length(cnms);
  lst = list();
  for (i in 1:ncs) lst[[cnms[i]]] = cnms[i];
  rw = tibble::as_tibble(lst); #--the template for a single row, all columns character-valued
  #--process dataframe
  advance();
  ln = stringr::str_trim(lns[iln]);
  rw_ctr = 1;
  doBreak=FALSE;
  while(!stringr::str_starts(stringr::str_trim(ln),">EOD<")){
    while((nchar(ln)==0)||stringr::str_starts(ln,"#")){
      advance();
      ln = stringr::str_trim(lns[iln]);
      if (verbose) cat("in parseDataframe",iln,ln,"\n")
      if (stringr::str_starts(ln,">EOD<")) doBreak=TRUE;
    }
    if (doBreak) break;
    splt = stringr::str_trim(as.character(stringr::str_split_1(ln,"[:blank:]+")));
    rw_new = rw;
    for (i in 1:length(splt)) rw_new[1,i] = splt[i];
    if (length(splt) < ncs)
      for (i in (length(splt)+1):ncs) rw_new[1,i] = NA;
    out[[rw_ctr]] = rw_new;
    rw_ctr = rw_ctr+1;
    advance();
    ln = stringr::str_trim(lns[iln]);
    if (verbose) cat("in parseDataframe: ",iln,ln,"\n")
  }
  dfr = dplyr::bind_rows(out);
  #--TODO: change to numeric values where appropriate
  return(dfr);
}

#'
#' @title Parse multiple text lines into a matrix
#' @description Function to parse multiple text lines into a matrix until reaching an ">EOM<" marker.
#' @param lns - the text lines parse
#' @param verbose - flag (T/F) to print debugging info
#' @return matrix of parsed values
#' @details
#'
#' IMPORTANT: A counter with the name `iln` should be defined in the global environment
#' while this function runs. The parsing of lines in `lns` uses the value of iln
#' at the time the function is called as the starting index in `lns`.
#' When the function exits, the value of `iln` indicates the index of the last element
#' examined.
#'
#' @importFrom stringr str_split_1
#'
parseMatrix<-function(lns,verbose=FALSE){
  out = list();
  ln = stringr::str_trim(lns[iln]);
  #--strategy here is to parse a dataframe, then convert to matrix
  #--set up tibble template
  vls  = stringr::str_trim(as.character(stringr::str_split_1(ln,"[:blank:]+")));
  ncls = length(vls);
  cnms = rep("V",ncls);#--column names do not exist, otherwise should be a dataframe
  for (i in 1:ncls) cnms[i] = paste0("V",i);
  lst = list();
  for (i in 1:ncls) lst[[cnms[i]]] = vls[i];
  rw = tibble::as_tibble(lst); #--the template for a single row, all columns character-valued
  rw_ctr = 1;
  out[[rw_ctr]] = rw;
  #--process remaining rows
  rw_ctr = 2;
  advance();
  ln = stringr::str_trim(lns[iln]);
  doBreak=FALSE;
  while(!stringr::str_starts(stringr::str_trim(ln),">EOM<")){
    while((nchar(ln)==0)||stringr::str_starts(ln,"#")){
      advance();
      ln = stringr::str_trim(lns[iln]);
      if (verbose) cat("in parseMatrix:",iln,ln,"\n")
      if (stringr::str_starts(ln,">EOM<")) doBreak=TRUE;
    }
    if (doBreak) break;
    splt = stringr::str_trim(as.character(stringr::str_split_1(ln,"[:blank:]+")));
    rw_new = rw;
    for (i in 1:length(splt)) rw_new[1,i] = splt[i];
    if (length(splt) < ncls)
      for (i in (length(splt)+1):ncls) rw_new[1,i] = NA;
    out[[rw_ctr]] = rw_new;
    rw_ctr = rw_ctr+1;
    advance();
    ln = stringr::str_trim(lns[iln]);
    if (verbose) cat("in parseMatrix",iln,ln,"\n");
  }
  dfr = dplyr::bind_rows(out);
  #--convert to numeric matrix
  nr =nrow(dfr);
  nc = ncol(dfr);
  mt = matrix(as.numeric(as.matrix(dfr)),
              nrow=nr,ncol=nc,byrow=FALSE);
  return(mt);
}

#'
#' @title Parse a 1-element text vector into a vector
#' @description Function to parse a 1-element text vector into a vector.
#' @param lns - the text lines parse
#' @param verbose - flag (T/F) to print debugging info
#' @return vector of parsed values
#' @details
#'
#' IMPORTANT: A counter with the name `iln` should be defined in the global environment
#' while this function runs. The line `lns[iln]` is parsed to create the vector.
#' The value of `iln` is not incremented (because there was no search).
#'
#' @importFrom stringr str_trim
#'
parseVector<-function(lns,verbose=FALSE){
  ln=stringr::str_trim(lns[iln]);
  if (verbose) cat("in parseVector: ",iln,ln,"\n");
  vls = parseVal(ln);
  return(vls);
}
