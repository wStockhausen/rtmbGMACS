#'
#' @title Extract character vector from section of longer vector
#' @description Function to extract character vector from section of longer vector.
#' @param txt - long character vector
#' @param start - string to regard as marking the start of section (not included in the section)
#' @param end - string to regard as marking the end of the section (not included in the section)
#' @return a character vector of the "enclosed" section
#' @details The text lines at which the `start` and `end` are found are not included in
#' the "section" returned by this function.
#' @import stringr
#' @md
#' @export
#'
extractLines<-function(txt,start,end){
  s = stringr::str_which(txt,stringr::regex(paste0("^",start),ignore_case=TRUE));
  e = stringr::str_which(txt,stringr::regex(paste0("^",end),ignore_case=TRUE));
  if ((s+1)<=(e-1)) return(txt[(s+1):(e-1)]);
  return(vector(mode="character"));
}

removeCommentLines<-function(strv,comment="#"){
  return(strv[stringr::str_starts(stringr::str_trim(strv),comment,negate=TRUE)]);
}

#' @title Identify first non-comment line in a character vector
#' @description Function to identify the index of the first non-comment line in a character vector.
#' @param txt - character vector to search
#' @param start - index into `txt` to regard as the start of the search
#' @param comment - character indicating the following text in the string is a comment
#' @return the index of the first line in `txt` that does not start with a comment character
#' @details Each element in `txt` is regarded as a line of text (as in a text file).
#' The function searches for the first non-comment line and returns its index.
#' @examples
#' txt = c("#this is a comment",
#'         "  #so is this,even thought the line starts with blank space",
#'         "this is not a comment",
#'         "this is not either, but the first is line 4.");
#' extractLines(txt,1); #--should be 3
#' @import stringr
#' @md
#' @export
#'
skipCommentLines<-function(txt,start,comment="#"){
  idx = start;
  while(stringr::str_starts(txt[idx],comment)){idx<-idx+1;}
  return(idx);
}

#' @title Extract a text section with only non-comment lines from a character vector
#' @description Function to extract a text section with only non-comment lines from a character vector.
#' @param txt - character vector to search
#' @param n - number of non-comment lines to extract (default=length(txt))
#' @param start - index into `txt` to regard as the start of the search (default=1)
#' @param comment - character indicating the following text in the string is a comment
#' @return list with elements `txt` and `end` (see details)
#' @details Each element in `txt` is regarded as a line of text (as in a text file).
#' The function searches for the first `n` non-comment lines and extracts them as a new vector.
#' The returned list has elements
#' \itemize{
#'   \item{txt - character vector of length `n` with comment lines removed (i.e., a subset of `txt`)}
#'   \item{end - index of the last element extracted in the input `txt` vector}
#' }
#' @importFrom stringr str_starts
#' @md
#' @export
extractTextSection<-function(txt,n=length(txt),start=1,comment="#"){
  #--extract a text section with `n` non-commented lines
  idx = 0; idxp = 0;
  itxt = vector("integer",length=length(txt));#--set to max length
  while ((idxp<n)&&(idxp<=length(txt))){
    while(stringr::str_starts(txt[start+idx+idxp],comment)) idx = idx+1;
    itxt[idxp+1] = start+idx+idxp; idxp = idxp + 1;
  }
  txtp = txt[itxt[1:idxp]];
  #--discard trailing comments
  idx=1;
  while(idx<=idxp){
    txtp[idx] = stringr::str_trim(stringr::str_split_1(txtp[idx],"#"))[1];
    idx = idx+1;
  }
  return(list(txt=txtp,end=itxt[idxp]));
}

#' @title Extract a list from a character vector
#' @description Function to extract a list from a character vector.
#' @param strv - character vector to parse
#' @param split - character(s) to use to split lines into name and value (default="<-")
#' @param verbose - flag to print diagnostic info
#' @return list with named elements
#' @details Each element in `strv` is an equation defining the
#' elements of a dimension.
#'
#' @examplesIf FALSE
#' # example code
#' str=paste(
#'   'MODEL_DIMS
#'   y <- 2020:2024;                           #--years
#'   s <- 1;                                   #--seasons
#'   r <- "EBS";                               #--regions
#'   x <- c("male","female");                  #--sex classes
#'   m <- c("immature","mature");              #--maturity state classes
#'   p <- c("new_shell","old_shell");          #--post-molt ages
#'   zc <- seq(55.5,104.5,5);                  #--size bin cutpoints
#'   f <- c("TCF","SCF","NMFS");               #--fleets
#'   END');
#' strv = stringr::str_split_1(str,"\\n") |> extractLines("MODEL_DIMS","END");
#' lstDims = parseStrAsList(strv);
#' @import stringr
#' @md
#' @export
parseStrAsList<-function(strv,split="<-",verbose=FALSE){
  strp = extractTextSection(strv)$txt;
  ns   = length(strp);
  lst = list();
  i = 1;
  while(i<=ns){
    strpp = strp[i] |> stringr::str_remove_all(" ") |> stringr::str_split_1(split);
    i = i+1;
    test = !any(strpp |> stringr::str_ends(";"));
    while((i<=ns)&&test){
      strpp = c(strpp,
                strp[i] |> stringr::str_remove_all(" ") |> stringr::str_split_1(split));
      i = i+1;
      test = !any(strpp |> stringr::str_ends(";"));
      if (verbose) cat(test,strpp,"\n")
    }
    if (verbose) cat(strpp,"\n")
    lst[[strpp[1]]] = eval(parse(text=strpp[2:length(strpp)]),envir=lst);
  }
  return(lst);
}

#' @title Extract a dataframe (a [tibble::tibble()]) from a character vector
#' @description Function to a dataframe (a [tibble::tibble()]) from a character vector.
#' @param txt - character vector to search
#' @param n - number of non-comment lines to extract
#' @param start - index into `txt` to regard as the start of the search
#' @param comment - character indicating the following text in the string is a comment
#' @return list with elements `txt` and `end` (see details)
#' @details Each element in `txt` is regarded as a line of text (as in a text file).
#' The dataframe structure and column types are determined by [readr::read_table()].
#' \itemize{
#'   \item{txt - character vector of length `n` with comment lines removed (i.e., a subset of `txt`)}
#'   \item{end - index of the last element extracted in the input `txt` vector}
#' }
#' @importFrom readr read_table
#' @md
#' @export
extractDataframe<-function(txt,n,start,comment="#"){
  idx = skipCommentLines(txt,start);
  dfr = readr::read_table(txt[idx:(idx+n)],col_names=TRUE,comment=comment,skip_empty_rows=TRUE);
  return(dfr);
}

#' @title Parse strings in a character vector to digits
#' @description Function to parse strings in a character vector to digits.
#' @param x - character vector to parse
#' @return integer vector
#' @details For each element in `x`, everything except for digits is removed and
#' the resulting character value is converted to an integer.
#' @examples
#' parseToDigits(c("ab23","cd42","a42 b23"))
#'
#' @importFrom stringr str_remove_all
#' @md
#' @export
parseToDigits<-function(x){
  str = stringr::str_remove_all(x,stringr::regex("[^[:digit:]]",ignore_case=TRUE));
  return(as.integer(str))
}

#' @title Extract text between first and last parentheses
#' @description Function to extract text between first and last parentheses
#' @param txt - character vector from which to extract strings
#' @return character vector of extracted strings
#' @details For each element in the vector, this extracts the string between the
#' *first* "(" and *last* ")" from each element of the input
#' character vector, or `NA` if no matching parentheses.
#'
#' @examples
#' # example code
#' extractTextBetweenParens("fred(x,y)=0");
#' extractTextBetweenParens(c("fred(x,y)=0","b(z0)"));
#' extractTextBetweenParens("fred(x,y(z))=0");
#'
#' @importFrom stringr str_extract
#' @export
#'
extractTextBetweenParens<-function(txt){
  return(stringr::str_extract(txt,"(?<=\\().*(?=\\))"));
}

#' @title Extract text between first and last brackets
#' @description Function to extract text between first and last brackets
#' @param txt - character vector from which to extract strings
#' @return character vector of extracted strings
#' @details For each element in the vector, this extracts the string between the
#' *first* "\[" and *last* "\]" from each element of the input
#' character vector, or `NA` if no matching brackets.
#'
#' @examples
#' # example code
#' extractTextBetweenBrackets("fred[x,y]=0");
#' extractTextBetweenParens(c("fred[x,y]=0","b[z0]"));
#' extractTextBetweenParens("fred[x,y[z]]=0");
#'
#' @importFrom stringr str_extract
#' @export
#'
extractTextBetweenBrackets<-function(txt){
  return(stringr::str_extract(txt,"(?<=\\[).*(?=\\])"));
}

#' @title Extract text after a string
#' @description Function to extract text after a string
#' @param txt - character vector from which to extract strings
#' @param str - string to extract after (see NOTE under details)
#' @return character vector of extracted strings
#' @details For each element in the vector, this extracts the string following
#' the `after` string from each element of the input
#' character vector, or `NA` if no `after` string found.
#'
#' NOTE: Characters with special regex meanings (e.g., ")","]") need to
#' be double-escaped to use in `str`.
#'
#' @examples
#' # example code
#' extractTextAfterString("test:fred[x,y]==0","==");
#' extractTextAfterString(c("fred[x,y]==0: FALSE","b[z0]: 1"),":");
#' extractTextAfterString("fred[x,y[z]]=0","=");
#' extractTextAfterString("fred[x,y[z]]=0","\\[");
#' extractTextAfterString("fred[x,y[z]]=0","fred\\[");
#'
#' @importFrom stringr str_extract
#' @export
#'
extractTextAfterString<-function(txt,str){
  return(stringr::str_extract(txt,paste0("(?<=",str,").*")));
}

#' @title Concatenate text
#' @description Function to concatenate text
#' @param ... - text values to concatenate
#' @param concat - string to use as concatenator
#' @return string
#' @details The elements in ... are concateated into a single string using the
#' `concat` string separated by spaces.
#'
#' @examples
#' # example code
#' concatText(1,2,3)
#' concatText(1,NA,3)
#' concatText(1,2,NULL,4)
#' export
#'
concatText<-function(...,concat="+"){
  n = ...length();
  str = ...elt(1);
  ctr = 2;
  while (ctr<=n){
    if (!is.null(...elt(ctr))) str = paste(str,concat,...elt(ctr));
    ctr = ctr+1;
  }
  return(str);
}

