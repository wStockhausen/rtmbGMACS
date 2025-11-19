
#' @title Get model dimension names
#' @description Function to get model dimension names.
#' @return character vector with model dimension names
#' @export
getDimNames<-function(){
  return(c("f","y","s","r","x","m","p","z"));
}

#' @title Get which model dimensions are identified in a character vector
#' @description Function to get which model dimensions are identified in a character vector.
#' @param char_vec - character vector to check for model dimension names
#' @return character vector with identifed model dimension names
#' @export
whichDims<-function(char_vec){
  return(char_vec[char_vec %in% getDimNames()]);
}

#'
#' @title Split text string into a character vector
#' @description Function to split text string to a character vector.
#' @param txt - text string (1-elenent character vector) to split
#' @param split - string used to split text (default="\\n)
#' @return a character vector
#' @details The text is plit into a character vector using [stringr::str_split_1()].
#' @importFrom stringr str_split_1
#' @md
#' @export
#'
splitText<-function(txt,split="\\n"){
  v = stringr::str_split_1(txt,split);
  return(v);
}
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
  txt = stringr::str_trim(txt);
  s = stringr::str_which(txt,stringr::regex(paste0("^",start),ignore_case=TRUE));
  e = stringr::str_which(txt,stringr::regex(paste0("^",end),ignore_case=TRUE));
  if ((length(s)>0)&&(length(e)>0))
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

#'
#' @title Evaluate text as lines of code
#' @description Function to evaluate text as lines of code.
#' @param strv - character vector with lines of code to evaluate.
#' @return nothing
#' @details The lines of `code` in strv are evaluated in the frame specified by
#' `frame`. If `frame` is non-negative, it specifies the parent frame relative to the
#' caller (so `frame=0`, the deault, is evaluated in the caller's frame). If `frame` is negative,
#' the lines of code are evaluated in the environment specified by its absolute value.
#' @examplesIf  FALSE
#' str = "
#'   v1  = 5;
#'   v2 = 1:v1;
#'   sum = list(x="contr.sum",m="contr.sum");
#'       ";
#' # example code
#'
#' @export
evalTextAsCode<-function(strv,frame=0){
  if (is.environment(frame)){
    eval(parse(text=strv),envir=frame);
  } else if (is.numeric(frame)){
    if (frame>=0){
      eval.parent(parse(text=strv),n=frame+2);#--caller frame is up 2
    } else {
      eval(parse(text=strv),envir=abs(frame));
    }
  } else if (is.character(frame)){
    eval(parse(text=strv),envir=frame);
  }
}

#' @title Evaluate a list from a character vector
#' @description Function to evaluate a list from a character vector.
#' @param strv - character vector to parse/evaluate
#' @param split - character(s) to use to split lines into name and value (default="<-")
#' @param verbose - flag to print diagnostic info
#' @return list with named elements
#' @details Each quantity assigned using `<-` in `strv` is interpreted as an equation defining an element of the returned list.
#' Quantities assigned using `=` outside a list structure are evaluated locally for use in constructing a list element of
#' the returned list. Attributes for a list element of the returned list must be assigned (see example code 2) immediately after
#' the list element.
#'
#' @examplesIf FALSE
#' # example code 1
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
#' strv = str |> splitText() |> extractLines("MODEL_DIMS","END");
#' lstDims = evalTextAsList(strv);
#' # example code 2
#'   str=paste(
#'   'MODEL_DIMS
#'     r <- "EBS";
#'     miz = as.character(seq(25,75,5)); #--size bin midpoints for immature males
#'     mmz = as.character(seq(50,75,5)); #--size bin midpoints for mature males
#'     fiz = as.character(seq(25,65,5)); #--size bin midpoints for immature females
#'     fmz = as.character(seq(30,65,5)); #--size bin midpoints for mature females
#'     x<-list(male=list(imm=list(`new`=miz,
#'                                  `old`=miz),
#'                         mat=list(`new`=mmz,
#'                                  `old`=mmz,
#'                                  `very`=mmz)
#'                       ),
#'               female=list(imm=list(`new`=fiz,
#'                                    `old`=fiz),
#'                           mat=list(`new`=fmz,
#'                                    `old`=fmz)
#'                          )
#'              );
#'     attr(x,"dmnms") = c("x","m","p","z");#--define order of dimensions in nested list
#'   END'
#'   )
#'   strv = str |> splitText() |> extractLines("MODEL_DIMS","END");
#'   lstDms = evalTextAsList(strv);
#' @import stringr
#' @md
#' @export
evalTextAsList<-function(strv,split="<-",verbose=FALSE){
  strp = extractTextSection(strv)$txt;
  ns   = length(strp);
  lst = list();
  i = 1;
  while(i<=ns){
    strpp = strp[i] |> stringr::str_remove_all(" ") |> stringr::str_split_1(split) |> stringr::str_trim();
    i = i+1;
    test = !any(strpp |> stringr::str_ends(";"));
    while((i<=ns)&&test){
      strpp = c(strpp,
                strp[i] |> stringr::str_remove_all(" ") |> stringr::str_split_1(split));
      i = i+1;
      test = !any(strpp |> stringr::str_ends(";"));
      if (verbose) cat(test,strpp,"\n")
    }
    if (verbose) cat(strpp,"\n");
    if (length(strpp)==1){
      if (stringr::str_starts(strpp,"attr")){
        #--apply attribute to "last" element in lst
        txt = stringr::str_replace(strpp,last,paste0("lst[['",last,"']]"));
        if (verbose) cat("Applying attributes:\n\t",txt,"\n");
        eval(parse(text=txt));
      } else {
        eval(parse(text=strpp));
      }
    } else {
      eval(parse(text=paste0(strpp[1],"=",paste0(strpp[2:length(strpp)],collapse=""))));
      lst[[strpp[1]]] = eval(parse(text=paste0(strpp[2:length(strpp)],collapse="")));
      last=strpp[1];
    }
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
#'
#' This function differs from [evalTextAsDataframe()] in that
#' the input character vector may be a much larger section of text than what defines the
#' dataframe to be extracted, but the number of rows (`n`) to extract is known, as is the index of
#' the element to begin with (`start`). In [evalTextAsDataframe()], after comments are removed are
#' removed from the input character vector `char_vec`, the resulting
#' character vector is regarded as defining the dataframe.
#'
#' @importFrom readr read_table
#' @md
#' @export
extractDataframe<-function(txt,n,start,comment="#"){
  idx = skipCommentLines(txt,start,comment=comment);
  dfr = readr::read_table(txt[idx:(idx+n)],col_names=TRUE,comment=comment,skip_empty_rows=TRUE);
  return(dfr);
}

#' @title Evaluate a character vector as a dataframe (a `tbl_df`)
#' @description Function to evaluate a character vector as a dataframe (a `tbl_df`).
#' @param char_vec - character vector to parse/evaluate as a dataframe
#' @param comment - string indicating that what follows is a comment
#' @return a dataframe (actually a `tbl_df`, see [tibble::tibble()])
#' @details Each element in `char_vec` is regarded as a line of text (as in a text file).
#' The dataframe structure and column types are determined by [readr::read_table()].
#'
#' This function differs from [extractDataframe()] in that, after comments are removed are
#' removed from the input character vector `char_vec`, the resulting
#' character vector is regarded as defining the dataframe. In [extractDataframe()],
#' the input character vector may be a much larger section of text than what defines the
#' dataframe to be extracted, but the number of rows (`n`) to extract is known as, is the index of
#' the element to begin with (`start`).
#'
#' @examplesIf FALSE
#' # example code
#' str<-paste0('
#'   id        function  frame   params              description
#'   pwrLaw1   pwrLaw1   mfALL   pA,pB               w(z)_=_pA_*_(z^pB)
#'   pwrLaw2   pwrLaw2   mf2024  pLnA,pLnS,pB,pZ0    w(z)_=_exp(pLnA+pLnS+pB*ln(z/pZ0))
#' ')
#' dfr = txt |> splitText() |> evalTextAsDataframe();
#'
#' @importFrom readr read_table
#' @md
#' @export
evalTextAsDataframe<-function(char_vec,comment="#"){
  txtp = (char_vec |> extractTextSection(comment=comment))$txt;
  dfr = readr::read_table(txtp,col_names=TRUE,comment=comment,skip_empty_rows=TRUE);
  return(dfr);
}

#' @title Evaluate strings in a character vector as digits
#' @description Function to evaluate strings in a character vector as digits.
#' @param x - character vector to parse/evaluate
#' @return integer vector
#' @details For each element in `x`, everything except for digits is removed and
#' the resulting character value is converted to an integer.
#' @examples
#' evalTextAsDigits(c("ab23","cd42","a42 b23"))
#'
#' @importFrom stringr str_remove_all
#' @md
#' @export
evalTextAsDigits<-function(x){
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

