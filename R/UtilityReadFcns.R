rmChars<-function(str,pattern=""){
  strp = stringr::str_remove_all(str,pattern=pattern)
  return(strp);
}

advance<-function(){
  iln<<-iln+1;
  invisible(NULL)
}

isKeyWord<-function(str){
  kws = c("list","vector","matrix","dataframe");
  return(stringr::str_trim(str) %in% kws)
}

parseList<-function(lns,verbose=FALSE){
  out = list();
  ln = stringr::str_trim(lns[iln]);
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

parseVector<-function(lns,verbose=FALSE){
  #advance();
  ln=stringr::str_trim(lns[iln]);
  if (verbose) cat("in parseVector: ",iln,ln,"\n");
  vls = parseVal(ln);
  return(vls);
}
