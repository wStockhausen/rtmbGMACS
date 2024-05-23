readCTL<-function(fn){
  if (!file.exists(fn)) stop("File ",fn," does not exist!");
  lns = readLines(con=fn);
  n = 1;
  nlns = length(lns);
  while (n<=nlns){
    ln = lns[n];
    cat("Reading line",n,"\n","\t",ln,"\n")
    if (stringr::str_starts(ln,"[:blank:]*#")){
      #--skip line
    } else if (stringr::str_starts(ln,"[:blank:]*INITIAL ABUNDANCE")){
        lst = readCTL_InitialAbundance(lns,n+1);
        n = lst$n;
        lst_ini = lst$lst;
    } else if (stringr::str_starts(ln,"[:blank:]*RECRUITMENT")){
        lst = readCTL_Recruitment(lns,n+1);
        n = lst$n;
        lst_rec = lst$lst;
    } else if (stringr::str_starts(ln,"[:blank:]*GROWTH")){
        lst = readCTL_Growth(lns,n+1);
        n = lst$n;
        lst_grw = lst$lst;
    } else if (stringr::str_starts(ln,"[:blank:]*MATURITY")){
        lst = readCTL_Maturity(lns,n+1);
        n = lst$n;
        lst_mat = lst$lst;
    } else if (stringr::str_starts(ln,"[:blank:]*NATURAL MORTALITY")){
        lst = readCTL_NaturalMortality(lns,n+1);
        n = lst$n;
        lst_nm = lst$lst;
    } else if (stringr::str_starts(ln,"[:blank:]*SELECTIVITY/RETENTION")){
        lst = readCTL_Selectivity(lns,n+1);
        n = lst$n;
        lst_sel = lst$lst;
    } else if (stringr::str_starts(ln,"[:blank:]*FISHING MORTALITY")){
        lst = readCTL_FishingMortality(lns,n+1);
        n = lst$n;
        lst_fm = lst$lst;
    } else if (stringr::str_starts(ln,"[:blank:]*SURVEY INDICES")){
        lst = readCTL_SurveyIndices(lns,n+1);
        n = lst$n;
        lst_si = lst$lst;
    }
    n = n+1;
  }
  return(0);
}

parseVal<-function(str,type=NULL){
  #--split by whitespace and drop anything after comment character
  strp1 = scan(text=str,what=character(),comment.char="#",quiet=TRUE,strip.white=TRUE);
  #--parse remaining string
  val = readr::parse_guess(strp1);
  return(val);
}

parseParamInfo<-function(str,dims="",hasLabel=TRUE){
  #--split by whitespace and drop anything after comment character
  # strp1 = scan(text=str,what=character(),comment.char="#",quiet=TRUE,strip.white=TRUE);
  # csv   = paste(strp1,collapse=",");
  nms  = c("id",dims,"scale","units","ival","lb","ub","phz","type","p1","p2");
  if (hasLabel) nms = c(nms,"label");
  res = readr::read_table(I(paste0(str,"\n")),col_names=nms,comment="#")
  return(res);
}

parseTable<-function(txt){

}

readCTL_InitialAbundance<-function(lns,n){
  type = "INITIAL ABUNDANCE";
  cat("reading",type,"section starting on line",n-1,"\n");
  nlns = length(lns);
  lst = list();
  while (n<=nlns){
    ln = lns[n];
    if (stringr::str_starts(ln,"[[:blank:]*#]")){
      #--skip line
    } else if (stringr::str_starts(ln,paste0("[:blank:]*END ",type))){
      #--end of section
      cat("section ended on line",n,"\n")
      return(list(n=n,lst=lst));
    } else {
      #--read section
      optInitAbund = parseVal(ln);
      if (optInitAbund==2){
        #--read R0 parameter specification
      } else if (optInitAbund==3){
        #--read Rini parameter specification
      } else if (optInitAbund==4){
        #--read free parameters specification
      }
    }
    n = n+1;
  }
}

readCTL_Recruitment<-function(lns,n){
  type = "RECRUITMENT";
  cat("reading",type,"section starting on line",n-1,"\n");
  nlns = length(lns);
  lst = list();
  while (n<=nlns){
    ln = lns[n];
    if (stringr::str_starts(ln,"[[:blank:]*#]")){
      #--skip line
    } else if (stringr::str_starts(ln,paste0("[:blank:]*END ",type))){
      #--end of section
      cat("section ended on line",n,"\n")
      return(list(n=n,lst=lst));
    } else {
      #--read section
    }
    n = n+1;
  }
}

readCTL_Growth<-function(lns,n){
  type = "GROWTH";
  cat("reading",type,"section starting on line",n-1,"\n");
  nlns = length(lns);
  lst = list();
  while (n<=nlns){
    ln = lns[n];
    if (stringr::str_starts(ln,"[[:blank:]*#]")){
      #--skip line
    } else if (stringr::str_starts(ln,paste0("[:blank:]*END ",type))){
      #--end of section
      cat("section ended on line",n,"\n")
      return(list(n=n,lst=lst));
    } else {
      #--read section
    }
    n = n+1;
  }
}

readCTL_Maturity<-function(lns,n){
  type = "MATURITY";
  cat("reading",type,"section starting on line",n-1,"\n");
  nlns = length(lns);
  lst = list();
  while (n<=nlns){
    ln = lns[n];
    if (stringr::str_starts(ln,"[[:blank:]*#]")){
      #--skip line
    } else if (stringr::str_starts(ln,paste0("[:blank:]*END ",type))){
      #--end of section
      cat("section ended on line",n,"\n")
      return(list(n=n,lst=lst));
    } else {
      #--read section
    }
    n = n+1;
  }
}

readCTL_NaturalMortality<-function(lns,n){
  type = "NATURAL MORTALITY";
  cat("reading",type,"section starting on line",n-1,"\n");
  nlns = length(lns);
  lst = list();
  while (n<=nlns){
    ln = lns[n];
    if (stringr::str_starts(ln,"[[:blank:]*#]")){
      #--skip line
    } else if (stringr::str_starts(ln,paste0("[:blank:]*END ",type))){
      #--end of section
      cat("section ended on line",n,"\n")
      return(list(n=n,lst=lst));
    } else {
      #--read section
    }
    n = n+1;
  }
}

readCTL_Selectivity<-function(lns,n){
  type = "SELECTIVITY/RETENTION";
  cat("reading",type,"section starting on line",n-1,"\n");
  nlns = length(lns);
  lst = list();
  while (n<=nlns){
    ln = lns[n];
    if (stringr::str_starts(ln,"[[:blank:]*#]")){
      #--skip line
    } else if (stringr::str_starts(ln,paste0("[:blank:]*END ",type))){
      #--end of section
      cat("section ended on line",n,"\n")
      return(list(n=n,lst=lst));
    } else {
      #--read section
    }
    n = n+1;
  }
}

readCTL_FishingMortality<-function(lns,n){
  type = "FISHING MORTALITY";
  cat("reading",type,"section starting on line",n-1,"\n");
  nlns = length(lns);
  lst = list();
  while (n<=nlns){
    ln = lns[n];
    if (stringr::str_starts(ln,"[[:blank:]*#]")){
      #--skip line
    } else if (stringr::str_starts(ln,paste0("[:blank:]*END ",type))){
      #--end of section
      cat("section ended on line",n,"\n")
      return(list(n=n,lst=lst));
    } else {
      #--read section
    }
    n = n+1;
  }
}

readCTL_SurveyIndices<-function(lns,n){
  type = "SURVEY INDICES";
  cat("reading",type,"section starting on line",n-1,"\n");
  nlns = length(lns);
  lst = list();
  while (n<=nlns){
    ln = lns[n];
    if (stringr::str_starts(ln,"[[:blank:]*#]")){
      #--skip line
    } else if (stringr::str_starts(ln,paste0("[:blank:]*END ",type))){
      #--end of section
      cat("section ended on line",n,"\n")
      return(list(n=n,lst=lst));
    } else {
      #--read section
    }
    n = n+1;
  }
}
