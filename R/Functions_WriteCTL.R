#--write functions for CTL file
writeVal<-function(val){
  str = "NA";
  if (!is.null(val))
    str = paste(val);
  return(str);
}
writeCTL_Table<-function(lst){
  str = "###-------------------------------------------------------------------------------------";
  if (!is.null(lst)){
    strp = readr::format_delim(lst,delim=" ",quote="none");
    strp = stringr::str_replace(strp,"\\n$","");
    str = c(str,strp,
                "END TABLE");
  }
  str = c(str,"###-------------------------------------------------------------------------------------");
  return(str);
}

writeCTL<-function(lst,conn,justText=TRUE){
  if (!justText)
    if(!isOpen(conn)) conn = file(conn,open="w");
  str = c("#--tmbGMACS CTL file",
          "#--model population dynamics",
          "################################################################################"
          );
  str = c(str,writeCTL_SeasonalIntegration(lst[["SeasonalIntegration"]],conn,justText=TRUE));
  str = c(str,"################################################################################");
  str = c(str,writeCTL_InitialAbundance(lst[["InitAbd"]],conn,justText=TRUE));
  str = c(str,"################################################################################");
  str = c(str,writeCTL_Allometry(lst[["allom"]],conn,justText=TRUE));
  str = c(str,"################################################################################");
  str = c(str,writeCTL_Recruitment(lst[["rec"]],conn,justText=TRUE));
  str = c(str,"################################################################################");
  str = c(str,writeCTL_NaturalMortality(lst[["nm"]],conn,justText=TRUE));
  str = c(str,"################################################################################");
  str = c(str,writeCTL_Growth(lst[["grw"]],conn,justText=TRUE));
  str = c(str,"################################################################################");
  str = c(str,writeCTL_Maturity(lst[["mat"]],conn,justText=TRUE));
  str = c(str,"################################################################################");
  str = c(str,writeCTL_Selectivity(lst[["sel"]],conn,justText=TRUE));
  str = c(str,"################################################################################");
  str = c(str,writeCTL_FishingMortality(lst[["fsh"]],conn,justText=TRUE));
  str = c(str,"################################################################################");
  str = c(str,writeCTL_SurveyIndices(lst[["srv"]],conn,justText=TRUE));

  str = c(str,"END CTL FILE");
  if (!justText) cat(str,sep="\n",file=conn);
  return(paste(str,collapse="\n"));
}

writeCTL_SeasonalIntegration<-function(lst,conn,justText=TRUE){
  if (!justText)
    if(!isOpen(conn)) conn = file(conn,open="w");

  str = c("SEASONAL INTEGRATION--##########################################################",
          "##--values (fraction of year) for category 'duration' must be defined for",
          "##----each year and season and sum to 1 for each year",
          "##",
          ,"#              duration (fraction of year)",
          writeCTL_Table(lst[["seasonal_integration"]])
          );
  str = c(str,"END SEASONAL INTEGRATION SECTION","\n");

  if (!justText) cat(str,sep="\n",file=conn);
  return(paste(str,collapse="\n"));
}

writeCTL_InitialAbundance<-function(lst,conn,justText=TRUE){
  if (!justText)
    if(!isOpen(conn)) conn = file(conn,open="w");

  str = c("INITIAL ABUNDANCE--##########################################################",
          "##--options",
          "###--1. no initial population (build up over time from recruitment)",
          "###--2. based on unfished equilibrium and recruitment (R0)",
          "###--3. based on equilibrium with initial year fishing mortality and recruitment (Rini)",
          "###--4. free parameters (relative to reference class(es))",
          "##--",
          paste(writeVal(lst$option),"    #--initial abundance option"),
          "###--Option 1: no parameters",
          "###--Option 2: parameters (R0)",
          writeCTL_Table(lst[["R0_tbl"]]),
          "###--Option 3: parameters (Rini)",
          writeCTL_Table(lst[["Rini_tbl"]]),
          "###--Option 4: free parameters",
          "###--reference classes",
          writeCTL_Table(lst[["ref_pars_tbl"]]),
          "###--nonparametric initial abundance offset options",
          "####--1. initial values as abundance in class",
          "####--2. initial values as ln-scale offsets from specified class",
          "####--3. initial values as ln-scale deviations",
          "###--",
          paste(writeVal(lst$opt4_vals),"     #--initial value scaling option"),
          "###--initial abundance offsets",
          writeCTL_Table(lst[["off_pars_tbl"]])
          );
  str = c(str,"END INITIAL ABUNDANCE SECTION","\n");

  if (!justText) cat(str,sep="\n",file=conn);
  return(paste(str,collapse="\n"));
}

writeCTL_Allometry<-function(lst,conn,justText=TRUE){
  if (!justText)
    if(!isOpen(conn)) conn = file(conn,open="w");

  str = c("ALLOMETRY--##########################################################");

  str = c(str,"##--weight-at-size function types",
              "###--1. PL: power-law (alpha*Z^beta)",
              "###--2. NP: nonparametric",
              writeCTL_Table(lst[["fcns"]]),
              " ",
              "##--weight-at-size reference parameters",
              writeCTL_Table(lst[["ref_pars"]]),
              " ",
              "##--weight-at-size non-parametric reference parameters",
              writeCTL_Table(lst[["ref_nonpars"]]));

  str = c(str,"END ALLOMETRY SECTION","\n");

  if (!justText) cat(str,sep="\n",file=conn);
  return(paste(str,collapse="\n"));
}

writeCTL_Recruitment<-function(lst,conn,justText=TRUE){
  if (!justText)
    if(!isOpen(conn)) conn = file(conn,open="w");

  str = c("RECRUITMENT--##########################################################");
  str = c(str,
          "##--size-at-recruitment function types",
          "###--1. gamma (alpha, beta parameters)",
          "###--2. np    (nonparametric)",
          "##--annual rec types",
          "###--1. noSR1   (no stock-recruit relationship, ln-scale offsets)",
          "###--2. noSR-RW (no stock-recruit relationship, ln-scale RW offsets)",
          "##--sex ratio function types",
          "###--1. constant: estimated constant",
          " ",
          "##--recruitment functions",
           writeCTL_Table(lst[["rec_fcns"]]),
          " ",
          "##--total recruitment reference parameters",
           writeCTL_Table(lst[["rec_ref_pars"]]),
          "##--total recruitment offset parameters",
           writeCTL_Table(lst[["rec_off_pars"]]),
          "##--total recruitment parameter covariates",
           writeCTL_Table(lst[["rec_covars"]]),
          " ",
          "##--size function reference parameters",
           writeCTL_Table(lst[["rec-z_ref_pars"]]),
          "##--size function offset parameters",
           writeCTL_Table(lst[["rec-z_off_pars"]]),
          "##--size function parameter covariates",
           writeCTL_Table(lst[["rec-z_covars"]]),
          " ",
          "##--sex ratio reference parameters",
           writeCTL_Table(lst[["sexrat_ref_pars"]]),
          "##--sex ratio offset parameters",
           writeCTL_Table(lst[["sexrat_off_pars"]]),
          "##--total recruitment parameter covariates",
           writeCTL_Table(lst[["sexrat_covars"]]));

  str = c(str,"END RECRUITMENT SECTION","\n");

  if (!justText) cat(str,sep="\n",file=conn);
  return(paste(str,collapse="\n"));
}

writeCTL_NaturalMortality<-function(lst,conn,justText=TRUE){
  if (!justText)
    if(!isOpen(conn)) conn = file(conn,open="w");

  str = c("NATURAL MORTALITY--##########################################################",
          "##--lnM = lnM_ref + lnM_off + lnM_cov + ln(f(z))",
          "##--size function (f(z)) types",
          "###--1. constant (1)",
          "###--2. Lorenzen (1/z)",
          writeCTL_Table(lst[["fcns"]]),
          " ",
          "##--natural mortality reference parameters",
          writeCTL_Table(lst[["ref_pars"]]),
          "##--natural mortality offset parameters",
          writeCTL_Table(lst[["off_pars"]]),
          "##--natural mortality parameter covariates",
          writeCTL_Table(lst[["covars"]]));
  str = c(str,"END NATURAL MORTALITY SECTION","\n");

  if (!justText) cat(str,sep="\n",file=conn);
  return(paste(str,collapse="\n"));
}

writeCTL_Growth<-function(lst,conn,justText=TRUE){
  if (!justText)
    if(!isOpen(conn)) conn = file(conn,open="w");

  str = c("GROWTH--##########################################################");
  str = c(str,"END GROWTH SECTION","\n");

  if (!justText) cat(str,sep="\n",file=conn);
  return(paste(str,collapse="\n"));
}

writeCTL_Maturity<-function(lst,conn,justText=TRUE){
  if (!justText)
    if(!isOpen(conn)) conn = file(conn,open="w");

  str = c("MATURITY--##########################################################");
  str = c(str,"END MATURITY SECTION","\n");

  if (!justText) cat(str,sep="\n",file=conn);
  return(paste(str,collapse="\n"));
}

writeCTL_Selectivity<-function(lst,conn,justText=TRUE){
  if (!justText)
    if(!isOpen(conn)) conn = file(conn,open="w");

  str = c("SELECTIVITY--##########################################################");
  str = c(str,"END SELECTIVITY SECTION","\n");

  if (!justText) cat(str,sep="\n",file=conn);
  return(paste(str,collapse="\n"));
}

writeCTL_FishingMortality<-function(lst,conn,justText=TRUE){
  if (!justText)
    if(!isOpen(conn)) conn = file(conn,open="w");

  str = c("FISHING MORTALITY--##########################################################");
  str = c(str,"END FISHING MORTALITY SECTION","\n");

  if (!justText) cat(str,sep="\n",file=conn);
  return(paste(str,collapse="\n"));
}

writeCTL_SurveyIndices<-function(lst,conn,justText=TRUE){
  if (!justText)
    if(!isOpen(conn)) conn = file(conn,open="w");

  str = c("SURVEY INDICES--##########################################################");
  str = c(str,"END SURVEY INDICES SECTION","\n");

  if (!justText) cat(str,sep="\n",file=conn);
  return(paste(str,collapse="\n"));
}

#--for testing
#cat(writeCTL(tst,"",justText=TRUE),sep="\n");
