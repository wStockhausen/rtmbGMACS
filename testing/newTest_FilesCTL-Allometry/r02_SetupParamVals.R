#--set up param values for allometry----
##--set up environment-----
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"testing/newTest_FilesCTL-Allometry");
require(Matrix);
require(reformulas);
require(rtmbGMACS);

##--set up model configuration based on censored data----
###--this creates dfrData and dfrLims
source(file.path(dirThs,"r01_SetupDims.R"));

#--rtmbGMACS allometry specifications----
#>-PROCESS_ALLOMETRY
  #--OPTIONS
  # pwrLaw1: w(z|pA,pB)       = pA*z^pB
  # pwrLaw2: w(z|pLnA,pB,pZ0) = exp(pLnA+pB*log(z/pZ0))
  #--parameter relationships:
  ##--pLnA = log(pA)  + pB*log(pZ0)
  ##--pA   = exp(pLnA - pB*log(pZ0))
  #>-CODE
    mfYXM  = dfrData |> dplyr::select(y,x,m) |> dplyr::distinct();
    contr_none  = "contr.none";
    contr_sum   = list(x="contr.sum",m="contr.sum");
  #>-END_CODE
  #>-PARAM_EQS
    # par_id  par_idx   frame  link_fcn  formaula             contrasts    rePDF
    # pLnA       1      mfYXM   ident    ~0+x:m+diag(1|y)     contr_none   normal  #--REs: intercept by year
    # pB         1      mfYXM   ident    ~0+x:m               contr_none   NA      #--FEs only
    # pZ0        1      mfYXM   ident    ~1                   contr_none   NA      #--really just a constant: FE intercept only, will turn off estimation)
pi_pLnA = createParamInfo(formula=~0+x:m+diag(1|y),
                            model_frame=mfYXM,
                            contrasts = list(x=contr_none,m=contr_none),
                            drop.unused.levels=TRUE,
                            offset = NULL, sparse=TRUE);
pi_pB = createParamInfo(formula=~0+x:m,
                            model_frame=mfYXM,
                            contrasts = list(x=contr_none,m=contr_none),
                            drop.unused.levels=TRUE,
                            offset = NULL, sparse=TRUE);
pi_pZ0 = createParamInfo(formula=~1,
                          model_frame=mfYXM,
                          contrasts = NULL,
                          drop.unused.levels=TRUE,
                          offset = NULL, sparse=TRUE);
txt1 = makeFileTemplate_Param("pLnA_1",pi_pLnA,tab_size=4,debug=FALSE);
txt2 = makeFileTemplate_Param("pB_1",  pi_pB,  tab_size=4,debug=FALSE);
txt3 = makeFileTemplate_Param("pZ0_1", pi_pZ0, tab_size=4,debug=FALSE);
if (FALSE){ #--print to console
  cat(txt1,"\n",txt2,"\n",txt3,sep="\n");
} else {   #--print to file
  cat(txt1,"\n",txt2,"\n",txt3,file=file.path(dirThs,"txt_SetupParamVals.txt"));
}
  #>-END_PARAM_EQS
rm(txt1,txt2,txt3);
