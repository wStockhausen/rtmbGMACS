#--set up param values for allometry----
##--set up environment-----
dirPrj = rstudioapi::getActiveProject();
require(Matrix);
require(reformulas);
require(rtmbGMACS);

##--set up model configuration----
source(file.path(dirPrj,"testing/newTest_FilesCTL-Allometry/r01_SetupDims.R"));

#--rtmbGMACS allometry specifications----
#>-PROCESS_ALLOMETRY
  #--OPTIONS
  # pwrLaw1: w(z|pA,pB)       = pA*z^pB
  # pwrLaw2: w(z|pLnA,pB,pZ0) = exp(pLnA+pB*log(z/pZ0))
  #--parameter relationships:
  ##--pLnA = log(pA)  + pB*log(pZ0)
  ##--pA   = exp(pLnA - pB*log(pZ0))
  #>-CODE
    ###--probably don't need to keep original sparse indices with any derived model frames
    mfYXM  = dmsYSC |> dplyr::select(sp_idx=sparse_idx,y,x,m) |> dplyr::distinct() |>
                dplyr::mutate(yn=as.numeric(y),
                              yblk=ifelse(yn<2019,"<2019",
                                          ifelse(yn<2022,"<2022",">=2022")));
    contr_none  = "contr.none";
    contr_sum   = list(x="contr.sum",m="contr.sum");
  #>-END_CODE
  #>-PARAM_EQS
    # par_id  par_idx   frame  link_fcn  formaula             contrasts    rePDF
    # pLnA       1      mfYXM   ident    ~0+x:m+diag(1|y)     contr_none   normal  #--fixed effects: sex x maturity state, REs: intercept by year
    # pB         1      mfYXM   ident    ~0+x:m               contr_none   NA      #
    # pZ0        1      mfYXM   ident    ~0+x:m               contr_none   NA
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
txt = makeFileTemplate_Param("pLnA_1",pi_pLnA);
txt = makeFileTemplate_Param("pB_1",pi_pLnA);
txt = makeFileTemplate_Param("pZ0_1",pi_pLnA);
  #>-END_PARAM_EQS
