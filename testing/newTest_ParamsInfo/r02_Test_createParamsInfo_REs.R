#--test createParamsInfo_REs functions
dirPrj = rstudioapi::getActiveProject();
require(Matrix);
require(reformulas);
require(RTMB);

#--source files
source(file.path(dirPrj,"R/MiscConstants_enums.R"));
source(file.path(dirPrj,"R/MiscFunctions.R"));
source(file.path(dirPrj,"R/MiscFunctions_Matrices.R"));
source(file.path(dirPrj,"R/MiscFunctions_Formulas.R"));
source(file.path(dirPrj,"R/MiscFunctions_CovQs.R"));
source(file.path(dirPrj,"R/MiscFunctions_REs.R"));
source(file.path(dirPrj,"R/createParamsInfo_REs.R"));
source(file.path(dirPrj,"R/createParamsInfo_FEs.R"));

#--setup model dimensions
source(file.path(dirPrj,"testing/newTest_ParamsInfo/r00_SetupDims.R"));

#--
contr_none  = "contr.none";
contr_sum   = list(x="contr.sum",m="contr.sum");
##--IMPORTANT NOTE:----
###--need to keep original sparse indices with any derived model frames
mfYXPx  = dmsYSC |> dplyr::filter(x=="male")   |> dplyr::select(sp_idx=sparse_idx,y,x,p) |> dplyr::distinct() |>
            dplyr::mutate(yn=as.numeric(y),
                          yblk=ifelse(yn<2019,"<2019",
                                      ifelse(yn<2022,"<2022",">=2022")));
mfYXMf  = dmsYSC |> dplyr::filter(x=="female") |> dplyr::select(sp_idx=sparse_idx,y,x,m) |> dplyr::distinct() |>
            dplyr::mutate(yn=as.numeric(y));

##--some fake formulas to test code
f_m = ~0 + #--no intercept
      (0+p|y) +    #--level of p for each level of y
      (1|p:yblk) + #--intercept w random value for each level of p:yblk interactions ()
      (1+p||y);    #--uncorrelated random intercept and level of p for each level of y

f_f = ~0 + m +       #--fixed effect for each level of m
       ar1(yn|yblk); #--ar1 process for yn (numeric y) w/in each level of yblk

#--just test the REs creation here (need to finish this function; works to creating reTerms)
pi_pLnA_m = createParamsInfo_REs(formula=f_m,
                                 model_frame=mfYXPx,
                                 contrasts=list(yblk=contr_none),
                                 sparse=TRUE);
pi_tmplt = makeFileTemplate_ParamsRE(formula=f_m,
                                 model_frame=mfYXPx,
                                 contrasts=list(yblk=contr_none));
#--just test the FEs creation here
pi_pLnA_f = createParamsInfo_FEs(formula=f_f,
                                 model_frame=mfYXMf,
                                 contrasts=list(m=contr_none),
                                 sparse=TRUE);
View(as.matrix(pi_pLnA_f$X));
View(pi_pLnA_f$model_matrix);


