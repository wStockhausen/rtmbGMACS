#--test createParamsInfo functions
dirPrj = rstudioapi::getActiveProject();
require(Matrix);
require(reformulas);

#--source files
source(file.path(dirPrj,"R/MiscConstants_enums.R"));
source(file.path(dirPrj,"R/MiscFunctions.R"));
source(file.path(dirPrj,"R/MiscFunctions_Formulas.R"));
source(file.path(dirPrj,"R/MiscFunctions_Matrices.R"));
source(file.path(dirPrj,"R/MiscFunctions_CovQs.R"));
source(file.path(dirPrj,"R/MiscFunctions_FEs.R"));
source(file.path(dirPrj,"R/MiscFunctions_REs.R"));
source(file.path(dirPrj,"R/createParamInfo_FEs.R"));
source(file.path(dirPrj,"R/createParamInfo_REs.R"));
source(file.path(dirPrj,"R/createParamInfo.R"));
source(file.path(dirPrj,"R/makeFileTemplates.R"));

#--setup model dimensions
source(file.path(dirPrj,"testing/newTest_ParamsInfo/r00_SetupDims.R"));

#--
contr_none  = "contr.none";
contr_sum   = list(x="contr.sum",m="contr.sum");
##--IMPORTANT NOTE:----
###--need to keep original sparse indices with any derived model frames
mfYXPa  = dmsYSC |> dplyr::select(sp_idx=sparse_idx,y,x,m,p) |> dplyr::distinct() |>
            dplyr::mutate(yn=as.numeric(y),
                          yblk=ifelse(yn<2019,"<2019",
                                      ifelse(yn<2022,"<2022",">=2022")));
mfYXPx  = dmsYSC |> dplyr::filter(x=="male")   |> dplyr::select(sp_idx=sparse_idx,y,x,p) |> dplyr::distinct() |>
            dplyr::mutate(yn=as.numeric(y),
                          yblk=ifelse(yn<2019,"<2019",
                                      ifelse(yn<2022,"<2022",">=2022")));
mfYXMf  = dmsYSC |> dplyr::filter(x=="female") |> dplyr::select(sp_idx=sparse_idx,y,x,m) |> dplyr::distinct() |>
            dplyr::mutate(yn=as.numeric(y));

f_a = ~0 + x:m +   #--fixed effects
      (0+x:m|y) +  #--REs: level of x:m for each level of y
      (1|p:yblk) + #--REs: intercept w random value for each level of p:yblk interactions ()
      (1+p||y);    #--REs: uncorrelated random intercept and level of p for each level of y

pi_pLnA_a = createParamInfo(formula=f_a,
                            model_frame=mfYXPa,
                            contrasts = list(yblk=contr_none,x=contr_none,m=contr_none),
                            drop.unused.levels=TRUE,
                            offset = NULL, sparse=TRUE);
txt = makeFileTemplate_Param("pLnA",pi_pLnA_a);

