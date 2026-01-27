#--test createParamsInfo_FEs functions
dirPrj = rstudioapi::getActiveProject();
require(reformulas);

#--source files
source(file.path(dirPrj,"R/MiscConstants_enums.R"));
source(file.path(dirPrj,"R/MiscFunctions.R"));
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

f_m = ~0 + yblk;
f_f = ~0 + m;

pi_pLnA_m = createParamsInfo_FEs(formula=f_m,
                                 model_frame=mfYXPx,
                                 contrasts = list(yblk=contr_none),
                                 offset = NULL, sparse=TRUE);
View(as.matrix(pi_pLnA_m$X));
View(pi_pLnA_m$model_matrix);

pi_pLnA_f = createParamsInfo_FEs(formula=f_f,
                                 model_frame=mfYXMf,
                                 contrasts = list(m=contr_none),
                                 offset = NULL, sparse=TRUE);
View(as.matrix(pi_pLnA_f$X));
View(pi_pLnA_f$model_matrix);


