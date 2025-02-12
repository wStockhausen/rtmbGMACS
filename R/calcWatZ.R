#'
#' @title Get allometry info as a list
#' @description Function to get allometry info as a list. Used in the objective function.
#' @param dims -
#' @param opts -
#' @param params - (full) parameters list
#' @param verbose - flag (T/F) to print diagnostics
#' @return list with elements
#' \itemize{
#'   \item{wAtZ - 3-d array (y,s,c) with weight-at-size by year, season, pop class}
#' }
getAllometry<-function(dims,opts,params,verbose=FALSE){
  nYs = length(dims$y); #--number of years
  nSs = length(dims$s); #--number of seasons
  nCs = nrow(dims$dmsC);#--number of population classes
  if (opts$allometry$type=="data"){
    wAtZ = AD(array(data=0,dim=c(nYs,nSs,nCs)));
    for (iY in 1:nYs){
      for (iS in 1:nSs){
        for (iC in 1:nCs){
# TODO:          wAtZ[[iY,iS,iC]] = ??;
        }
      }
    }
  }
}
#'
#' @title "Calculate" weight-at-size based on "data" inputs
#' @description Function to "calculate" weight-at-size based on "data" inputs.
#'
#' @param dims - model dimensions tibble at which to "calculate" weight-at-size
#' @param dfrDims2Pars - tibble map from model dims to parameter indices
#' @param par_vals - vector of potential parameter (i.e., weight-at-size) values
#' @return vector of weights-at-size for the given model dimensions in `dims`
#'
#' @details To use this function, values for weight-at-size (i.e., "data") have
#' been given and expanded to all relevant model dimensions. Consequently, the function
#' merely looks up (hence the quotes around "calculate" in the title and description)
#' and returns the values associated with the given `dims` tibble.
#'
#' Sample code
#' ```
#' dirPrj = rstudioapi::getActiveProject();
#' conn=file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.data.txt");
#' res = readParamInfo_Allometry(conn,verbose=FALSE);
#' res1 = extractParameters_Allometry(res,dims$dms_yrxmaz,verbose=FALSE);
#' sel_dims = (dims$dms_yrxmaz |> dplyr::filter(y==2021))[1:3,];#--just extract values at 3 dimension combinations
#' sel_pars = calcWatZ_data(sel_dims,res1$dfrDims2Pars,res1$params);
#' dplyr::bind_cols(sel_dims,par_val=sel_pars)
#' # A tibble: 3 × 8
#' #     sparse_idx y    r     x      m        a         z       par_vals
#' #        <int> <fct> <fct> <fct>  <fct>    <fct>     <fct>      <dbl>
#' # 1      11777 2021  EBS   female immature new_shell 27    0.00000605
#' # 2      11778 2021  EBS   female immature new_shell 32    0.00000976
#' # 3      11779 2021  EBS   female immature new_shell 37    0.0000147
#' ````
#' @importFrom dplyr inner_join
#' @md
#' @export
#'
calcWatZ_data<-function(dims,dfrDims2Pars,par_vals){
  #--
  idxs = (dims |> dplyr::inner_join(dfrDims2Pars))$pidx;
  return(par_vals[idxs]);
}
# #--test it
# dirPrj = rstudioapi::getActiveProject();
# conn=file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.data.txt");
# res = readParamInfo_Allometry(conn,verbose=FALSE);
# res1 = extractParameters_Allometry(res,dims$dms_yrxmaz,verbose=FALSE);
# sel_dims = (dims$dms_yrxmaz |> dplyr::filter(y==2021))[1:3,];
# sel_pars = calcWatZ_data(sel_dims,res1$dfrDims2Pars,res1$params);
# dplyr::bind_cols(sel_dims,par_vals=sel_pars)

calcWatZ_function<-function(dims,
                            lstAllomInfo,
                            mpPars,
                            opPars,
                            dvPars,
                            pecPars,
                            fecPars){
  #
}
