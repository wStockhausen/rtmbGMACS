#' @title Get parameter information for fixed effects, random effects, and smooths in a formula
#' @description Function to get parameter information for fixed effects, random effects, and smooths in a formula.
#' @param formula - a formula from which to extract parameter info
#' @param model_frame - the model frame
#' @param contrasts - a list of contrasts for fixed effects
#' @param offset - the name of the column in `model_frame` to be used as a fixed offset (or NULL)
#' @param sparse - (logical) return sparse model matrices?
#'
#' @details The
#' @return a list composed of
#' \itemize{
#'  \item{parameter}{the name of the model parameter}
#'  \item{formula}{the input formula}
#'  \item{model_frame - the input model frame}
#'  \item{contrasts - list of contrasts used to determine model matrices}
#'  \item{drop.unused.levels - flag (T/F) to drop unused levels from design matrices}
#'  \item{offset - the name of the model_frame column used as a fixed offset}
#'  \item{piFEs}{parameter info list for fixed effects from [createParamsInfo_FEs()]}
#'  \item{piREs}{parameter info list for random effects from [createParamsInfo_REs()]}
#'  \item{piSMs}{parameter info list for smooths from [createParamsInfo_SMs()]}
#' }
#'
#' @export
#'
createParamInfo<-function(formula,
                          model_frame,
                          contrasts = NULL,
                          drop.unused.levels=FALSE,
                          offset = NULL,
                          sparse=TRUE) {

  piFEs = createParamInfo_FEs(formula,
                              model_frame,
                              contrasts=contrasts,
                              drop.unused.levels=drop.unused.levels,
                              offset=offset,
                              sparse=sparse);
  piREs = createParamInfo_REs(formula,
                              model_frame,
                              contrasts=contrasts,
                              drop.unused.levels=drop.unused.levels,
                              sparse=sparse);
  piSMs = NULL; #--TODO: create the function createParamInfo_SMs(...)
  return(namedList(formula,            #--input formula
                   model_frame,        #--input model frame
                   contrasts,          #--list of contrasts
                   drop.unused.levels, #--flag that missing levels were dropped from design matrices
                   offset,             #--name of column in `model_frame` used as a fixed offset
                   piFEs,              #--parameter info list for fixed effects
                   piREs,              #--parameter info list for random effects
                   piSMs               #--parameter info list for smooths
          ));
} #--CreateParamInfo

