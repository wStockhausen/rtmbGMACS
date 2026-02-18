#'
#' @title Create the parameters information list for a set of RE terms
#' @description Function to create the parameters information list for a set of RE terms.
#' @param formula formula, possibly containing random effects
#' @param model_frame full model frame
#' @param contrasts a list of contrasts
#' @param sparse (logical) return sparse model matrix?
#' @return a list, or NULL (see details)
#' @details
#'
#' If `formula` does not contain any RE terms, NULL is returned. Otherwise, a list is returned with the following elements:
#' \itemize{
#'  \item{reTerms}{RE terms as character strings}
#'  \item{Z}{model matrix for the REs}
#'  \item{re_pvs_idx}{a named vector with a cumulative index for the RE parameter values}
#'  \item{re_pvs_lbls}{a list if labels for the RE parameter values}
#'  \item{lstQtp}{a list of term-associated precision matrix templates (see [mkReTermInfoList()])}
#'  \item{lstQfill}{a list of term-associated functions to fill the associated precision matrix (see [mkReTermInfoList()])}
#' }
#'
#' @export
#'
createParamsInfo_REs<-function(formula,
                                 model_frame,
                                 contrasts=NULL,
                                 drop.unused.levels=TRUE,
                                 sparse=TRUE) {
  #--get number of "observation" (rows) in model frame
  nobs <- nrow(model_frame);
  #--extract only RE parts of a formula
  re_form <- getREs(formula);
  if (is.null(re_form)){
    return(NULL);
  }

  #--split RE formula into component RE terms
  ##--returns a list with elements reTrms, reTrmFormulas, reTrmAddArgs, reTrmCovTypes
  sp_form <- splitForm_RE(re_form);

  #--create a list of output from mkReTermInfoList for each RE term
  reTermsInfoList = lapply(sp_form$reTrms, mkReTermInfoList, model_frame,
                           drop.unused.levels, sparse = sparse);

  #--combine info (TODO: modify `model_frame` for output somehow?)
  reTerms <- names(reTermsInfoList);
  Ztlist  <- lapply(reTermsInfoList, `[[`, "termZt"); #--extract termZt's into a list
  Zt      <- do.call(rbind, Ztlist)                   #--eq. 7, JSS lmer paper
  Z       <- t(Zt);
  q_re         <- ncol(Z); #--total number of RE parameter values (number of columns in Z)
  re_pvs_idx   <- cumsum(sapply(reTermsInfoList,`[[`,"npvs")); #--cumulative index for RE parameter values
  re_pvs_lbls  <- lapply(reTermsInfoList,`[[`,"pv_lbls");      #--parameter value labels
  if (q_re != rev(re_pvs_idx)[1]) {
    stop(paste0("{q_re = ",q_re,"} != {last(re_pvs_idx) = ",rev(re_pvs_idx)[1],"}\n"));
  }
  if (rev(re_pvs_idx)[1] != length(unlist(re_pvs_lbls))) {
    stop(paste0("{last(re_pvs_idx) = ",rev(re_pvs_idx)[1],"} != {length(re_pvs_lbls) = ",length(re_pvs_lbls),"}\n"));
  }

  covtype  = lapply(reTermsInfoList, `[[`, "covstr"); #--list of covariance types
  lstQtp   = lapply(reTermsInfoList, `[[`, "Qtp");    #--list of term-associated precision matrix templates
  lstQfill = lapply(reTermsInfoList, `[[`, "fillQ");  #--list of term-associated functions to fill the precision matrices

  return(list(reTerms=reTerms,          #--RE terms as character strings
              model_frame=model_frame,  #--model frame (TODO: want to relate this to Z)
              Z=Z,                      #--RE model matrix (see eq. 2 in JSS lmer paper)
              re_pvs_idx=re_pvs_idx,    #--named vector with cumulative index for RE parameter values
              re_pvs_lbls=re_pvs_lbls,  #--list of labels for RE parameter values
              covtype = covtype,        #--list of covariance types
              lstQtp=lstQtp,            #--list of term-associated precision matrix templates
              lstQfill=lstQfill         #--list of term-associated functions to fill the precision matrices
             ));

}#--createParamsInfo_REs
