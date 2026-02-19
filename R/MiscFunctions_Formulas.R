#' @title: Create "better" labels for parameter values.
#' @param fe_form - formula or term from formula
#' @param colnms - column names from the associated design matrix
#' @return dataframe with original (`orig`) and revised (`revd`) names.
#' @details
#' The original parameter names (`colnms`) from the design matrix for those columns involving factors are a combination of
#' the factor name and a factor level, with interactions indicated by a colon (":").
#' The revised labels for these add a "@" between a factor name and its level; interactions are still indicated with a ":".
#'
#' @importFrom rlang is_formula
#' @importFrom stringr str_detect str_extract str_split_1 str_remove
#' @importFrom tibble tibble
#' @export
#'
createParamValLabels<-function(form,colnms,debug=TRUE){
  if (!rlang::is_formula(form)){
    form = as.formula(paste0("~",deparse(form)));
  }
  trms = terms(form);
  vars= rownames(attr(trms,"factors"));
  rgx = paste0("^",vars,collapse="|");
  dfrNames <- tibble::tibble(orig=colnms,revd="");
  extract <- function(o){
      if (stringr::str_detect(o,"\\(")){
        r = o;
      } else {
        m = stringr::str_extract(o,rgx)
        r = stringr::str_remove(o,rgx);
        r = paste0(m,"@",r);
      }
    return(r);
  }
  for (i in 1:nrow(dfrNames)){
    o = stringr::str_split_1(dfrNames$orig[i],stringr::fixed(":"));
    if (debug) cat("o: '",paste(o,collapse=" "),"'\n",sep="");
    if (length(o)==1) {
      dfrNames$revd[i] = extract(o);
    } else {
      r = "";
      for (j in 1:(length(o)-1)){
        r = paste0(r,extract(o[j]),":"); #--reinsert the ":" for interactions
      }
      j = length(o);
      r = paste0(r,extract(o[j]));
      dfrNames$revd[i] = r;
    }
  }
  return(dfrNames);
}

#'
#' @title Create the model frame from a parameter formula and base dataframe
#' @description Function to create the model frame from a parameter formula and base dataframe.
#' @param formula - RHS formula describing model for which to create the model frame
#' @param data - base dataframe from which to extract model frame
#' @return A data.frame
#' @details The formula is processed to substitute '+' for RE grouping symbols ('|' and '||')
#' and drop "specials" (see [.valid_covstruct] and [.valid_smooths]) from the formula.
#' Column names for derived columns (e.g., `cos(x)`) are a back-quoted version of the deriving formula.
#'
#' @importFrom reformulas inForm noSpecials sub_specials
#' @export
#'
createModelFrame <- function(formula,data){

    mf_call <- match.call();

    if (reformulas::inForm(formula, quote(`$`))) {
        warning("use of the ", sQuote("$"), " operator in formulas is not recommended");
    }

    environment(formula) <- parent.frame();

    ## now work on evaluating model frame
    mf_call$drop.unused.levels <- TRUE;
    mf_call[[1]] <- as.name("model.frame");
    mf_call$data <- data ## propagate ..offset modification?

    ## want the model frame to contain the union of all variables used in any of the terms
    ## combine all formulas
    f <- reformulas::noSpecials(reformulas::sub_specials(formula),
                                delete=FALSE,
                                specials = c(.valid_covstruct, .valid_smooths))
    environment(f) <- environment(formula);

    mf_call$formula <- f;
    mf <- eval(mf_call,envir=environment(f),enclos=parent.frame());

    return(mf);
}

#' @title Create a RHS formula from a call or a (list of) language object(s)
#' @param terms - a call or a (list of) language object(s)
#' @return a RHS formula
#' @examplesIf FALSE
#' # example code
#' f = ~ 1 + x + s(x, bs = "tp") + s(y, bs = "gp") + s(z)
#' bars = reformulas::findbars_x(f,target="s",default.special=NULL);
#' makeRHSFormula(bars);
#'
#' @importFrom reformulas makeOp sumTerms
#'
makeRHSFormula<-function(terms){
  if (is.list(terms))
    terms = reformulas::sumTerms(terms);
  return(as.formula(reformulas::makeOp(terms,as.symbol("~"))))
}

