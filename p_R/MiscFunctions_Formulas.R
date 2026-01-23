#'
#'@title Extract `mgcv` "smooths" from text defining a formula
#'@description Function to extract `mgcv` "smooths" from text defining a formula.
#'@param txt - length 1 character vector with text describing a formula
#'@return character vector with smooth functions
#'
#'@details This function is used when processing a formula (as text) to extract any
#'`mgcv` smooth functions `s` and `ti` from the text. `mgcv` smooth functions must be
#'handled separately in terms of creating a model matrix from a formula and a model frame.
#'
#'@examplesIf FALSE
#'# example code
#'txt = "~1 + yblk + s(z,k=10) + s(m,k=20,by='yblk') +ti(z)";
#'res = extractSmoothsFromFormulaText(txt); #--res == c("s(z,k=10)","s(m,k=20,by='yblk')","ti(z)")
#'
#'@importFrom stringr fixed
#'@importFrom stringr str_remove
#'@importFrom stringr str_remove_all
#'@export
#'
extractSmoothsFromFormulaText<-function(txt){
  #--see https://www.spsanderson.com/steveondata/posts/2024-06-25/ for regexp
  res = stringr::str_extract_all(txt,"s\\(.*?\\)",simplify=TRUE)[1,];
  res = c(res,stringr::str_extract_all(txt,"ti\\(.*?\\)",simplify=TRUE)[1,]);
  return(res);
}

#'
#'@title Remove `mgcv` "smooths" from text defining a formula
#'@description Function to remove `mgcv` "smooths" from text defining a formula.
#'@param txt - length 1 character vector with resulting text describing a formula
#'@return length 1 character vector with smooth functions removed
#'
#'@details This function is used when processing a formula (as text) to remove any
#'`mgcv` smooth functions `s` and `ti` from the text. `mgcv` smooth functions must be
#'handled separately in terms of creating a model matrix from a formula and a model frame.
#'
#'@examplesIf FALSE
#'# example code
#'txt = "~1 + yblk + s(z,k=10) + s(m,k=20,by='yblk') +ti(z)";
#'res = removeSmoothsFromFormulaText(txt); #--res == "~1+yblk";
#'
#'@importFrom stringr fixed
#'@importFrom stringr str_remove
#'@importFrom stringr str_remove_all
#'@export
#'
removeSmoothsFromFormulaText<-function(txt){
  txtp = txt;
  smths = extractSmoothsFromFormulaText(txtp);
  if (length(smths)>0)
    for (i in 1:length(smths))
      txtp = stringr::str_remove(txtp,stringr::fixed(smths[i]));
  txtp = stringr::str_remove_all(txtp," ");   #--remove any white space
  txtp = stringr::str_remove(txtp,"^(\\+)+"); #--remove any forward "+"'s
  txtp = stringr::str_remove(txtp,"(\\+)+$"); #--remove any trailing "+"'s
  return(txtp);
}

#' @title Extract the variables associated with a 1-d "smooth" function
#' @title Function to extract the variables associated with a 1-d "smooth" function.
#' @param txt - text string defining the "smooth" function call
#' @return character vector with the name of the variables
#' @details The function text must have the form "smooth(...,k=...)", where "smooth" is the
#' `mgcv` smooth function (`s` or `ti`); the first "..." indicates the comma-delimited variables for the smooth,
#' where the first variable is numeric and any other (optional) variable must be a factor; and the second ...
#' is the number of knots followed by any other (comma-delimited) inputs to the smooth function.
#' @examplesIf FALSE
#' # example code: smooth specification with "by" variable
#' txt = "s(x,k=10,by='r')"; #--x is numeric
#' res = extractSmoothVar(txt); # res == "x"
#' # example code: factor smooths specification
#' txt = "s(x,y,k=10,bs='fs',by='r')"; #--x is numeric, y is a factor
#' res = extractSmoothVar(txt); # res == c("x","y")
#'
#' @export
extractSmoothVars<-function(txt){
  res = stringr::str_extract_all(txt,"\\(.*?,k",simplify=TRUE)[1,];
  res = stringr::str_remove(res,"\\(");
  res = stringr::str_remove(res,",k");
  res = stringr::str_split_1(res,",");
  return(res);
}

#' @title Calculate a model matrix for non-covariate fixed effects and smooths
#' @description Function to calculate a non-covariate fixed effects model matrix (and associated information)
#' @param txt - model formula as text
#' @param dfrMF - basis dataframe for model frame
#' @param ctrs - contrasts for linear fixed effects factors
#' @param verbose - flag to print debugging info
#' @return list with number of parameters (`npars`) and, if `npars`>0, the model matrix and other information.
#' @details This function calculates a model matrix based on the input formula (as text) and
#' model frame basis. The returned list has elements
#' \items{
#'   \item{mtxRedMM - the reduced model matrix (for distinct factor and variable levels)}
#'   \item{dfrRedMF - the reduced model (data)frame with distinct factor and variable levels; `row_idx` indicates corresponding row in mtxRedMM}
#'   \item{dfrMdMtx - dataframe combining MdMtx and dfrRedMF}
#'   \item{dfrMF - the basis dataframe for the model frame, with an appended column with each value indicating the corresponding row in mtxRedMM}
#'   \item{f - the input formula (as text)}
#'   \item{ns - the number of non-covariate fixed effects parameters}
#'   \item{ks - an integer vector with the number of knots for each smooth (if any)}
#'   \item{npars - total number of RTMB parameters (this includes smooths parameters)}
#'   \item{SmCons - a list of the smooth construct objects for each smooth (if any)}
#' }
#' The model matrix includes columns for both non-covariate fixed effects (if any) and smooths (if any),
#' with fixed effect terms coming first. The number of columns reflects the number of RTMB (i.e., actual)
#' parameters associated with the model matrix. The model matrix has the same number of rows as the
#' model frame and `dfrMF`. The model frame is a dataframe consisting of a subset of the columns of `dfrMF` that
#' includes only the variables included in the formula,
#' with the columns corresponding to any smooth variables converted to numeric type and possibly preceded
#' by a column of 1's with name '(Intercept)' if the formula includes an intercept.
#'
#' A version of the model matrix with only distinct rows is returned as `mtxRedMM`. `dfrRedMF` is a dataframe with
#' the rows corresponding to the (distinct) combination of variable values in the associated row in `mtxRedMM`, and
#' `dfrMdMtx` is the concatenated dataframe version of these. Given an appropriate RTMB parameter vector (`p`, say),
#' the values of the model parameter vector derived for each distinct combination of variable values in `dfrRedMF` can be
#' calculated as `mtxRedMM` %*% p. The basis dataframe (`dfrMF`) is included as part of the
#' returned list to allow
#'
#' If `txt` is not a valid formula (as a character string), a list with only `npars`= 0 is returned.
#'
#' @examplesIf FALSE
#' # example code
#' #--define mdoel DimsMap
#' vRs="EBS";
#' miz = as.character(seq(25,75,5));
#' mmz = as.character(seq(50,75,5));
#' fiz = as.character(seq(25,65,5));
#' fmz = as.character(seq(30,65,5));
#' vXs=list(male=list(imm=list(`new`=miz,
#'                              `old`=miz),
#'                     mat=list(`new`=mmz,
#'                              `old`=mmz,
#'                              `very`=mmz)
#'                   ),
#'           female=list(imm=list(`new`=fiz,
#'                                `old`=fiz),
#'                       mat=list(`new`=fmz,
#'                                `old`=fmz)
#'                      )
#'          ); attr(vXs,"dmnms")<-c("x","m","p","z");
#' ys = as.character(2015:2024); names(ys) = ys;
#' ss = as.character(1);         names(ss) = ss;
#' dmsYSC  = createSparseDimsMap(y=ys,s=ss,r=vRs,x=vXs);  #--DimsMap for combination
#' #--define time block for NMFS Q's
#' lstNMFS<-list(preCV =as.character(c(2015:2019)),       #--pre gear change
#'               preGC =as.character(c(2021:2022)),       #--pre gear change
#'               postGC=as.character(c(2023:2025)));      #--post gear change
#' tbNMFS = lstNMFS |> defineYearBlocks();
#' ##--expand to dimensions map for NMFS Q's
#' dmsNMFS<-tbNMFS |>
#'          defineYearBlocks(dimsMdl=dmsYSC) |>
#'          dplyr::mutate(f="NMFS");
#' ###--string defining a formula with fixed effects and a smooth on size (`z`)
#' txt = "~0 + yblk + s(z,k=10)";
#' ##--calculate model matrix and associated information
#' ctrs = list(yblk="contr.sum");
#' res = calcModelMatrixFEs(txt,dmsNMFS,ctrs);
#'
#'@import dplyr
#'@import Formula
#'@import glue
#'@importFrom Matrix cbind2
#'@importFrom Matrix Matrix
#'@import mgcv
#'
#'@export
#'
calcModelMatrixFEs<-function(txt,dfrMF,ctrs=NULL,verbose=FALSE){
  ###--check for NA or non-formula
  if ((is.na(txt))||(!(txt |> stringr::str_starts("~")))) return(list(npars=0));

  ###--define symbol for model matrix
  MdMtx = NULL;
  ###--define symbol for model frame
  MdFrm = NULL;
  ###--extract any fixed effects
  if (verbose) cat("txt = ",txt,"\n");
  f  = removeSmoothsFromFormulaText(txt);
  if (verbose) cat("f = ",f,"\n");
  ns = vector("integer",0);
  fs = vector("character",0);
  vs = vector("character",0);
  if (f!="~"){
    ###--create fixed factors model formula
    frmla = eval(parse(text=f));
    ###--create the model frame from the formula and the basis for the model frame
    MdFrm = model.frame(frmla,data=dfrMF);
    ###--revise formula by dropping any 1-value variables
    frmla = reviseFormula(frmla,MdFrm);
    ###--convert character columns to factors and apply contrasts to model frame
    dfr   = convertToDimsFactors(MdFrm,NULL,contrasts=ctrs,verbose=verbose);
    #--create the model matrix
    MdMtx = model.matrix(frmla,data=dfr,rhs=NULL);
    ns = ncol(MdMtx); #--number of fixed effects parameters (including covariates)
    fs = c(getColumnsByType(MdFrm,"character"),
           getColumnsByType(MdFrm,"factor"));#--names of factor columns
    vs = getColumnsByType(MdFrm,"numeric");  #--names of covariate columns
  } #--otherwise no fixed factors

  ###--extract any smooths
  sms  = extractSmoothsFromFormulaText(txt);
  ks  = vector("integer",0);
  scs = list();
  if (length(sms)>0){
    for (i in 1:length(sms)){
      sso = eval(parse(text=sms[i])); #--smooth specification formula (requires package mgcv)
      vrs = extractSmoothVars(sms[i]);
      vr = vrs[1]; #--numeric variable
      dfrTmp = dfrMF |> dplyr::mutate("{vr}":=as.numeric(!!rlang::sym(vr))); #--basis for model frame
      if (length(vrs)>1){
        for (j in 2:length(vrs)){
          vr = vrs[j];
          dfrTmp = dfrTmp |> dplyr::mutate("{vr}":=factor(!!rlang::sym(vr))); #--basis for model frame
        }
      }
      sc  = mgcv::smoothCon(sso,data=dfrTmp,knots=NULL);             #--smooth construct
      terms = sc[[1]]$term; #--terms (variables) for smooths
      byvar = sc[[1]]$by;   #--by variable
      scs[[f]] = sc;
      if (is.null(MdMtx)) {
        MdMtx  = sc[[1]]$X;   #--model matrix for (first) smooth term
        ks = ncol(MdMtx);     #--number of columns in the model matrix
        MdFrm = NULL;
        for (term in terms)
          MdFrm = dplyr::bind_cols(MdFrm,dfrTmp |> dplyr::select({{term}}));
        if (!(is.na(byvar)||(byvar=="NA"))) MdFrm = dplyr::bind_cols(MdFrm,dfrTmp |> dplyr::select({{byvar}}));
      } else {
        Xf = sc[[1]]$X;       #--model matrix for (first) smooth term
        ks = c(ks,ncol(Xf));  #--vector of number of terms in smooths
        MdMtx = Matrix::cbind2(MdMtx,Xf);   #--combined model matrix
        for (term in terms)
          if (!(term %in% names(MdFrm))) MdFrm = dplyr::bind_cols(MdFrm,dfrTmp |> dplyr::select({{term}}));
      }
    }#--i loop
  }

  #--combine model frame and model matrix, then reduce to unique combinations of elements
  ##--need to exclude covariate columns in MdFrm because they are included in MdMtx
  nfs   = length(fs);  #--number of factors
  nvs   = length(vs);  #--number of simple covariates
  nps   = ncol(MdMtx); #--length of parameter vector
  MdFrmp = tibble::as_tibble(MdFrm);
  drop   = names(MdFrmp)[(names(MdFrmp) %in% colnames(MdMtx))];
  dfrMM  = MdFrmp |> dplyr::select(!dplyr::any_of(drop)) |> dplyr::bind_cols(MdMtx) |> dplyr::distinct();
  ##--extract reduced model matrix (only model matrix columns)
  #redMM = Matrix::Matrix(as.matrix(dfrMM |> dplyr::select(!dplyr::any_of(c(fs,vs)))));
  redMM = Matrix::Matrix(as.matrix(dfrMM |> dplyr::select(!dplyr::any_of(c(fs)))));
  #--extract reduced model frame w/ corresponding reduced model matrix row index as first column
  redMF = dfrMM |> dplyr::select(dplyr::any_of(c(fs,vs))) |> tibble::rowid_to_column("row_idx");
  ##--augment basis dataframe with corresponding row index of reduced model matrix
  if (length(colnames(MdFrm))>0){
    dfrMFp = dfrMF |> dplyr::left_join(redMF,by=colnames(MdFrm));
  } else {
    dfrMFp = dfrMF |> dplyr::cross_join(redMF);
  }
  if ("(Intercept)" %in% names(dfrMFp)) {
    warning("calcModelMatrixFE: dropping intercept term from dfrMFp.")
    dfrMFp = dfrMFp |> dplyr::select(!`(Intercept)`);
  }

  facnms = fs; #--should this include covariates?
  #--get variable (covariate) names from model matrix/frame
  varnms = vs;
  if (verbose) {
    cat("facnms =",facnms,"\n");
    cat("varnms =",varnms,"\n");
  }
  ##--calculate reduced model matrix corresponding to model parameter factor levels
  colnms=colnames(redMM);
  dfrRedT = dfrMFp |>
              #dplyr::bind_cols(as.matrix(lst$mtxRedMM)) |>
              dplyr::inner_join(dfrMM,by=c(facnms,varnms)) |>
              dplyr::group_by(dplyr::pick(dplyr::any_of(facnms))) |>
              dplyr::summarize(dplyr::across(dplyr::all_of(c(colnms,varnms)),mean)) |> dplyr::ungroup();
  mtxRedT = Matrix::Matrix(as.matrix(dfrRedT |> dplyr::select(!dplyr::any_of(c(facnms)))));

  #--create output list
  res = list(mtxRedMM=mtxRedT,mtxRedMM0=redMM,dfrRedMF=redMF,dfrModMtx=dfrMM,dfrMF=dfrMFp,
             f=txt,factors=fs,variables=vs,npars=nps);
  return(res);
}

#' @title Calculate a model matrix for "environmental" covariates
#' @description Function to calculate a model matrix (and associated information) for "environmental" covariates
#' @param txt - model formula for environmental covariates as text
#' @param dfrMF - basis dataframe for model frame
#' @param ctrs - contrasts for any covariate factors
#' @param verbose - flag to print debugging info
#' @return list with number of parameters (`npars`) and, if `npars`>0, the model matrix and other information.
#' @details This function calculates a model matrix based on the input formula (as text) and
#' model frame basis. The returned list has elements
#' \items{
#'   \item{mtxRedMM - the reduced model matrix (for distinct factor and variable levels)}
#'   \item{dfrRedMF - the reduced model (data)frame with distinct factor and variable levels; `row_idx` indicates corresponding row in mtxRedMM}
#'   \item{dfrMdMtx - dataframe combining MdMtx and dfrRedMF}
#'   \item{dfrMF - the basis dataframe for the model frame, with an appended column with each value indicating the corresponding row in mtxRedMM}
#'   \item{f - the input formula (as text)}
#'   \item{ns - the number of covariate parameters}
#'   \item{npars - total number of RTMB parameters (same as `ns`)}
#' }
#' The model matrix includes columns for both fixed effects (if any) and smooths (if any),
#' with fixed effect terms coming first. The number of columns reflects the number of RTMB (i.e., actual)
#' parameters associated with the model matrix. The model matrix has the same number of rows as the
#' model frame and `dfrMF`. The model frame is a dataframe consisting of a subset of the columns of `dfrMF` that
#' includes only the variables included in the formula,
#' with the columns corresponding to any smooth variables converted to numeric type and possibly preceded
#' by a column of 1's with name '(Intercept)' if the formula includes an intercept.
#'
#' A version of the model matrix with only distinct rows is returned as `mtxRedMM`. `dfrRedMF` is a dataframe with
#' the rows corresponding to the (distinct) combination of variable values in the associated row in `mtxRedMM`, and
#' `dfrMdMtx` is the concatenated dataframe version of these. Given an appropriate RTMB parameter vector (`p`, say),
#' the values of the model parameter vector derived for each distinct combination of variable values in `dfrRedMF` can be
#' calculated as `mtxRedMM` %*% p. The basis dataframe (`dfrMF`) is included as part of the
#' returned list to allow
#'
#' If `txt` is not a valid formula (as a character string), a list with only `npars`= 0 is returned.
#'
#' @examplesIf FALSE
#' # example code
#' #--define mdoel DimsMap
#' vRs="EBS";
#' miz = as.character(seq(25,75,5));
#' mmz = as.character(seq(50,75,5));
#' fiz = as.character(seq(25,65,5));
#' fmz = as.character(seq(30,65,5));
#' vXs=list(male=list(imm=list(`new`=miz,
#'                              `old`=miz),
#'                     mat=list(`new`=mmz,
#'                              `old`=mmz,
#'                              `very`=mmz)
#'                   ),
#'           female=list(imm=list(`new`=fiz,
#'                                `old`=fiz),
#'                       mat=list(`new`=fmz,
#'                                `old`=fmz)
#'                      )
#'          ); attr(vXs,"dmnms")<-c("x","m","p","z");
#' ys = as.character(2015:2024); names(ys) = ys;
#' ss = as.character(1);         names(ss) = ss;
#' dmsYSC  = createSparseDimsMap(y=ys,s=ss,r=vRs,x=vXs);  #--DimsMap for combination
#' #--define time block for NMFS Q's
#' lstNMFS<-list(preCV =as.character(c(2015:2019)),       #--pre gear change
#'               preGC =as.character(c(2021:2022)),       #--pre gear change
#'               postGC=as.character(c(2023:2025)));      #--post gear change
#' tbNMFS = lstNMFS |> defineYearBlocks();
#' ##--expand to dimensions map for NMFS Q's
#' dmsNMFS<-tbNMFS |>
#'          defineYearBlocks(dimsMdl=dmsYSC) |>
#'          dplyr::mutate(f="NMFS");
#' ###--string defining a formula with fixed effects and a smooth on size (`z`)
#' txt = "~0 + yblk + s(z,k=10)";
#' ##--calculate model matrix and associated information
#' ctrs = list(yblk="contr.sum");
#' res = calcModelMatrixFEs(txt,dmsNMFS,ctrs);
#'
#'@import dplyr
#'@import Formula
#'@import glue
#'@importFrom Matrix cbind2
#'@importFrom Matrix Matrix
#'@import mgcv
#'
#'@export
#'
calcModelMatrixECs<-function(txt,dfrMF,ctrs=NULL,verbose=FALSE){
  ###--check for NA or non-formula
  if ((is.na(txt))||(!(txt |> stringr::str_starts("~")))) return(list(npars=0));

  ###--define symbol for model matrix
  MdMtx = NULL;
  ###--define symbol for model frame
  MdFrm = NULL;
  ###--extract covariates
  if (verbose) cat("txt = ",txt,"\n");
  f  = removeSmoothsFromFormulaText(txt); #--should not have any smooths
  if (verbose) cat("f = ",f,"\n");
  ns = vector("integer",0);
  fs = vector("character",0);
  vs = vector("character",0);
  if (f!="~"){
    ###--create fixed factors model formula
    #Frmla = Formula::Formula(eval(parse(text=f)));
    frmla = eval(parse(text=f));
    ###--create the model frame from the formula and the basis for the model frame
    #MdFrm = Formula:::model.frame.Formula(Frmla,data=dfrMF);
    MdFrm = model.frame(frmla,data=dfrMF);
    ###--revise formula by dropping any 1-value variables
    frmla = reviseFormula(frmla,MdFrm);
    ###--convert character columns to factors and apply contrasts to model frame
    dfr   = convertToDimsFactors(MdFrm,NULL,contrasts=ctrs,verbose=verbose);
    #--create the model matrix
    #MdMtx = Formula:::model.matrix.Formula(Frmla,data=dfr,rhs=NULL);
    MdMtx = model.matrix(frmla,data=dfr,rhs=NULL);
    ns = ncol(MdMtx); #--number of fixed effects parameters (including covariates)
    fs = c(getColumnsByType(MdFrm,"character"),
           getColumnsByType(MdFrm,"factor"));#--factor columns
    vs = getColumnsByType(MdFrm,"numeric")   #--covariate columns
  } #--otherwise no fixed factors

  #--combine model frame and model matrix, then reduce to unique combinations of elements
  ##--need to exclude covariate columns in MdFrm because they are included in MdMtx
  nfs   = length(fs);  #--number of factors
  nvs   = length(vs);  #--number of simple covariates
  nps   = ncol(MdMtx); #--length of parameter vector
  MdFrmp = tibble::as_tibble(MdFrm);
  drop   = names(MdFrmp)[(names(MdFrmp) %in% colnames(MdMtx))];
  dfrMM  = MdFrmp |> dplyr::select(!dplyr::any_of(drop)) |> dplyr::bind_cols(MdMtx) |> dplyr::distinct();
  ##--extract reduced model matrix (only model matrix columns)
  redMM = Matrix::Matrix(as.matrix(dfrMM |> dplyr::select(!dplyr::any_of(c(fs,vs)))));
  #--extract reduced model frame w/ corresponding reduced model matrix row index as first column
  redMF = dfrMM |> dplyr::select(dplyr::any_of(c(fs,vs))) |> tibble::rowid_to_column("row_idx");
  ##--augment basis dataframe with corresponding row index of reduced model matrix
  dfrMFp = dfrMF |> dplyr::left_join(redMF,by=colnames(MdFrm));
  if ("(Intercept)" %in% names(dfrMFp)) {
    str = paste0("The equation '",txt,", for the formula for an environmental covariate has an intercept. ",
                 "If you did something like 'v*(1+f)', where v is a numeric covariate and f is a factor, ",
                 "replace '*' with ':' (receall that 'a*b' expands the formula to 'a+b+a:b'. ",
                 "If required, the intercept should be included in formula for the fixed effects.");
    stop(str);
  }

  #--create output list
  res = list(mtxRedMM=redMM,dfrRedMF=redMF,dfrModMtx=dfrMM,dfrMF=dfrMFp,f=txt,
             factors=fs,variables=vs,npars=nps);
  return(res);
}

#TODO: fill in for random effects model matrix and associated info
#' @title Calculate a model matrix for random effects
#' @description Function to calculate a random effects model matrix (and associated information)
#' @param txt - model formula as text
#' @param dfrMF - basis dataframe for model frame
#' @param cov_type - covariance type
#' @param verbose - flag to print debugging info
#' @return list with number of parameters (`npars`) and, if `npars`>0, the model matrix and other information.
#' @details TODO: fill in.
#'@import dplyr
#'@import Formula
#'@import glue
#'@importFrom Matrix cbind2
#'@import mgcv
#'
#'@export
#'
calcModelMatrixREs<-function(txt,dfrMF,cov_type,verbose=FALSE){
  ###--check for NA or non-formula
  if ((is.na(txt))||(!(txt |> stringr::str_starts("~")))) return(list(npars=0));

  if (tolower(cov_type)=="iid"){

  } else
  if (tolower(cov_type)=="ar1"){

  } else {
    error(paste0("unrecognized `cov_type` in `calcModelMatrixREs`: '",cov_type,"'."))
  }

  ###--TODO: fill this in!
  res = list(npars=0);#--for now!
  return(res);
}

#'
#' @title Convert a formula to a text string
#' @description Function to convert a formula to a text string .
#' @param f - formula object
#' @return a text string representation of the formula
#' @details None.
#' @examplesIf FALSE
#' # example code
#' ##--one-sided formula
#' getFormulaAsText(~1+x*m+s(y));
#' ##--two-sided formula
#' getFormulaAsText(z~1+x*m+s(y));
#'
#' @export
#'
getFormulaAsText<-function(f){
  txt = as.character(f);
  if (attr(terms(f),"response")==0){
    #--one-sided formula
    txt = paste(txt,collapse="");
  } else {
    #--two-sided formula
    txt = paste(txt[2],"~",txt[3])
  }
  return(txt);
}

#'
#' @title Remove all single-value variables included in a model frame from a formula
#' @description Function to remove all single-value variables included in a model frame from a formula.
#' @param f - model formula or equivalent text string (1-element character vector)
#' @param dfrMF - basis dataframe from which to extract the model frame for `f`
#' @result the (possibly revised) formula with all 1-value variables removed
#' @details `R` formulae cannot be used with variables with only one value in the model frame due to the apparent lack
#' of a method to create 1-element contrasts when creating the model matrix (`R` throws an error when any
#' 1-value variables are included in a formula). Based on the input model frame (`dfrMF`), this function removes any
#' 1-value variables included in the formula `f` and returns the revised formula (or the original, if no variables in the model frame had only one value)
#'
#' This function throws an error if all variables in `f` have only one value and `f` does not include an intercept term.
#'
#' @import dplyr
#' @export
#'
reviseFormula<-function(f,dfrMF){
  #--need to revise formula for any 1-value columns in the model frame
  if (rlang::is_formula(f)){
    txt = getFormulaAsText(f);#--create text version of formula
  } else if (rlang::is_character(f)){
    txt = f;                  #--copy text version of formula
    f = eval(parse(text=txt));#--create formula
  }
  #--create model frame consistent with original formula (contains all variables included in f)
  mf = model.frame(f,data=dfrMF,drop.unused.levels=FALSE,xlev=NULL);#--using defaults for drop... and xlev
  #--determine distinct combinations of variable values
  dmf = mf |> dplyr::distinct();
  ##--determine number of unique values for each variable in formula
  ndmf = dmf |> dplyr::group_by() |> dplyr::summarise(dplyr::across(tidyr::everything(),~length(unique(.x))));
  ##--determine which variables have only one value (these variables cannot be included in the final formula)
  cols = colnames(ndmf);      #--variables (column names) in formula
  vars = vector("character"); #--place
  for (col in cols) if (ndmf[[col]][1]==1) vars = c(vars,col);
  ##--eliminate all 1-value variables from the model formula
  ###--NOTE: 1-value variables cause errors when evaluating a model matrix due to inability to form contrasts
  if (length(vars)==0){
    #--no 1-value variables: nothing to do
    ##--return the original formula, the equivalent text, the model frame with all variables in the original formula,
    ###--and the "distinctivised" model frame with all variables in the original formula
    return(f);
  } else {
    ###--remove 1-value variables from all formula terms
    trms = attr(terms(f),"term.labels");#--text versions of individual terms in formula
    itcp = attr(terms(f),"intercept");  #--0 if no intercept in formula, 1 otherwise
    ###--loop over 1-value variables
    for (var in vars){
      trms = trms |> stringr::str_remove_all(paste0("\\b",var,"\\b")) |>
                     stringr::str_remove("^:") |>
                     stringr::str_remove(":$");
    }
    ###--reconstitute model formula (without possible intercept) keeping only unique terms
    trms1 = paste0(unique(trms[trms!=""]),collapse="+");
    ###--add in possible intercept
    txt1 = paste0("~",itcp);
    if (trms1!="") {
      txt1  = paste(txt1,trms1,sep="+");
    } else if (itcp==0) {
      stop(paste0("Formula '",txt,"' has no intercept and only 1-value variables. Revise to include an intercept."))
    }
    ##--return the revised formula
    return(eval(parse(text=txt1)));
  }
  #--shouldn't get here
  return(NULL);
}
if (FALSE){
model.matrix(reviseFormula(~yblk,dfrMF),dfrMF);
model.matrix(reviseFormula(~yblk,dfrMF),dfrMF,contrasts.arg=list(yblk="contr.sum"));
model.matrix(reviseFormula(~yblk-1,dfrMF),dfrMF);
model.matrix(reviseFormula(~yblk-1,dfrMF),dfrMF,contrasts.arg=list(yblk="contr.sum"));
model.matrix(reviseFormula(~yn,          dfrMF),dfrMF,contrasts.arg=list(yblk="contr.treatment"));
model.matrix(reviseFormula(~yn-1,        dfrMF),dfrMF,contrasts.arg=list(yblk="contr.treatment"));
model.matrix(reviseFormula(~yn+yn:yblk,  dfrMF),dfrMF,contrasts.arg=list(yblk="contr.treatment"));
model.matrix(reviseFormula(~yn+yn:yblk-1,dfrMF),dfrMF,contrasts.arg=list(yblk="contr.treatment"));
model.matrix(reviseFormula(~yn+yn:yblk-1,dfrMF),dfrMF,contrasts.arg=list(yblk="contr.sum"));
model.matrix(reviseFormula(~   yn*yblk-1,dfrMF),dfrMF,contrasts.arg=list(yblk="contr.sum"));
model.matrix(reviseFormula(~   yn*yblk,  dfrMF),dfrMF,contrasts.arg=list(yblk="contr.sum"));
model.matrix(reviseFormula(~   yn:yblk,  dfrMF),dfrMF,contrasts.arg=list(yblk="contr.treatment"));
mm = model.matrix(reviseFormula(~   yn %in% yblk,  dfrMF),dfrMF,contrasts.arg=list(yblk="contr.sum")); mm;

mm = model.matrix(reviseFormula(~yn:(yblk*p)-1,dfrMF),dfrMF,contrasts.arg=list(yblk="contr.sum",p="contr.sum")); mm;

vIVs = c(0.0027,0.0027,0.0027,0.0027,0.0027,0.0027);
varnms = c("x","yblk","p");
parnms=colnames(mm);
redT = dfrMF |> dplyr::bind_cols(mm) |> dplyr::group_by(dplyr::pick(dplyr::any_of(varnms))) |>
         dplyr::summarize(dplyr::across(dplyr::all_of(parnms),mean)) |> dplyr::ungroup();
mtxRedT = as.matrix(redT |> dplyr::select(!dplyr::any_of(varnms)));
lst = svd(mtxRedT);
vP  = as.vector((lst$v %*% diag(1/lst$d,nrow=length(lst$d)) %*% t(lst$u)) %*% vIVs);
mtxRedT %*% vP;
}
