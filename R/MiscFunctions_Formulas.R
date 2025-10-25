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

#' @title Calculate a model matrix for fixed effects
#' @description Function to calculate a fixed effects model matrix (and associated information)
#' @param txt - model formula as text
#' @param dfrMF - basis dataframe for model frame
#' @param ctrs - contrasts for linear fixed effects factors
#' @param verbose - flag to print debugging info
#' @return list with model matrix and other information, or NULL if not `txt` is not a formula
#' @details This function calculates a model matrix based on the input formula (as text) and
#' model frame basis. The returned list has elements
#' \items{
#'   \item{MdMtx - the reduced model matrix (for distinct factor and variable levels)}
#'   \item{dfrMdMtx - dataframe with distinct factor and variable levels and reduced model matrix}
#'   \item{fullMdMtx - the full model matrix}
#'   \item{ModFrm - the full model frame}
#'   \item{f - the input formula (as text)}
#'   \item{ns - the number of fixed effects parameters}
#'   \item{ks - an integer vector with the number of knots for each smooth (if any)}
#'   \item{npars - the length of the resulting parameter vector}
#'   \item{SmCons - a list of the smooth construct objects for each smooth (if any)}
#' }
#' The model matrix includes columns for both fixed effects (if any) and smooths (if any),
#' with fixed effect terms coming first. The number of columns reflects the number of actual
#' parameters associated with the model matrix. The model matrix has the same number of rows as the
#' model frame and `dfrMF`. The model frame is a dataframe consisting of a subset of the columns of `dfrMF`,
#' with the columns corresponding to any smooth variables converted to numeric type and possibly preceded
#' by a column of 1's with name '(Intercept)' if the formula includes an intercept.
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
  if (f!="~"){
    ###--create fixed factors model formula
    Frmla = Formula::Formula(eval(parse(text=f)));
    ###--create the model frame from the formula and the basis for the model frame
    MdFrm = Formula:::model.frame.Formula(Frmla,data=dfrMF);
    ###--convert character columns to factors and apply contrasts to model frame
    dfr   = convertToDimsFactors(MdFrm,NULL,contrasts=ctrs,verbose=verbose);
    #--create the model matrix
    MdMtx = Formula:::model.matrix.Formula(Frmla,data=dfr,rhs=NULL);
    ns = ncol(MdMtx); #--number of fixed effects parameters
  } #--otherwise no fixed factors

  ###--extract any smooths
  fs  = extractSmoothsFromFormulaText(txt);
  ks  = vector("integer",0);
  scs = list();
  if (length(fs)>0){
    for (i in 1:length(fs)){
      sso = eval(parse(text=fs[i])); #--smooth specification formula (requires package mgcv)
      vrs = extractSmoothVars(fs[i]);
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
  #--get model matrix reduced to unique combinations of elements
  dfrMM = MdFrm |> dplyr::bind_cols(MdMtx) |> dplyr::distinct();
  nps   = ncol(MdMtx); #--length of parameter vector
  nfs   = ncol(MdFrm); #--number of factors/variables (should be ns+ks)
  redMM = Matrix::Matrix(as.matrix(dfrMM[,(nfs+1):(nfs+nps)]));
  #--create output list
  res = list(ModMtx=redMM,dfrModMtx=dfrMM,fullMM=MdMtx,ModFrm=MdFrm,f=txt,
             ns=ns,ks=ks,npars=nps,SmCons=scs);
  return(res);
}

#TODO: fill in for random effects model matrix and associated info
#' @title Calculate a model matrix for random effects
#' @description Function to calculate a random effects model matrix (and associated information)
#' @param txt - model formula as text
#' @param dfrMF - basis dataframe for model frame
#' @param cov_type - covariance type
#' @param verbose - flag to print debugging info
#' @return list with model matrix and other information, or NULL if not `txt` is not a formula
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
  if ((is.na(txt))||(!(txt |> stringr::str_starts("~")))) return(npars=0);

  if (tolower(cov_type)=="iid"){

  } else
  if (tolower(cov_type)=="AR1"){

  } else {
    error(paste0("unrecognized `cov_type` in calcModelMatrixREs"))
  }

  ###--TODO: fill this in!
  res = list(naprs=0);#--for now!
  return(res);
}


