#--test reading new CTL file format
#--read in data file---
dirPrj =  rstudioapi::getActiveProject();
#dirPrj= getwd();
source(file.path(dirPrj,"R","MiscFunctions_ReadParse.R"));
source(file.path(dirPrj,"R","readADMB_DataFileFunctions.R"));
fn = file.path(dirPrj,"testing/inputFiles",
               "/TannerCrab/tcData_AllFisheriesNMFS-2SexCatchData.dat");
lstData = readADMB_DataFile(fn);

#--create model dimensions----
source(file.path(dirPrj,"R","MiscFunctions_Text.R"));
source(file.path(dirPrj,"R","DimensionsFunctions.R"));
source(file.path(dirPrj,"R","DimensionsUtilities.R"));
## The following string would be in the "dims"" file:
str=paste(
'MODEL_DIMS
  y <- 1948:2024;                           #--years
  s <- 1;                                   #--seasons
  r <- "EBS";                               #--regions
  x <- c("male","female");                  #--sex classes
  m <- c("immature","mature");              #--maturity state classes
  p <- c("new_shell","old_shell");          #--post-molt ages
  zc <- seq(24.5,184.5,5);                  #--size bin cutpoints
  z  <- 0.5*(zc[2:(length(zc))]+zc[1:(length(zc)-1)]); #--size bin midpoints
END'
)
strv = stringr::str_split_1(str,"\\n") |> extractLines("MODEL_DIMS","END");
lstDims = parseStrAsList(strv,verbose=TRUE);
dimsMdl = createSparseDimsMap(!!!(lstDims[names(lstDims)!="zc"]));#--drop zc (size bin cutpoints)
dimsZBC = tibble::tibble(lft=lstDims$zc[1:(length(lstDims$zc)-1)],#--size bin dims (left, center,right)
                         mid=lstDims$z,
                         rgt=lstDims$zc[2:(length(lstDims$zc))]);

#--define fleets
f = c("TCF","SCF","NMFS");               #--fleets

#--create year blocks----
source(file.path(dirPrj,"R","defineYearBlocks.R"));
str=paste(
'YEAR_BLOCKS
tbAll <- list(all=1949:2022) |>
         defineYearBlocks(dimsMdl=dimsMdl);
tbRec <- list(spinup=1949:1974,
              normal=1975:2022) |>
         defineYearBlocks(dimsMdl=dimsMdl);
tbNM  <- list(normal=c(1949:1979,1985:2022),
              heightened=1980:1984) |>
         defineYearBlocks(dimsMdl=dimsMdl);
tbTCFRetM<-list(early=c(1949:1983,1987:1989),
             mid  =c(1990:1996),
             PR   =c(2005:2009,2013:2015,2017:2018,2020:2022)) |>
           defineYearBlocks(dimsMdl=dimsMdl) |>
           dplyr::filter(x=="male");
tbTCFSelM<-list(early=c(1949:1983,1987:1989),
             mid  =c(1990:1996),
             PR   =c(2005:2009,2013:2015,2017:2018,2020:2022)) |>
           defineYearBlocks(dimsMdl=dimsMdl) |>
           dplyr::filter(x=="male") |> dplyr::mutate(f="TCF");
tbTCFSelF<-list(all=c(1949:2022)) |>
           defineYearBlocks(dimsMdl=dimsMdl) |>
           dplyr::filter(x=="female") |> dplyr::mutate(f="TCF");
tbNMFS<-list(preGC=c(1975:1981),
             postGC=c(1982:2019,2021:2022)) |>
        defineYearBlocks(dimsMdl=dimsMdl) |>
        dplyr::mutate(f="NMFS");
END')
strv = stringr::str_split_1(str,"\\n") |> extractLines("YEAR_BLOCKS","END");
lstYBlks = parseStrAsList(strv,verbose=FALSE);
for (nm in names(lstYBlks)){
  eval(parse(text=paste0(nm,"=lstYBlks[['",nm,"']]")));
}

#--create function/parameter definitions----
##--selectivity functions----
str=paste('
FUNCTIONS
fcn_idx  function     type          D    fcn_params  description                            #--What is "D"??
   1     selTCFF      ascnormal    TCFF    pZMd,pWd   size_at_mode,_width
   2     selTCM1      logistic     TCFM1   pZ50,pWd   size_at_50%_selected,_width
   4     selTCFM3     smooth       TCFM3     k=6      smooth_with_6_degrees_of_freedom
  -4     selTCFM4     smooth       TCFM4     k=6      mirrors_fcn_idx_4
   5     selSCM1      fxd_vector   SCFM1    values    vector_of_fixed_values
   6     selSCM2      fxd_matrix   SCFM2    values    matrix of fixed values
END'  #--terminator for functions block
);

### FUNCTION
### id - function id
### function - function name
### MF - model frame for fixed effects (including function-level covariates) and function-level random effects
### pars - fixed effects parameters (separater by commas, no spaces)
### REs - one-sided formula for function-level random effects (e.g. ~(1|yblk))
### dispREs - one-sided formula defining RE dispersion parameters
### linkREs - link function for random effects (add, mult, exp)
### description - function label (no spaces)
### PARAMETERS
### id - parameter id
### IV - initial value for leading parameter
### LB - lower bound on function scale
### UB - upper bound on function scale
### phz - parameter "phase"
### PriorType - name of prior function
### Pr1 - prior function mean parameter
### Pr2 - prior function dispersion parameter
### FEs - formula for fixed effects
### linkFEs - link function for fixed effects (add, exp)
### REs - one-sided formula for parameter-level random effects (e.g. ~(1|yblk))
### dispREs - one-sided formula defining RE dispersion parameters
### linkREs - link function for random effects (add, mult, exp)
### description - parameter label label (no spaces)

str=paste('
FUNCTION  #--must start section with "FUNCTION"
   # a comment line
# another  comment line
id           function      MF        pars     REs   dispREs  linkREs   description
selNMFS      ascnormal   tbNMFS    pZMd,pWd   ~0      ~1       add     ascending_normal
PARAMETERS
id    IV   LB     UB    phz   PriorType   Pr1    Pr2    FEs         contrasts        linkFEs       REs         dispREs  linkREs   description
pZMd  100   5    150     5    uniform      NA     NA    ~yblk          ident          ident       ~0            ~1     add     size_at_mode
pWd    20   5     50    -5    uniform      NA     NA    ~0+yblk+x   yblk=contr.sum    ident     ~(y|yblk)     ~yblk    add     width
INITIAL_VALUES
END')
strp = stringr::str_split_1(str,"\\n") |> removeCommentLines();
strv = strp |> extractLines("FUNCTION","PARAMETERS");
dfrFcn = readr::read_table(strv);
strv = strp |> extractLines("PARAMETERS","INITIAL_VALUES");
dfrPars = readr::read_table(strv);
strv = strp |> extractLines("INITIAL_VALUES","END") |> removeCommentLines();

for (rf in 1:nrow(dfrFcn)){
  #--testing: rf = 1;
  MF = dfrFcn[["MF"]][rf];
  dfrMF = eval(parse(text=MF));
  pars = stringr::str_split_1(dfrFcn[["pars"]][rf],",");
  for (rp in 1:length(pars)){
    #--testing: rp = 1;
    par = dfrPars[["id"]][rp];
    ctrs = dfrPars[["contrasts"]][rp];
    Frmla = Formula::Formula(formula(dfrPars[["FEs"]][rp]));
    MdFrm = Formula:::model.frame.Formula(Frmla,data=dfrMF);
    MdMtx = Formula:::model.matrix.Formula(Frmla,data=convertToDimsFactors(MdFrm,dimsMdl,contrasts=ctrs),rhs=NULL);
  }
}


