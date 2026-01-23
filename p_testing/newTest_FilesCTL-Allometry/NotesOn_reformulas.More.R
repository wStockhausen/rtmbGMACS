#--test defining formulas as text

reformulas:::.valid_covstruct #--named number vector
#   diag      us      cs     ar1      ou     exp     gau     mat    toep      rr homdiag
#      0       1       2       3       4       5       6       7       8       9      10
glmmTMB:::.valid_covstruct   #--adds a couple to extra to above
#   diag      us      cs     ar1      ou     exp     gau     mat    toep      rr homdiag  propto  hetar1   homcs homtoep
#      0       1       2       3       4       5       6       7       8       9      10      11      12      13      14

reformulas:::findReTrmClasses(); #--character vector (c(names(.valid_covstruct),"s")
# [1] "diag"    "us"      "cs"      "ar1"     "ou"      "exp"     "gau"     "mat"     "toep"    "rr"      "homdiag" "s"
##--Probably want to add other "specials", e.g. "ti"
glmmTMB:::.

#--`lme4`-style formulas with `reformulas` functions
##--fixed effects only formula
f1s = "~1+x*m";
f1 = as.formula(f1s);
str(f1);
f1.sf = reformulas::splitForm(f1)
str(f1.sf);

##--mixed effects formula
f2s = "~1+x*m+(y|z)";
f2 = as.formula(f2s);
f2;
str(f2);
f2.sf = reformulas::splitForm(f2)
str(f2.sf);
f2.sf$fixedFormula;
f2.sf$reTrmFormulas;
f2.sf$reTrmAddArgs;
f2.sf$reTrmClasses;

##--mixed effects formula specifying covstruct for RE term
f2s = "~1+x*m+diag(y|z)";
f2 = as.formula(f2s);
f2;
str(f2);
f2.sf = reformulas::splitForm(f2)
str(f2.sf);
f2.sf$fixedFormula;
f2.sf$reTrmFormulas;
f2.sf$reTrmAddArgs;
f2.sf$reTrmClasses;

##--mixed effects formula specifying covstruct for RE terms and double bars for last RE term
###--single bar  implies RE terms left of "|" are correlated   (covariance matrix is not diagonal)
###--double bars implies RE terms left of "||" are independent (covariance matrix is (block?) diagonal)
f2s = "~1+x*m+cs(y|z) + (1+a||b)";
f2 = as.formula(f2s);
f2;
str(f2);
f2.sf = reformulas::splitForm(f2)
str(f2.sf);
f2.sf$fixedFormula;
f2.sf$reTrmFormulas;
f2.sf$reTrmAddArgs;
f2.sf$reTrmClasses;

##--mixed effects formula specifying a "special' ("s"), a covstruct for the first RE term, and double bars for last RE term
###--single bar  implies RE terms left of "|" are correlated   (covariance matrix is not diagonal)
###--double bars implies RE terms left of "||" are independent (covariance matrix is (block?) diagonal)
f2s = "~1 + x*m + s(q)+ cs(y|z) + (1+a|b) + (c||d) + cs(1+e||f)";
f2 = as.formula(f2s);
f2;
str(f2);
###--splitForm
####--returns list with elements fixedFormula, reTrmFormulas, reTrmArgs, reTrmClasses
####--fixedFormula: `formula` object (with default `specials`, does not include s(...) terms)
####--reTrmFormulas: list with RE terms as language objects (with default `specials`, includes terms inside s(...) terms)
####--reTrmAddArgs:
####--reTrmClasses
####--Presumably s(...) should be included in fixedFormula
f2.sf = reformulas::splitForm(f2);
str(f2.sf);
f2.sf$fixedFormula;  #--does not include s(...) terms ("s" is identified as a "special" in reformulas::findReTrmClasses)
f2.sf$reTrmFormulas; #--list of "ordinary" RE terms and terms within "specials" (e.g., s(...))
f2.sf$reTrmAddArgs;  #--list of covstruct functions
f2.sf$reTrmClasses;  #--vector of covstruct class names
###--`findbars`
f2.fb = reformulas::findbars(f2);#--returns RE terms as list of language objects with no covstruct info
str(f2.fb)

#--formula processing for rtmbGMACS
##--set up model dimensions----
dirPrj = rstudioapi::getActiveProject();
source(file.path(dirPrj,"testing/newTest_FilesCTL-Allometry","r01_SetupDims.R"));
head(dmsYSC);
names(dmsYSC);
###--create base model frame dmsYSC by adding a couple of columns
mf = dmsYSC |> dplyr::mutate(yn=as.numeric(y),
                             zn=as.numeric(z),
                            yblk=ifelse(yn<2019,"<2019",
                                        ifelse(yn<2022,"<2022",">=2022")));
##--example FE formulas
###--numerical values for cos(zn) are created as part of the `call` creating the model frame
mfp = createModelFrame(as.formula("~1 + x:m + cos(zn)"),mf);
##--example RE formulas
###--random slope: numeric covariate `yn` within `x` grouping
mfp = createModelFrame(as.formula("~(yn|x)"),mf);
###--REs for (yn|yblk) specified as having an AR1 covariance structure
mfp = createModelFrame(as.formula("~ar1(1|x)"),mf);
##--mixed effects parameter formula with yn, cos(zn) as numeric covariates, yblk and x as RE grouping factors
f = as.formula("~1 + x:m + s(zn) + ar1(yn|yblk) + cos(zn) + (yn||x)"); #--fake formula
###--substitute "safe" characters ('+') for "specials" in order to create a model frame
fp = reformulas::sub_specials(f,
                              specials=c("|","||","s","ar1"),              #--names of specials to process ("ar1" not in default: add to a "specials" enum?)
                              keep_args=c(2L,2L,NA_integer_,NA_integer_)); #--number of arguments to retain (by special)
#  fp = ~1 + x:m + +zn + +(yn + yblk) + cos(zn) + (yn + x) #--note s and ar1 have been removed but cos(zn) has been retained
source(file.path(dirPj,"R/MiscFunctions_parameters.R"));
mfp = createModelFrame(fp,mf);
mfp = createModelFrame(f,mf);

##--create model matrices using reformulas::splitForm and mkRetrms----
###--split fe and RE terms----
f  = as.formula("~1 + x:m + s(zn) + ar1(1+yn|yblk) + cos(zn) + (yn||x)"); #--fake formula
fp = reformulas::splitForm(f);
f_fe  = fp$fixedFormula;
f_res = list();
for (i in seq_along(fp$reTrmFormulas)){
  f_res[[i]] = list(f=fp$reTrmFormulas[[i]],addargs=fp$reTrmAddArgs[[i]],class=fp$reTrmClasses[i])
}
###--make FE model matrix----
mf_fe = createModelFrame(f_fe,mf);
mm_fe = model.matrix.default(f_fe,mf_fe,contrasts.arg=NULL);
colnames(mm_fe);
###--RE terms----
bars = reformulas::findbars(f);
mf_re = createModelFrame(bars[1],mf);
mf_re = createModelFrame(f_res[[1]]$f,mf);

###--using mkReTrms----
mf_re <- model.frame(reformulas::sub_specials(f,specials=c("|","||","s",names(reformulas:::.valid_covstruct))), data = mf);
bars = reformulas::findbars(f);
reTrm = mkReTrms(bars[1],fr=mf_re)
head(t(reTrm$Zt)); #--head of model matrix (Z)
head(t(reTrm$Ztlist$`1 + yn | yblk`)) #--head of model matrix for individual RE terms
reTrm$Lambdat; #--transpose of sparse relative covariance factor
reTrm$theta;#--initial values of covariance parameters
reTrm$lower;#--lower bounds on covariance parameters
for (i in seq_along(reTrm$flist)) {
  cat(names(reTrm$flist)[i],"\n",levels(reTrm$flist[[i]]),"\n");#--ith element of list of grouping factors
  cat(reTrm$cnms[[i]],"\n");                                    #--ith element of list of column names for RE by grouping factor
}
reTrm$Gp;
reTrm$nl;
reTrm$ord;

##--using getXReTrms (from glmmmTMB)----
###--set up model dimensions----
dirPrj = rstudioapi::getActiveProject();
source(file.path(dirPrj,"testing/newTest_FilesCTL-Allometry","r01_SetupDims.R"));
head(dmsYSC);
names(dmsYSC);
####--create base model frame dmsYSC by adding a couple of columns
fr = dmsYSC |> dplyr::mutate(yn=as.numeric(y),
                             zn=as.numeric(z),
                            yblk=ifelse(yn<2019,"<2019",
                                        ifelse(yn<2022,"<2022",">=2022")));
###--set up inputs to my and glmmTMB versions of getXReTerms
mf = NULL;             #--call necessary for glmmTMB version
ranOK=TRUE;            #--flag indicating random effects allowed here
type="";               #--label for model type
contrasts=NULL;        #--a list of contrasts (see ?glmmTMB)
sparse=FALSE;          #--flag to return sparse model matrix?
old_smooths = NULL;    #--smooth information from a prior model fit (for prediction)
require(mgcv); require(Matrix);
###--create formulas and inputs for getXReTrms----
####--glmmTMB: formula without smooth----
formula  = as.formula("resp~1 + p + x:m + ar1(0+y|yblk) + cos(zn) + (1+yn||x)");
data = dplyr::bind_cols(fr,resp=0.0);
mf       = match.call(glmmTMB,call("glmmTMB",formula,data=data));
res_nsm = glmmTMB:::getXReTrms(formula, mf, data, ranOK=ranOK, type=type,
                               contrasts=contrasts, sparse=sparse, old_smooths = old_smooths);
View(res_nsm);
View(res_nsm$X);
View(as.matrix(res_nsm$Z));
View(as.matrix(res_nsm$reTrms$Lambdat)); #--relative covariance factor matrix
res_nsm$reTrms$Lambdat@x;                #--relative covariance factor matrix values as vector
res_nsm$reTrms$theta;                    #--initial covariance parameters values
res_nsm$reTrms$Lind;                     #--indices mapping theta to Lambdat@x
res_nsm$reTrms$theta[res_nsm$reTrms$Lind];   #--mapped values
res_nsm$reTrms$Lambdat@x<-res_nsm$reTrms$theta[res_nsm$reTrms$Lind];#--assign initial theta to relative covariance factor
View(as.matrix(res_nsm$reTrms$Lambdat));#--initial relative covariance factor
res_nsm$reTrms$theta<-as.double(seq_along(res_nsm$reTrms$theta)); #--revise initial values to better see mapping to RCF
res_nsm$reTrms$Lambdat@x<-res_nsm$reTrms$theta[res_nsm$reTrms$Lind];#--assign initial theta to relative covariance factor
View(as.matrix(res_nsm$reTrms$Lambdat));#--initial relative covariance factor

####--my version: formula without smooth----
formula  = as.formula("~1 + p + x:m + ar1(0+y|yblk) + cos(zn) + cos(zn):x + (1+yn||x)");
res_nsm = getXReTrms(formula, fr, ranOK=ranOK, type=type,
                       contrasts=contrasts, sparse=sparse, old_smooths = old_smooths);
View(res_nsm);
View(res_nsm$X);
View(as.matrix(res_nsm$Z));
View(as.matrix(res_nsm$reTrms$Lambdat)); #--relative covariance factor matrix
res_nsm$reTrms$Lambdat@x;                #--relative covariance factor matrix values as vector
res_nsm$reTrms$theta;                    #--initial covariance parameters values
res_nsm$reTrms$Lind;                     #--indices mapping theta to Lambdat@x
res_nsm$reTrms$theta[res_nsm$reTrms$Lind];   #--mapped values
res_nsm$reTrms$Lambdat@x<-res_nsm$reTrms$theta[res_nsm$reTrms$Lind];#--assign initial theta to relative covariance factor
View(as.matrix(res_nsm$reTrms$Lambdat));#--initial relative covariance factor
res_nsm$reTrms$theta<-as.double(seq_along(res_nsm$reTrms$theta)); #--revise initial values to better see mapping to RCF
res_nsm$reTrms$Lambdat@x<-res_nsm$reTrms$theta[res_nsm$reTrms$Lind];#--assign initial theta to relative covariance factor
View(as.matrix(res_nsm$reTrms$Lambdat));#--initial relative covariance factor

####--my version: formula with REs only----
formula  = as.formula("~ar1(0+y|yblk) + (1+yn||x)");
res_res = getXReTrms(formula, fr, ranOK=ranOK, type=type,
                       contrasts=contrasts, sparse=sparse, old_smooths = old_smooths);
View(res_res);
View(res_res$X);
View(as.matrix(res_res$Z));
View(as.matrix(res_res$reTrms$Lambdat)); #--relative covariance factor matrix
res_res$reTrms$Lambdat@x;                #--relative covariance factor matrix values as vector
res_res$reTrms$theta;                    #--initial covariance parameters values
res_res$reTrms$Lind;                     #--indices mapping theta to Lambdat@x
res_res$reTrms$theta[res_res$reTrms$Lind];   #--mapped values
res_res$reTrms$Lambdat@x<-res_res$reTrms$theta[res_res$reTrms$Lind];#--assign initial theta to relative covariance factor
View(as.matrix(res_res$reTrms$Lambdat));#--initial relative covariance factor
res_res$reTrms$theta<-as.double(seq_along(res_res$reTrms$theta)); #--revise initial values to better see mapping to RCF
res_res$reTrms$Lambdat@x<-res_res$reTrms$theta[res_res$reTrms$Lind];#--assign initial theta to relative covariance factor
View(as.matrix(res_res$reTrms$Lambdat));#--initial relative covariance factor

####--formula with smooth only----
formula  = as.formula("~0+ s(zn)");
res_wsm = getXReTrms(formula, fr, ranOK=ranOK, type=type,
                       contrasts=contrasts, sparse=sparse, old_smooths = old_smooths);
View(res_wsm);
View(res_wsm$X);
View(as.matrix(res_wsm$Z));
#####--TODO: need to add theta, Lambdat, Lind to getXReTerms for smooths
# View(as.matrix(res_wsm$reTrms$Lambdat)); #--relative covariance factor matrix
# res_wsm$reTrms$Lambdat@x;                #--relative covariance factor matrix values as vector
# res_wsm$reTrms$theta;                    #--initial covariance parameters values
# res_wsm$reTrms$Lind;                     #--indices mapping theta to Lambdat@x
# res_wsm$reTrms$theta[res_wsm$reTrms$Lind];   #--mapped values
# res_wsm$reTrms$Lambdat@x<-res_wsm$reTrms$theta[res_wsm$reTrms$Lind];#--assign initial theta to relative covariance factor
# View(as.matrix(res_wsm$reTrms$Lambdat));#--initial relative covariance factor
# res_wsm$reTrms$theta<-as.double(seq_along(res_wsm$reTrms$theta)); #--revise initial values to better see mapping to RCF
# res_wsm$reTrms$Lambdat@x<-res_wsm$reTrms$theta[res_wsm$reTrms$Lind];#--assign initial theta to relative covariance factor
# View(as.matrix(res_wsm$reTrms$Lambdat));#--initial relative covariance factor

####--formula w/ REs and smooths----
formula  = as.formula("~1 + x:m + s(zn) + ar1(y|yblk) + cos(zn) + (1+yn||x)"); #--fake formula
###--get FE and RE terms
res_all = getXReTrms(formula, fr, ranOK=ranOK, type=type,
                       contrasts=contrasts, sparse=sparse, old_smooths = old_smooths);
View(res_all);
View(res_all$X);
View(as.matrix(res_all$Z));
#####--TODO: need to add theta, Lambdat, Lind to getXReTerms for smooths
View(as.matrix(res_all$reTrms$Lambdat)); #--relative covariance factor matrix
res_all$reTrms$Lambdat@x;                #--relative covariance factor matrix values as vector
res_all$reTrms$theta;                    #--initial covariance parameters values
res_all$reTrms$Lind;                     #--indices mapping theta to Lambdat@x
res_all$reTrms$theta[res_all$reTrms$Lind];   #--mapped values
res_all$reTrms$Lambdat@x<-res_all$reTrms$theta[res_all$reTrms$Lind];#--assign initial theta to relative covariance factor
View(as.matrix(res_all$reTrms$Lambdat));#--initial relative covariance factor
res_all$reTrms$theta<-as.double(seq_along(res_all$reTrms$theta)); #--revise initial values to better see mapping to RCF
res_all$reTrms$Lambdat@x<-res_all$reTrms$theta[res_all$reTrms$Lind];#--assign initial theta to relative covariance factor
View(as.matrix(res_all$reTrms$Lambdat));#--initial relative covariance factor



