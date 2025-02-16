#'
#' @title Calculate the objective function for a gmacs model
#' @description Function to calculate the objective function for a gmacs model.
#' @param params - named list of parameter values (p=1 for now)
#'
#' @return The objective function (negative log-likelihood) value given the data and parameters
#'
#' @details Input data and options should be organized in the calling environment as a list of lists
#' with the name "inputs". `inputs` should have the structure
#' \itemize{
#'   \item{dims - a list with model dimensions maps}
#'   \item{data - a named list with data-related entries}
#'   \item{options - a named list with selected options}
#'   \item{testing - a logical to indicate the function is run in a testing environment}
#' }
#'
#' Example
#' ```{r, eval=FALSE}
#' inputs = list(dims=NULL,data=NULL,options=NULL,testing=TRUE);
#' params = list(dummy=0.0);
#' obj = MakeADFun(objfun,params,silent=TRUE);
#' opt = nlminb(obj$par, obj$fn, obj$gr);
#' sdreport(obj);
#' ````
#'
#' @import RTMB
#'
#' @export
#' @md
#'
#--Programming notes:----
##--Improving speed----
###--1. When building array one element at a time via loops, use x[[i,j,k]]<-y rather than x[i,j,k]<-y.
###--   use double bracket sub-assignment x[[i]]<-y rather than single bracket. Example:
###--  f <- function(x) {
###--      n <- 100
###--      a <- AD(array(0, c(n,n,n)))
###--      for (i in 1:n)
###--          for (j in 1:n)
###--              for (k in 1:n)
###--                  a[[i,j,k]] <- x
###--      sum(x)
###--  }
###--  system.time(F <- MakeTape(f, 0))
###--2. Combine inner loop vectorization and outer loop in-place insertion by using an array of lists.
###--   Then convert to normal array after construction:
###--  library(RTMB)
###--  n <- 100
###--  f <- function(x) {
###--      a <- array( list() , c(n, n) )
###--      for (i in 1:n)
###--          for (j in 1:n)
###--              a[[i,j]] <- x
###--      ## convert to normal array
###--      a <- do.call("c", a)
###--      dim(a) <- c(n, n, n)
###--      sum(a)
###--  }
###--  system.time(F <- MakeTape(f, numeric(n)))

##
obj_fun<-function(params){
  #--expand the list of parameters----
  #getAll(inputs,params); #--don't do this? more convenient not
  if (verbose) cat(objects(all.names=TRUE),"\n");

  #--inputs has elements
  #----dims, data, options, testing, verbose
  dims = intputs$dims;
  data = inputs$data;
  opts = inputs$options;

  #--get model indices-----
  nYs = dims$nYs;    #--number of years
  nSs = dims$nSs;    #--number of seasons
  nCs = dims$nCs;    #--number of population categories (size of population vector)
  nFs = dims$nFs     #--number of fleets
  I_c  = double(nCs);                #--identity vector for population state
  I_cc = diag(x=1,nrow=nCs,ncol=nCs);#--identity matrix for population state

  #--get options
  trackCohorts = opts$trackCohorts;#--T/F to track cohorts
  maxCohortAge = opts$maxCohortAge;#--max age for a cohort

  #--set up model processes----
  ##--initial N's----
  nInit = getInitN(dims,opts,params,verbose);
  REPORT(nInit);

  if (FALSE){  ##--remove FALSE when ready
    ##--Allometry----
    lstAlm = getAllometry(dims,opts,params,verbose);

    ##--Recruitment----
    lstRec = getRec(dims,options,params,verbose)

    ##--Natural mortality----
    lstNMRs = getNatMortRates(dims,options,params,verbose)

    ##--Growth----
    ###--will calculate all growth transition matrices
    ###--includes molt probabilities, growth matrices,
    ###--maturity transition matrices, probability of terminal molt
    lstGTMs = getGrowthTMs(dims,options,params,verbose)

    ##--Selectivity/Retention----
    ###--will calculate all selectivity curves
    lstSelFcns = getSelFcns(dims,options$sel_fcns,sel_params,verbose);

    ##--Fishery characteristics----
    ###--will calculate capture, handling mortality, and retention rates
    ###--for all fisheries
    lstFRs = getFisheryRates(dims,options,params,lstSelFcns,verbose);

    ##--Survey characteristics----
    ###--will calculate selectivity, q for all surveys
    lstSQs = getSurveyCatchability(dims,options,params,lstSelFcns,verbose);

    ##--Movement transition matrices----
    ###--will calculate all movement transition matrices
    ###--(to be implemented at some point in future)
    lstMTMs = getMovementTMs(dims,options,params,verbose);

    #--set up population (and cohort, requested) structures----
    lstCohs = list();
    if (trackCohorts){
      ##--set up cohorts----
      for (y_ in names(dims$y)){ #--`y_` is model year represented as a string (i.e., "2020", not 2020)
        nYCs = min(maxCohortAge+1,as.numeric(dplyr::last(dims$y))-as.numeric(y_)+1);
        lstCohs[[y_]] = list(y     = as.numeric(y_)+(0:(nYCs-1)),#--years cohort exists (numerical)
                             n_ysc = AD(array(0,nYCs,nSs,nCs)),  #--cohort abundance at beginning of season `s` in cohort year `y`
                             b_ysc = AD(array(0,nYCs,nSs,nCs)),  #--cohort biomass   at beginning of season `s` in cohort year `y`
                             n_Final = AD(array(0,nCs)),         #--cohort abundance at start of `dims$y+1`
                             b_Final = AD(array(0,nCs)),         #--cohort biomass   at start of `dims$y+1`
                             srvN_ysc = AD(array(0,nFs,nYCs,nSs,nCs)),#--cohort survey abundance in season `s` in cohort year `y`
                             srvB_ysc = AD(array(0,nFs,nYCs,nSs,nCs)),#--cohort survey biomass   in season `s` in cohort year `y`
                             fcn_fysc = AD(array(0,nFs,nYCs,nSs,nCs)),#--cohort fishery catch abundance           in season `s` in cohort year `y`
                             frn_fysc = AD(array(0,nFs,nYCs,nSs,nCs)),#--cohort fishery retained abundance        in season `s` in cohort year `y`
                             fdn_fysc = AD(array(0,nFs,nYCs,nSs,nCs)),#--cohort fishery discard abundance         in season `s` in cohort year `y`
                             fmn_fysc = AD(array(0,nFs,nYCs,nSs,nCs)),#--cohort fishery catch mortality abundance in season `s` in cohort year `y`
                             ssb_ysc  = AD(array(0,nYCs,nSs,nCs)),    #--cohort spawning stock biomass at of season `s` in cohort year `y`
                            );
      }#--y_ loop
    } else {
      ##--no cohort structure----
      lstCohs[["all"]]=list(y = dims$y,                       #--population years (numerical, with names "yyyy" for yyyy)
                            n_ysc = AD(array(0,nYs,nSs,nCs)), #--population abundance at beginning season `s` in year `y`
                            b_ysc = AD(array(0,nYs,nSs,nCs)), #--population biomass   at beginning of season `s` in model year `y`
                            n_Final = AD(array(0,nCs)),         #--population abundance at start of `last(dims$y)+1`
                            b_Final = AD(array(0,nCs)),         #--population biomass   at start of `last(dims$y)+1`
                            srvN_ysc = AD(array(0,nFs,nYs,nSs,nCs)),#--population survey abundance in season `s` in model year `y`
                            srvB_ysc = AD(array(0,nFs,nYs,nSs,nCs)),#--population survey biomass   in season `s` in model year `y`
                            fcn_fysc = AD(array(0,nFs,nYs,nSs,nCs)),#--population fishery catch abundance           in season `s` in model year `y`
                            frn_fysc = AD(array(0,nFs,nYs,nSs,nCs)),#--population fishery retained abundance        in season `s` in model year `y`
                            fdn_fysc = AD(array(0,nFs,nYs,nSs,nCs)),#--population fishery discard abundance         in season `s` in model year `y`
                            fmn_fysc = AD(array(0,nFs,nYs,nSs,nCs)),#--population fishery catch mortality abundance in season `s` in model year `y`
                            ssb_ysc  = AD(array(0,nYs,nSs,nCs)),    #--population spawning stock biomass at of season `s` in model year `y`
                            );
    }#--cohort/no cohort structure

    #--integrate model over years by season----
    ##--assign initial abundance----
    if (trackCohorts){
      #--TBD: assign abundance to initial cohorts, if not building up pop by recruitment only
      ##--nInit could be generated by equilibrium considerations for first `maxCohortAge` cohorts
      ##--or (tediously) assigned values
    } else {
      lstCohs[["all"]]$n_ysc[1,1,] = nInit;
    }

    for (y_ in dims$y){ ##--loop over years----
      for (s_ in dims$s){ ###--loop over seasons within years----
        n_c = n_ysc[y_,s_];#--population state at start of y_, s_
        if (verbose) cat("y, s =",y_,s_,"\n");
        for (cht in names(lstCohs)){ ####--loop over cohorts----
          if (y_ %in% lstCohs[[cht]]$y){
            n_cp  = lstCohs[[cht]]$n_ysc;#--cohort state at start of y_, s_

            if (seasonLength[y_,s_]==0){  #####--instantaneous processes----
              ######---happen in order of
              ######----survey, growth, movement, and/or (possibly) fishery catches and mortality

              ######--surveys----
              for (v_ in 1:nSrvs){
                if (doSurveys_vys[v_,y_,s_]){
                  lstCohs[[cht]]$srvN_vysc[v_,y_,s_,] = lstSQs[[v_]]$q_ysc[y_,s_,] * n_cp;
                  lstCohs[[cht]]$srvB_vysc[v_,y_,s_,] = lstAlm$wAtZ_ysc[y_,s_,] * srvN_vysc[v_,y_,s_,];
                  srvN_vysc[v_,y_,s_,] = srvN_vysc[v_,y_,s_,] + lstCohs[[cht]]$srvN_vysc[v_,y_,s_,]
                  srvB_vysc[v_,y_,s_,] = srvB_vysc[v_,y_,s_,] + lstCohs[[cht]]$srvB_vysc[v_,y_,s_,]
                }
              }#--v_

              ######--growth----
              ######--might want to use [y_,s_] matrices with pointers to
              ######--actual growth process vectors/matrices to reduce storage reqts
              if (doGrowth_ys[y_,s_]){
                ######--determine which crab will/will not molt----
                nn_cp = 1-lstGTMs$prMolt[y,s,] * n_cp; #--crab that will skip the molt
                mn_cp =   lstGTMs$prMolt[y,s,] * n_cp; #--crab that will molt
                ######--increment post-molt age for non-molting crab----
                nn_cp = inputs$growth$tmPMAI %*% nn_cp;
                ######--determine which crab will undergo the molt-to-maturity based on pre-molt size
                if (tolower(opts$maturity$molt2maturity)=="pre-molt_size")
                  mn_cp = lstGTMs$prM2M[y,s,,] %*% mn_cp;
                ######--do molt: apply growth matrix to crab that will molt----
                mn_cp = lstGTMs$prGr[y,s,,] %*% mn_cp;
                ######--decrement post-molt age to 0 for molting crab----
                mn_cp = inputs$growth$tmPMAD %*% mn_cp;
                ######--determine which crab will undergo the molt-to-maturity based on pre-molt size
                if (tolower(opts$maturity$molt2maturity)=="post-molt_size")
                  mn_cp = lstGTMs$prM2M[y,s,,] %*% mn_cp;
              }

              ######--instantaneous movement----
              if (doInstMovement[y_,s_]){
                n_cp = lstMTMs$mtm_cc[y_,s_,,] %*% n_cp;
              }

              ######--instantaneous fishery catches and mortality----
              #--calculate total fishing capture, retention, discard mortality, and total mortality rates
              #######--calculate total fishing mortality rates (FMR) by population component----
              totFMR_c = AD(vector("numeric",nCs));#--total FMRs
              for (f_ in 1:nFsh){
                if (doInstFisheries_fys[f_,y_,s_]){
                  totFMR_c = totFMR_c + lstFMs[[f_]]$fmr_ysc[y_,s_,];#--add fishery-specific FMRs to total FMRs
                }
              }#--f_
              #######--calculate fishery-specific catch-related numbers----
              scl = (1.0-exp(-totFMR_c))/(totFMR_c+(1*(totFMR_c==0)));#--scaling factor (last bit excludes divides by 0)
              for (f_ in 1:nFsh){
                if (doInstFisheries_fys[f_,y_,s_]){
                  lstCohs[[cht]]$fcn_fysc[f_,y_,s_,] =  lstFMs[[f_]]$fcr_ysc[y_,s_,]*scl*n_cp;#--captured numbers
                  lstCohs[[cht]]$frn_fysc[f_,y_,s_,] =  lstFMs[[f_]]$frr_ysc[y_,s_,]*scl*n_cp;#--retained numbers
                  lstCohs[[cht]]$fdn_fysc[f_,y_,s_,] =  lstFMs[[f_]]$fdr_ysc[y_,s_,]*scl*n_cp;#--discard mortality numbers
                  lstCohs[[cht]]$fmn_fysc[f_,y_,s_,] =  lstFMs[[f_]]$fmr_ysc[y_,s_,]*scl*n_cp;#--retained + discard mortality numbers
                  fcn_fysc[f_,y_,s_,] = fcn_fysc[f_,y_,s_,] + lstCohs[[cht]]$fcn_fysc[f_,y_,s_,];
                  frn_fysc[f_,y_,s_,] = frn_fysc[f_,y_,s_,] + lstCohs[[cht]]$frn_fysc[f_,y_,s_,];
                  fdn_fysc[f_,y_,s_,] = fdn_fysc[f_,y_,s_,] + lstCohs[[cht]]$fdn_fysc[f_,y_,s_,];
                  fmn_fysc[f_,y_,s_,] = fmn_fysc[f_,y_,s_,] + lstCohs[[cht]]$fmn_fysc[f_,y_,s_,];
                }
              }#--instantaneous fisheries
              n_cp = n_cp * exp(-totFMR_c);#--pop numbers after fisheries

            } else { #####--continuous processes happen throughout season----

              ######--continuous process can include
              ######--natural mortality and continuous fishery catch and mortality
              #######--NOTE: Incorporating movement here requires exponentiating a matrix
              #######--to calculate changes to n_c. Not quite sure how to calculate fishery-specific
              #######--values (catch, mortality) in this case: is it directly analagous to the
              #######--case without movement??
              dt = seasonLength[y_,s_];

              ######--continuous movement rates----
              # if (doContMovement[y_,s_])
              #   mtm_cc = lstMTMs$mtm_yscc[y_,s_,,];#--c *x* c matrix

              ######--natural mortality rate----
              nmr_c = lstNMRs$nmr_ysc[y_,s_,];#---vector

              ######--continuous fishery catch rates
              #######--calculate total fishing mortality rates (FMR) by population component----
              totFMR_c = AD(vector("numeric",nCs));#--total FMRs (vector)
              for (f_ in 1:nFsh){
                if (doContFisheries_fys[f_,y_,s_]){
                  totFMR_c = totFMR_c + lstFMs[[f_]]$fmr_ysc[y_,s_,];#--add fishery-specific FMRs to total FMRs
                }
              }#--f_ loop
              #######--calculate fishery-specific catch-related numbers----
              #######--these may not be correct if continuous movement is included
              scl = (1.0-exp(-(totFMR_c+nmr_c)*dt))/(totFMR_c+nmr_c);#--scaling factor
              for (f_ in 1:nFsh){
                if (doContFisheries_fys[f_,y_,s_]){
                  lstCohs[[cht]]$fcn_fysc[f_,y_,s_,] =  lstFMs[[f_]]$fcr_ysc[y_,s_,]*scl*n_cp;#--captured numbers
                  lstCohs[[cht]]$frn_fysc[f_,y_,s_,] =  lstFMs[[f_]]$frr_ysc[y_,s_,]*scl*n_cp;#--retained numbers
                  lstCohs[[cht]]$fdn_fysc[f_,y_,s_,] =  lstFMs[[f_]]$fdr_ysc[y_,s_,]*scl*n_cp;#--discard mortality numbers
                  lstCohs[[cht]]$fmn_fysc[f_,y_,s_,] =  lstFMs[[f_]]$fmr_ysc[y_,s_,]*scl*n_cp;#--retained + discard mortality numbers
                  fcn_fysc[f_,y_,s_,] = fcn_fysc[f_,y_,s_,] + lstCohs[[cht]]$fcn_fysc[f_,y_,s_,];
                  frn_fysc[f_,y_,s_,] = frn_fysc[f_,y_,s_,] + lstCohs[[cht]]$frn_fysc[f_,y_,s_,];
                  fdn_fysc[f_,y_,s_,] = fdn_fysc[f_,y_,s_,] + lstCohs[[cht]]$fdn_fysc[f_,y_,s_,];
                  fmn_fysc[f_,y_,s_,] = fmn_fysc[f_,y_,s_,] + lstCohs[[cht]]$fmn_fysc[f_,y_,s_,];
                }
              }#--continuous fisheries

              ######--update N for continuous processes----
              #######--if continuous movement is implemented,
              #######--the following needs to be converted to an exponentiated matrix (rather thn vector)
              n_cp = exp(-(nmr_c+totFMR_c)*dt) * n_cp;
            }#--if (seasonLength[y_,s_])

            #####--Finally, add recruitment----
            if (doRecruitment[y_,s_]) {
              if (trackCohorts && (y_==lstCohs[[cht]]$y[1])){
                #--only cohort in first year has recruitment
                lstCohs[[cht]]$n_ysc[y_,s_,] = lstRec$rec_ysc[y_,s_,];
                n_cp = n_cp + lstRec$rec_ysc[y_,s_,];
              } else {
                #--"cohort" is really the population, so don't worry about "cohort" age
                lstCohs[[cht]]$n_ysc[y_,s_,] = lstCohs[[cht]]$n_ysc[y_,s_,] + lstRec$rec_ysc[y_,s_,];
                n_cp = n_cp + lstRec$rec_ysc[y_,s_,];
              }
            }#--doRecruitment

            #####--finished all processes in current season, n_cp is cohort (or pop) abundance at start of next season
            if (s<nSs) {#####--not at end of year----
              lstCohs[[cht]]$n_ysc[y_,s_+1,] = n_cp;                    #--cohort (or pop) state at start of y_, s_+1
              lstCohs[[cht]]$b_ysc[y_,s_+1,] = lstAlm$wAtZ[y_,s_,]*n_cp;#--biomass at start of y_, s_+1
            } else {    #####--at end of current year, so need to increment year if cohort age or model timeframe allows----
              if (trackCohorts){
                if (y_ < last(lstCohs[[cht]]$y)){ ######--have not reached last age
                  lstCohs[[cht]]$n_ysc[y_+1,s_,] = n_cp;                    #--cohort state at start of y_+1
                  lstCohs[[cht]]$b_ysc[y_+1,s_,] = lstAlm$wAtZ[y_,s_,]*n_cp;#--cohort biomass at start of y_+1
                } else if (y_ < last(dims$y)){ ######--reached last age, so cohort does not go forward
                  n_cp = 0*n_cp;
                } else { ######--in final year, save final cohort state
                  lstCohs[[cht]]$nFinal = n_cp;                    #--cohort state at start of y_+1
                  lstCohs[[cht]]$bFinal = lstAlm$wAtZ[y_,s_,]*n_cp;#--cohort biomass at start of y_+1
                }
              } else { #####--"cohort" is population
                if (y_ < last(dims$y)) { ######--save pop state and biomass at start of y_ + 1
                  lstCohs[[cht]]$n_ysc[y_+1,1,] = n_cp;                    #--pop state at start of y_+1
                  lstCohs[[cht]]$b_ysc[y_+1,1,] = lstAlm$wAtZ[y_,s_,]*n_cp;#--pop biomass at start of y_+1
                } else { ######--in final year, save final "cohort" state
                  lstCohs[[cht]]$nFinal = n_cp;                    #--cohort state at start of y_+1
                  lstCohs[[cht]]$bFinal = lstAlm$wAtZ[y_,s_,]*n_cp;#--cohort biomass at start of y_+1
                }
              }#--if trackCohorts or not
            }#--if s < nSs or not

            #####--add cohort changes to population changes----
            if (trackCohorts) n_c = n_c + n_cp;
          }#--y_ in lstCoh[[cht]]$y
        }#--end cohort loop

        if (s_<nSs) { ##--if (!trackChanges) these should be same as equivalents in lstCohs[["all]]
          n_ysc[y_,s_+1,] = n_c;                    #--save pop state (integrated over cohorts) at start of y_, s_+1
          b_ysc[y_,s_+1,] = lstAlm$wAtZ[y_,s_,]*n_c;#--biomass at start of y_, s_+1
        }
      } ####--end seasons loop (s_)----

      if (y_<nYs) { ##--if (!trackChanges) these should be same as equivalents in lstCohs[["all]]
        n_ysc[y_+1,1,] = n_c;                    #--save pop state as start of y_+1, s_
        b_ysc[y_+1,1,] = lstAlm$wAtZ[y_,s_,]*n_c;#--biomass at start of y_, s_+1
      }
    } ###--end years loop (y_)----

    ##--need to save final pop state (start of year nYs+1, season 1; integrated over cohorts)----
    ##--if (!trackChanges) these should be same as equivalents in lstCohs[["all]]
    nFinal = n_c;                      #--final pop abundance
    bFinal = lstAlm$wAtZ[nYs,nSs,]*n_c;#--final pop biomass
  }#--FALSE block

  #--calculate likelihoods----

  if (testing){
    nll = -dnorm(1,dummy,1,log=TRUE);
  }
  if (verbose) cat("end objective function\n");
  return(nll);
}#--end of objective_function

if (FALSE){
  #--for simple testing----
  require(rtmbGMACS);
  dirPrj = rstudioapi::getActiveProject();
  #--start opts list----
  opts = list(tracKCohorts=FALSE,
              maxCohortAge=10);
  #--set up dims----
  source(file.path(dirPrj,"testing","r_setupModelDimensions.TestA.R"));
  dims = setupModelDims();

  #--set up allometry
  source(file.path(dirPrj,"R","readParamInfo_Allometry.R"));
  source(file.path(dirPrj,"R","extractParamInfo_Allometry.R"));
  conn     = file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.function.txt");
  res      = readParamInfo_Allometry(conn,TRUE);
  lstAllom = extractParamInfo_Allometry(res,dims$dmsYSC);
}
