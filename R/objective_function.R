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
###--When building array one element at a time via loops, use x[[i,j,k]]<-y rather than x[i,j,k]<-y
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

  #--set up model indices-----
  nYs = length(dims$y);    #--number of years
  nSs = length(dims$s);    #--number of seasons
  nFs = length(dims$f)     #--number of fleets
  nRs = length(dims$r);    #--number of regions
  nXs = length(dims$x);    #--number of sex classes
  nMs = length(dims$m);    #--number of maturity state classes
  nAs = length(dims$a);    #--number of post-recruitment age classes
  nPs = length(dims$p);    #--number of post-molt age classes
  nZs = length(dims$z);    #--number of size classes
  nCs = nrow(dims$dmsN);   #--number of population categories (size of population vector)
  I_c  = double(nCs);                #--identity vector for population state
  I_cc = diag(x=1,nrow=nCs,ncol=nCs);#--identity matrix for population state

  #--set up model processes----
  ##--initial N's----
  nInit = getInitN(dims,opts,params,verbose);
  REPORT(nInit);

  if (FALSE)
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

    #--integrate model over years by season----
    n_ysc = array(nYs,nSs,nCs);#--population abundance at beginning season `s` in year `y`
    n_ysc[1,1,] = nInit;
    for (y_ in dims$y){ ##--loop over years----
      for (s_ in dims$s){ ###--loop over seasons within years----
        n_c = n_ysc[y_,s_];#--population state at start of y_, s_
        if (verbose) cat("y, s =",y_,s_,"\n");

        if (seasonLength[y_,s_]==0){  ####--instantaneous processes----
          #####---happen in order of
          #####----survey, growth, movement, and/or (possibly) fishery catches and mortality

          #####--surveys----
          for (v_ in 1:nSrvs){
            if (doSurveys_vys[v_,y_,s_]){
              srvN_vysc[[v_,y_,s_,]] = lstSQs[[v_]]$q_ysc[[y_,s_,]] * n_c;
              srvB_vysc[[v_,y_,s_,]] = lstWatz$wAtZ_ysc[[y_,s_,c_]] * srvN_vysc[[v_,y_,s_,]];
            }
          }#--v_

          #####--growth----
          if (doGrowth_ys[y_,s_]){
            n_c = lstGTMs$gtm_yscc[y_,s_,,] %*% n_c;
          }

          #####--instantaneous movement----
          if (doInstMovement[y_,s_]){
            n_c = lstMTMs$mtm_cc[y_,s_,,] %*% n_c;
          }

          #####--instantaneous fishery catches and mortality----
          #--calculate total fishing capture, retention, discard mortality, and total mortality rates
          ######--calculate total fishing mortality rates (FMR) by population component----
          totFMR_c = AD(vector("numeric",nCs));#--total FMRs
          for (f_ in 1:nFsh){
            if (doInstFisheries_fys[f_,y_,s_]){
              totFMR_c = totFMR_c + lstFMs[[f_]]$fmr_ysc[y_,s_,];#--add fishery-specific FMRs to total FMRs
            }
          }#--f_
          ######--calculate fishery-specific catch-related numbers----
          scl = (1.0-exp(-totFMR_c))/(totFMR_c+(1*(totFMR_c==0)));#--scaling factor (last bit excludes divides by 0)
          for (f_ in 1:nFsh){
            if (doInstFisheries_fys[f_,y_,s_]){
              fcn_f_ysc[f,y,s,] =  lstFMs[[f_]]$fcr_ysc[y_,s_,]*scl*n_c;#--captured numbers
              frn_f_ysc[f,y,s,] =  lstFMs[[f_]]$frr_ysc[y_,s_,]*scl*n_c;#--retained numbers
              fdn_f_ysc[f,y,s,] =  lstFMs[[f_]]$fdr_ysc[y_,s_,]*scl*n_c;#--discard mortality numbers
              fmn_f_ysc[f,y,s,] =  lstFMs[[f_]]$fmr_ysc[y_,s_,]*scl*n_c;#--retained + discard mortality numbers
            }
          }#--instantaneous fisheries
          n_c = n_c * exp(-totFMR_c);#--pop numbers after fisheries

        } else if (seasonLength[y_,s_]>0){ ####--continuous processes happen throughout season----

          #####--continuous process can include
          #####--natural mortality and continuous fishery catch and mortality
          ######--NOTE: Incorporating movement here requires exponentiating a matrix
          ######--to calculate changes to n_c. Not quite sure how to calculate fishery-specific
          ######--values (catch, mortality) in this case: is it directly analagous to the
          ######--case without movement??
          dt = seasonLength[y_,s_];

          #####--continuous movement rates----
          # if (doContMovement[y_,s_])
          #   mtm_cc = lstMTMs$mtm_yscc[y_,s_,,];#--c *x* c matrix

          #####--natural mortality rate----
          nmr_c = lstNMRs$nmr_ysc[y_,s_,];#---vector

          #####--continuous fishery catch rates
          ######--calculate total fishing mortality rates (FMR) by population component----
          totFMR_c = AD(vector("numeric",nCs));#--total FMRs (vector)
          for (f_ in 1:nFsh){
            if (doContFisheries_fys[f_,y_,s_]){
              totFMR_c = totFMR_c + lstFMs[[f_]]$fmr_ysc[y_,s_,];#--add fishery-specific FMRs to total FMRs
            }
          }#--f_
          ######--calculate fishery-specific catch-related numbers----
          ######--these may not be correct if continuous movement is included
          scl = (1.0-exp(-(totFMR_c+nmr_c)*dt))/(totFMR_c+nmr_c);#--scaling factor
          for (f_ in 1:nFsh){
            if (doContFisheries_fys[f_,y_,s_]){
              fcn_f_ysc[f,y,s,] =  lstFMs[[f_]]$fcr_ysc[y_,s_,]*scl*n_c;#--captured numbers
              frn_f_ysc[f,y,s,] =  lstFMs[[f_]]$frr_ysc[y_,s_,]*scl*n_c;#--retained numbers
              fdn_f_ysc[f,y,s,] =  lstFMs[[f_]]$fdr_ysc[y_,s_,]*scl*n_c;#--discard mortality numbers
              fmn_f_ysc[f,y,s,] =  lstFMs[[f_]]$fmr_ysc[y_,s_,]*scl*n_c;#--retained + discard mortality numbers
            }
          }#--continuous fisheries

          #####--update N for continuous processes----
          ######--if continuous movement is implemented,
          ######--the following needs to be converted to an exponentiated matrix (rather thn vector)
          n_c = exp(-(nmr_c+totFMR_c)*dt) * n_c;
        }#--if (seasonLength[y_,s_])

        ####--Finally, add recruitment----
        if (doRecruitment[y_,s_])
          n_c = n_c + lstRec$rec_ysc[y_,s_,];

        if (s_<nSs) {
          n_ysc[y_,s_+1,] = n_c;                    #--save pop state at start of y_, s_+1
          b_ysc[y_,s_+1,] = lstAlm$wAtZ[y_,s_,]*n_c;#--biomass at start of y_, s_+1
        }
      } #--end seasons loop (s_)

      if (y_<nYs) {
        n_ysc[y_+1,1,] = n_c;                    #--save pop state as start of y_+1, s_
        b_ysc[y_+1,1,] = lstAlm$wAtZ[y_,s_,]*n_c;#--biomass at start of y_, s_+1
      }
    } #--end years loop (y_)

    nFinal = n_c;                      #--save final pop state (start of year nYs+1, season 1)
    bFinal = lstAlm$wAtZ[nYs,nSs,]*n_c;#--final biomass
  }

  #--calculate likelihoods----

  if (testing){
    nll = -dnorm(1,dummy,1,log=TRUE);
  }
  if (verbose) cat("end objective function\n");
  return(nll);
}

