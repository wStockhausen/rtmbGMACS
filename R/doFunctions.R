#--TODO: fill in details
#--function to do estimation
doEstimation<-function(mdl,inputs){
  opt = nlminb(mdl,gradient=mdl$gradient,hessian=mdl$hessian);
  #--set mdl's current par to best par?
  return(mdl)
}

#--function to do projections with optimized model
doProjections<-function(mdl,inputs){
  #--???
}

#--function to do MCMC with optimized model
doMCMC<-function(mdl){
  #--???
}
