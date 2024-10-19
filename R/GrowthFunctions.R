#--growth functions
calcMnPostMoltZ1<-function(z,params,consts=NULL,debug=FALSE){
  grA = params[1];
  grB = params[2];
  mnZs = mfexp(grA+grB*log(z));
  names(mnZs) = as.character(z);
  return(mnZs);
}
calcMnPostMoltZ2<-function(z,params,consts,debug=FALSE){
  grA = params[1];
  grB = params[2];
  zGrA = consts[1]; #--pre-molt size corresponding to grA as mean post-molt size
  zGrB = consts[2]; #--pre-molt size corresponding to grB as mean post-molt size
  mnZs = grA*mfexp(log(grB/grA)/log(zGrB/zGrA)*log(z/zGrA));
  names(mnZs) = as.character(z);
  return(mnZs);
}
calcMnPostMoltZ3<-function(z,params,consts,debug=FALSE){
  grA = params[1];
  grB = params[2];
  zGrA = consts[1]; #--pre-molt size corresponding to grA as mean post-molt size
  mnZs = grA*mfexp(grB*log(z/zGrA));
  names(mnZs) = as.character(z);
  return(mnZs);
}
