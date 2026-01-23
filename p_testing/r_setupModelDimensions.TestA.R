#--code to create model dimensions for testing----
require(rtmbGMACS);

##--set up model dimensions----
###--temporal dims        : "y","s"----
###--population categories: "r","x","m","p","z"----
###--fleet dims           : "f"----
##--dimensions maps----
###--dmsC    : population categories (r *x* x *x* m *x* p *x* z)
###--dmsYS   : years *x* seasons
###--dmsYSC  : years *x* seasons *x* pop categories
###--dmsFYSC : fleets *x* years *x* seasons *x* pop categories
setupModelDims<-function(zcs=seq(55.5,104.5,5)){
  nzcs = length(zcs);

  dims = list();
  dims$y = 2020:2024;                           #--years
  dims$s = 1;                                   #--seasons
  dims$r = "EBS";                               #--regions
  dims$x = c("male","female");                  #--sex classes
  dims$m = c("immature","mature");              #--maturity state classes
  dims$p = c("new_shell","old_shell");          #--post-molt ages
  dims$zc = zcs;                                #--size bin cutpoints
  dims$zb =  0.5*(zcs[2:nzcs]+zcs[1:(nzcs-1)]); #--size bin centers
  dims$z =  dims$zb;                            #--size classes (bin centers, repeat of zb)
  dims$f = c("TCF","SCF","NMFS");               #--fleets
  for (dim in names(dims)){
    dimp = dims[[dim]];
    names(dimp) = as.character(dimp);
    dims[[dim]] = dimp;
  }
  rm(nzcs,zcs,dim,dimp);

  dims$dmsC    = createSparseDimsMap(r=dims$r,x=dims$x,m=dims$m,p=dims$p,z=dims$z);
  dims$dmsYS   = createSparseDimsMap(y=dims$y,s=dims$s);
  dims$dmsYSC  = createSparseDimsMap(y=dims$y,s=dims$s,r=dims$r,x=dims$x,m=dims$m,p=dims$p,z=dims$z);
  dims$dmsFYSC = createSparseDimsMap(f=dims$f,y=dims$y,s=dims$s,r=dims$r,x=dims$x,m=dims$m,p=dims$p,z=dims$z);

  dims$nCs = nrow(dims$dmsC); #--number of population categories
  dims$nFs = length(dims$f);  #--number of fleets
  dims$nYs = length(dims$y);  #--number of model years
  dims$nSs = length(dims$s);  #--number of model seasons
  return(dims);
}
if (FALSE){
  dims = setupModelDims();
  opts = list(trackCohorts=FALSE,
              maxCohortAge=10);
}
