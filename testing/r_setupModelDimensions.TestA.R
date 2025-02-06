#--code to create model dimensions for testing----
require(rtmbGMACS);

##--set up model dimensions----
###--temporal dims  : "y","s"----
###--population dims: "r","x","m","a","p","z"----
###--fleet dims     : "f"----
##--dimensions maps----
###--dmsYS : years *x* seasons
###--dmsN : population state (r *x* x *x* m *x* a *x* p *x* z)
###--dmsF : fleets
setupModelDims<-function(){
  dims = list();
  dims$y = 2020:2024;                           #--years
  dims$s = 1:4;                                 #--seasons
  dims$r = "EBS";                               #--regions
  dims$x = c("male","female");                  #--sex classes
  dims$m = c("immature","mature");              #--maturity state classes
  dims$a = "all";                               #--post-recruitment age classes
  dims$p = c("newshell","oldshell");            #--post-molt ages
  dims$zc = seq(24.5,184.5,5);
  zc = dims$zc; nzcs = length(zc);              #--size bin cutpoints
  dims$zb =  0.5*(zc[2:nzcs]+zc[1:(nzcs-1)]);   #--size bin centers
  dims$z =  dims$zb;                            #--size classes (bin centers)
  dims$f = c("TCF","NMFS");                     #--fleets
  for (dim in names(dims)){
    dimp = dims[[dim]];
    names(dimp) = as.character(dimp);
    dims[[dim]] = dimp;
  }

  dims$dmsYS   = createSparseDimsMap(y=dims$y,s=dims$s);
  dims$dmsN    = createSparseDimsMap(r=dims$r,x=dims$x,m=dims$m,a=dims$a,p=dims$p,z=dims$z);
  dims$dmsYSN  = createSparseDimsMap(y=dims$y,s=dims$s,r=dims$r,x=dims$x,m=dims$m,a=dims$a,p=dims$p,z=dims$z);
  dims$dmsFYSN = createSparseDimsMap(f=dims$f,y=dims$y,s=dims$s,r=dims$r,x=dims$x,m=dims$m,a=dims$a,p=dims$p,z=dims$z);
  return(dims);
}
#--dims = setupModelDims();
