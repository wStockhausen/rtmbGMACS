% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DimensionsFunctions.R
\name{createIntrinsicDims}
\alias{createIntrinsicDims}
\title{Create the intrinsic model dimensions from the user model dimensions}
\usage{
createIntrinsicDims(udfr)
}
\arguments{
\item{udfr}{- dataframe representing user dimensions map (sparse or full)}
}
\value{
a dataframe representing the intrinsic model dimensions.
}
\description{
Function to create the intrinsic model dimensions from the user model dimension.
}
\details{
The user dimensions defined by `udfr` are expanded to the full set of
intrinsic model dimensions ("y","s","r","x","m","a","p","z"), with default values
("all") added for intrinsic dimensions undefined in the user dimensions.
   * `y` - year index
   * `s` - season index
   * `r` - region index
   * `x` - sex index
   * `m` - maturity state index
   * `a` - age index
   * `p` - post-molt age (or shell condition) index
   * `z` - size index
}
\examples{
#--create user dimensions----
require(tibble);
#----years
vYs = 2001:2005;                  attr(vYs,"dmnms")<-"y";
#----seasons
vSs = c("fall","spring");         attr(vSs,"dmnms")<-"s";
#----regions
vRs = c("All EBS");               attr(vRs,"dmnms")<-"r";
#----sexes
vXs = c("MALE","FEMALE");         attr(vXs,"dmnms")<-"x";
#----post-molt ages nested in maturity state/sex (working outward)
vPs = list(MALE=  list(IMMATURE=c("NEW SHELL"),
                        MATURE=c("NEW SHELL","OLD SHELL")),
          FEMALE=list(IMMATURE=c("NEW SHELL"),
                        MATURE=c("NEW SHELL","OLD SHELL")));
attr(vPs,"dmnms")<-c("x","m","p");
#----sizes nested in age/shell condition/maturity state/sex (working outward)
vZs = list( MALE=  list(IMMATURE=seq(40,60,5),
                          MATURE=seq(80,120,10)),
          FEMALE=list(IMMATURE=seq(25,50,5),
                        MATURE=seq(50,100,10)));
attr(vZs,"dmnms")<-c("x","m","z");

#--create maps from the 1-d index to the dimension levels for various combinations
#----of the example dimensions.
require(rtmbGMACS);
dfrSparse = createSparseDimsMap(y=vYs,s=vSs,r=vRs,x=vXs,p=vPs,z=vZs);
dfrDense  = createDenseDimsMap(dfrSparse);

#--create intrinsic dimensions
dfrSprsI = createIntrinsicDims(dfrSparse);
dfrDensI = createIntrinsicDims(dfrDense);

}
