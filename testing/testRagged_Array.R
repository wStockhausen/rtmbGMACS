#--test ragged_array S3 class

dirPrj = rstudioapi::getActiveProject();
source(file.path(dirPrj,"R","DimensionsUtilities.R"),local=TRUE);
source(file.path(dirPrj,"R","DimensionsFunctions.R"),local=TRUE);
source(file.path(dirPrj,"R","ragged_array.R"),local=TRUE);

#--define dimensions
dZs = list( MALE=  list(IMMATURE=list(`NEW SHELL`=seq(25,180,5),
                                      `OLD SHELL`=seq(25,180,5)),
                          MATURE=list(`NEW SHELL`=seq(25,180,5),
                                      `OLD SHELL`=seq(25,180,5))),
          FEMALE=list(IMMATURE=list(`NEW SHELL`=seq(25,130,5),
                                    `OLD SHELL`=seq(25,130,5)),
                        MATURE=list(`NEW SHELL`=seq(25,130,5),
                                    `OLD SHELL`=seq(25,130,5))));
attr(dZs,"dmnms")<-c("x","m","s","z");
dfrDZs = createSparseDimsMap(z=dZs);

dfrDYs = createSparseDimsMap(y=c(`2021`=2021,`2022`=2022));

#--extract subset
dfrDZsp = dfrDZs |> dplyr::filter(x=="MALE");
dfrDZsp1 =  tibble::tibble(z=factor(c(25,50),levels=levels(dfrDZsp$z)));

#--create vector for testing: must be same size as dimensions object
x  = 10000*as.numeric(as.character(dfrDZs$z)) + dfrDZs$sparse_idx;

#--test ragged_vector creation and subsetting
y   = ragged_array(x,dfrDims=dfrDZs); #--convert a vector to a ragged array
z   = ragged_array(y);                #--copy a ragged_array
zp  = ragged_array(z,dfrDZsp);        #--subset a ragged array using a subset of its dimensions
zp1 = ragged_array(z,dfrDZsp1);       #--subset a ragged array using a subset of its dimensions
y[1];     #--subset a ragged array by index
as.vector(zp);
zpp=zp[20:30]; #--subset a ragged array by vector of index

#--convert a ragged_array to a tibble
#----TODO: include vector name somehow?
tblZP  = as_tibble.ragged_array(zp);
tblZP1 = as_tibble.ragged_array(zp1);
print(tblZP1);

#--test expanding a vector or ragged_array
yy = expand(zp1,dfrDYs);
as_tibble.ragged_array(yy);
str(yy);


