#
dirPrj = rstudioapi::getActiveProject();
source(file.path(dirPrj,"R","DimensionsUtilities.R"),local=TRUE);
source(file.path(dirPrj,"R","DimensionsFunctions.R"),local=TRUE);
source(file.path(dirPrj,"R","ragged_array.R"),local=TRUE);

#--define some dimensions
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
dfrDZsp  = dfrDZs |> dplyr::filter(x=="MALE"); #<- still has dfrDZs attributes
dfrDZsp1 = createSparseDimsMap(z=factor(c(25,50),levels=levels(dfrDZsp$z)));

#--want to model value for function as v1 = mean + p_x across dfrDZs
count_unique<-function(x){return(length(unique(x)));}
count_params<-function(dfrMM){
  sub1<-function(x){x-1}
  uP = (dplyr::distinct(dfrMM) |>
         dplyr::summarize(dplyr::across(tidyselect::everything(),count_unique)) |>
         dplyr::mutate(dplyr::across(tidyselect::everything(),sub1)) |>
         dplyr::rowwise() |> dplyr::transmute(s=sum(dplyr::c_across(tidyselect::everything()))))[[1]];
  if ("`(Intercept)`" %in% names(dfrMM)) uP = uP-1;
  return(uP);
}
mm_v1 = tibble::as_tibble(model.matrix(~1+x,data=dfrDZs));
np_v1 = ncol(mm_v1);

dTYs = list(TB1 = c(2005:2009),
            TB2 = c(2010:2014));
attr(dTYs,"dmnms")<-c("tb","y");
dfrTYs = createSparseDimsMap(tb=dTYs);
mm_v2 = tibble::as_tibble(model.matrix(~tb+y,data=dfrTYs));
np_v2 = ncol(mm_v1);

#--growth parameterization
lst = list(TB1=2015:2020,
           TB2=2021:2024);
attr(lst,"dmnms")<-c("tb","y");
dfrGrTB = createSparseDimsMap(y=lst);
lst = list(  MALE="IMMATURE",
           FEMALE="IMMATURE");
attr(lst,"dmnms")<-c("x","m");
dfrGrXM = createSparseDimsMap(m=lst);
