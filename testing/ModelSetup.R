#--model setup (see ModelParametersSetup.xlsx for ideas)----
##--(pseudo) cohort classes----
#----post-molt ages nested in maturity state/sex (working outward)
vXMPs = list(MALE=  list(IMMATURE=c("NEW SHELL"),
                         MATURE=c("NEW SHELL","OLD SHELL")),
             FEMALE=list(IMMATURE=c("NEW SHELL"),
                         MATURE=c("NEW SHELL","OLD SHELL")));
attr(vXMPs,"dmnms")<-c("x","m","p");
#----sizes nested in age/shell condition/maturity state/sex (working outward)
vXMZs = list( MALE=  list(IMMATURE=seq(40,60,5),
                          MATURE=seq(80,120,10)),
              FEMALE=list(IMMATURE=seq(25,50,5),
                          MATURE=seq(50,100,10)));
attr(vXMZs,"dmnms")<-c("x","m","z");
dfrMUS = createSparseDimsMap(p=vXMPs,z=vXMZs);

##--model time frame----
###--years----
vYs = 1982:2022;
attr(vYs,"dmnms") = "y";
###--seasons----
vSs = c("surveys","fall","fisheries","growth","spring","recruitment");
attr(vSs,"dmnms") = "s";
dfrTBMod = createSparseDimsMap(y=vYs,s=vSs);
####--season lengths----
vSs = c("surveys"=0.0,
        "fall"=0.625,
        "fisheries"=0.0,
        "mating"=0.0,
        "growth"=0.0,
        "spring"=0.375,
        "recruitment"=0.0);
####--season assignments for processes----
ssnPrcs=c(surveys     = "surveys",
          fisheries   = "fisheries",
          mating      = "mating",
          growth      = "growth",
          recruitment = "recruitment");

##--fleets----
fleets = c("TCF","NMFS");
fisheries = c("TCF");
infoFshs   = tibble::tibble(fleet="TCF",y=c(1982:1983,1987:1996,2005:2009,2013:2015,2017:2018,2020:2022),s="fisheries");
surveys   = c("NMFS");
infoSrvs   = tibble::tibble(fleet="NMFS",y=c(1982:2019,2021:2023),s="surveys");

##--time blocks----
TBs = list(`1`=vYs,
           `pre-rat`=1982:2004,
           `post-rat`=2005:max(vYs));

#--process/parameter specification----
##--initial abundance specification----
opts_InitN = list()

##--recruitment specification----


##--allometry----


##--growth----


##--natural mortality----
opts_LnM = list();
spcs_LnM = list(
  cats = c("x","m","TB"),
  process =
    tibble::tribble(
        ~idx,   ~x   ,    ~m    ,   ~TB  ,~base,  ~fcn, ~prType,   ~pr1   , ~pr2, ~peDist ,
          1 ,  "MALE","IMMATURE", "model",  0  ,"const","normal",log(0.23),  0.4, "normal",
          2 ,"FEMALE","IMMATURE", "model",  1  ,"mirr" ,   NA   ,   NA    ,  NA ,    NA   ,
          3 ,  "MALE",  "MATURE", "model",  1  ,"add" ,"normal",log(0.23) ,  0.4, "none"  ,
          4 ,"FEMALE",  "MATURE", "model",  1  ,"add" ,"normal",log(0.23) ,  0.4, "none"  ),
  pe_pars  = tibble::tribble(
      ~prc_idx,~idx, ~par , ~transf, ~ival , ~lb , ~ub, ~phs, ~`RE?`,
          1   ,  1 ,"mean", "none",   0   ,  -1 ,   1,   -1, FALSE,
          1   ,  2 ,"sd"  , "none",   1   , -10 ,  10,    5, FALSE),
  prc_lnks = tibble::tribble(
      ~prc_idx, ~idx,  ~TB  ,~eSeries,~ival, ~lb, ~ub, ~phs , ~prType, ~pr1,~pr2,~link,
          3   ,   1 ,"model", "WPDO" ,  0  , -5 ,  5 ,   4  , "normal",  0 ,  1 ,"add",
          3   ,   2 ,"model", "NPGO" ,  0  , -5 ,  5 ,   4  , "normal",  0 ,  1 ,"mult"
    ),
  params = tibble::tribble(
      ~prc_idx, ~idx , ~par ,   ~TB  , ~base, ~transf ,  ~ival  ,   ~lb   ,  ~ub   , ~phs , ~prType, ~pr1,~pr2, ~`RE?`,
          1,      1  ,  "p" , "model",   0  , "none"  ,log(0.23), log(0.1),log( 1.0),   4 ,"normal",  0  ,  1 , FALSE,
          3,      2  ,  "dp", "model",   0  , "none"  ,log(1.00), log(0.1),log(10.0),   4 ,"normal",  0  ,  1 , FALSE,
          4,      3  ,  "dp", "model",   0  , "none"  ,log(1.00), log(0.1),log(10.0),   4 ,"normal",  0  ,  1 , FALSE),
  par_lnks = tibble::tribble(),

);


##--selectivity----


##--fishing mortality----


##--survey catchability----


#--management-related quantities----


#--data and likelihoods----


#-projections----

