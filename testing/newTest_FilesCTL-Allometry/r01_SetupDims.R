#--set up model dimensions
#require(rtmbGMACS);
dirPrj = rstudioapi::getActiveProject();
source(file.path(dirPrj,"R/DimensionsFunctions.R"))
source(file.path(dirPrj,"R/DimensionsUtilities.R"))


vRs="EBS";
maz = as.character(seq(25,75,5));
fiz = as.character(seq(25,65,5));
fmz = as.character(seq(30,65,5));
vXs=list(male=list(imm=list(`new`=maz,
                             `old`=maz),
                   mat=list(`new`=maz,
                             `old`=maz)
                  ),
          female=list(imm=list(`new`=fiz,
                               `old`=fiz),
                      mat=list(`new`=fmz,
                               `old`=fmz)
                     )
         ); attr(vXs,"dmnms")<-c("x","m","p","z");
ys = as.character(2015:2024); names(ys) = ys;
ss = as.character(1:2);       names(ss) = ss;
dmsC    = createSparseDimsMap(r=vRs,x=vXs);            #--DimsMap for stock categories
dmsYS   = createSparseDimsMap(y=ys,s=ss);              #--DimsMap for years and seasons
dmsYSC  = createSparseDimsMap(y=ys,s=ss,r=vRs,x=vXs);  #--DimsMap for combination
dimsZBC = dmsC |> dplyr::mutate(lft=as.numeric(z)-2.5,mid=as.numeric(z),rgt=as.numeric(z)+2.5);
fleets = c("Directed","Bycatch","Survey"); #--fishery and survey fleets
