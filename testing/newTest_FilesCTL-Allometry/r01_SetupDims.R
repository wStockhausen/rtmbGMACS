#--set up model dimensions to fit allometry data
##--since fitting to data here, will use it to create the dimensions
require(rtmbGMACS);
chkData=FALSE;

#--set up paths----
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"testing/newTest_FilesCTL-Allometry");

#--get allometric data----
dfrZW  = wtsUtilities::getObj(file.path(dirThs,"rda_AllometryData.RData")) |>
           dplyr::rename(obs=w) |> tibble::rownames_to_column(var="obs_id") |>
           dplyr::mutate(dplyr::across(dplyr::everything(),as.character));
if (chkData){
  plt = ggplot(dfrZW,aes(x=as.numeric(z),y=as.numeric(obs),colour=x)) + geom_point();
  print(plt);
  rm(plt)
}

##--remove some clearly bad values----
dfrZWp = dfrZW |>
          dplyr::filter(!((as.numeric(z)>125)&(as.numeric(obs)<100))) |>
          dplyr::filter(!((as.numeric(z)<100)&(as.numeric(obs)>400))) |>
          dplyr::filter(!((as.numeric(z)< 40)&(as.numeric(obs)> 75))) |>
          dplyr::filter(!(dplyr::between(as.numeric(z),125, 150)&dplyr::between(as.numeric(obs),250,300))) |>
          dplyr::filter(as.numeric(obs)>0);
if (chkData){
  plt = ggplot(dfrZWp,aes(x=as.numeric(z),y=as.numeric(obs),colour=x)) + geom_point() +
          coord_cartesian(xlim=c(0,175),ylim=c(0,NA))
  print(plt);
  rm(plt);
}

##--save censored data to dataframe----
dfrData = dfrZWp;
rm(dfrZW,dfrZWp);

##--calculate observed size limits by sex, maturity state, and post-molt age
dfrLims = dfrData |> dplyr::group_by(x,m,p) |>
            dplyr::summarize(minZ=min(z,na.rm=TRUE),
                             maxZ=max(z,na.rm=TRUE)) |>
            dplyr::ungroup();

