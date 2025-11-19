#--get molt increment data for Tanner crab
require(ggplot2);
dirThs = dirname(rstudioapi::getSourceEditorContext()$path);

#--read current premolt, postmolt data
fn = "~/Work/StockAssessments-Crab/Data/MoltIncrementData/CurrentData/TannerCrab.AllGrowthData.csv";
dfrMID = readr::read_csv(fn) |>
           dplyr::filter(region=="EBS") |>
           dplyr::select(x=sex,y=year,TM=`TM?`,z=premolt,obs=postmolt) |>
           dplyr::mutate(y=as.character(y),
                         x=ifelse(x==1,"male",ifelse(x==2,"female",NA_character_)));


ggplot(dfrMID,aes(x=z,y=obs,colour=x)) + geom_point();

wtsUtilities::saveObj(dfrMID,file.path(dirThs,"rda_GrowthData.RData"));
readr::write_csv(dfrMID,file.path(dirThs,"csv_DataGrowth.csv"));

