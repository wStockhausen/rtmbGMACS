#--get weight-at-size data for Tanner crab
require(tcsamSurveyData);
dirThs = dirname(rstudioapi::getSourceEditorContext()$path);

#--for Tanner crab"
{species <- "BTC"; col.Size<-"WIDTH";strata_type="crabpack";}

#----read strata file
fn = "~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current/TannerCrab_SurveyStrata.csv";
dfrStrata<-readr::read_csv(file.path(fn)) |>
             dplyr::select(STATION_ID,DISTRICT,TOWS,SURVEY_YEAR,LATITUDE,LONGITUDE,STRATUM,TOTAL_AREA_SQ_NM);
#----select appropriate strata
dfrSD<-selectStrata.TrawlSurvey(dfrStrata,
                                species   =species,
                                strataType=strata_type,
                                export    =FALSE,
                                verbosity=0);
rm(dfrStrata);

#--read in haul data file
fn = "~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current/TannerCrab_HaulData.csv";
dfrHaulData<-readr::read_csv(file=file.path(fn),
                             guess_max=Inf,
                             skip=5);
#----select haul data
dfrHD<-selectHauls.TrawlSurvey(dfrSD,
                               tbl=dfrHaulData,
                               YearRange=c(1975,2025),
                               export   =FALSE,
                               verbosity=0);
#----select indiv data
dfrID<-selectIndivs.TrawlSurvey(dfrHD,
                                tbl             =dfrHaulData,
                                sex             ="ALL",
                                maturity        ="ALL",
                                shell_condition ="ALL",
                                calcMaleMaturity=FALSE,
                                noImmOldShell   =TRUE,
                                dropLevels      =NULL,
                                minSize         =0,
                                maxSize         =Inf,
                                export          =FALSE,
                                verbosity       =0);
rm(dfrHaulData);

#--select size, weight data
dfrZW = dfrHD |> dplyr::inner_join(dfrID,by="HAULJOIN") |>
        dplyr::select(y=YEAR,x=SEX,m=MATURITY,p=SHELL_CONDITION,z=SIZE,w=WEIGHT) |>
        dplyr::filter(!is.na(w),(x %in% c("MALE","FEMALE"))) |>
        dplyr::mutate(r="EBS",
                      x=tolower(x),m=tolower(m),p=tolower(p),  #--convert to lower case
                      m=ifelse(m=="undetermined","all",m),     #--convert "undetermined" to "all"
                      m=stringr::str_extract(m,"..."),         #--keep only first 3 letters (e.g., 'imm','mat','all')
                      p=stringr::str_extract(p,"..."),
                      yb=ifelse(y==1975,"1975",">1975"));      #--distinguish 1975 from other years

require(ggplot2);
p = ggplot(dfrZW,aes(x=z,y=w,shape=yb,colour=m)) +
      geom_point() + facet_grid(x~.,scales="free_y") +
      labs(x="size (mm CW)",y="weight (g)",colour="maturity",shape="year") +
      wtsPlots::getStdTheme();
p;

#--apply some additional filtering based on above plot
dfrZWp = dfrZW |> dplyr::filter(!(y==1975),!((x=="female")&(w>600)));
p %+% dfrZWp;

wtsUtilities::saveObj(dfrZWp,file.path("rda_AllometryData.RData"));
readr::write_csv(dfrZWp,file.path(dirThs,"csv_DataAllometry.csv"));

