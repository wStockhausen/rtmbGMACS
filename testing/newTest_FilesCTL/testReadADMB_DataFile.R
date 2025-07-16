#--test readDataFile(...)
dirPrj =  rstudioapi::getActiveProject();
source(file.path(dirPrj,"R","MiscFunctions_ReadParse.R"));
source(file.path(dirPrj,"R","readADMB_DataFileFunctions.R"));
fn = file.path(dirPrj,"testing/inputFiles",
               "/TannerCrab/tcData_AllFisheriesNMFS-2SexCatchData.dat");
lstData = readADMB_DataFile(fn);
