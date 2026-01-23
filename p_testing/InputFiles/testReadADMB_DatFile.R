#--test readPrjFile(...)
dirPrj =  rstudioapi::getActiveProject();
source(file.path(dirPrj,"R","MiscFunctions_ReadParse.R"));
source(file.path(dirPrj,"R","readADMB_DatFile.R"));
fn = file.path(dirPrj,"testing/inputFiles",
               "/TannerCrab/tcDat_24_06.dat");
lst = readADMB_DatFile(fn);
