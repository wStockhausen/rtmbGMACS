#--test readPrjFile(...)
dirPrj =  rstudioapi::getActiveProject();
source(file.path(dirPrj,"R","readADMB_PrjFile.R"));
fn = file.path(dirPrj,"testing/inputFiles",
               "/TannerCrab/tcPrj_24_06.prj");
lst = readADMB_PrjFile(fn,nFlts=9);
