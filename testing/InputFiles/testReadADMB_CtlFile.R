#--test readCtlFile(...)
dirPrj =  rstudioapi::getActiveProject();
source(file.path(dirPrj,"R","readADMB_CtlFileFunctions.R"));
fn = file.path(dirPrj,"testing/inputFiles",
               "/TannerCrab/tcCTL_24_06.ctl");
lst = readADMB_CtlFile(fn);
