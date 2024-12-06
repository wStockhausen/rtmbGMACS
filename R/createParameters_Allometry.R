#--create parameters for allometry
createParameters_Allometry<-function(lst){
  mpPars = lst$dfrMPs$

}
# #--test it
dirPrj = rstudioapi::getActiveProject();
conn=file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.function.txt");
res0 = readParamInfo_Allometry(conn,verbose=FALSE);
res1 = extractParameters_Allometry(res0,dims$dms_yrxmaz,verbose=FALSE);
res2 = createParameters_Allometry(res1)
