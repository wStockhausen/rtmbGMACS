#--examples for 3 dimensions, each with 5 levels
require(rtmbGMACS);
getFullDimsIndex(c(1,1,1),c(5,5,5));
getFullDimsIndex(c(5,1,1),c(5,5,5));
getFullDimsIndex(c(1,5,1),c(5,5,5));
getFullDimsIndex(c(1,1,5),c(5,5,5));
getFullDimsIndex(c(5,5,1),c(5,5,5));
getFullDimsIndex(c(1,5,5),c(5,5,5));
getFullDimsIndex(c(5,5,5),c(5,5,5));
getFullDimsIndex(c(1,1,1),c(5,5,5),TRUE);
getFullDimsIndex(c(5,1,1),c(5,5,5),TRUE);
getFullDimsIndex(c(1,5,1),c(5,5,5),TRUE);
getFullDimsIndex(c(1,1,5),c(5,5,5),TRUE);
getFullDimsIndex(c(5,5,1),c(5,5,5),TRUE);
getFullDimsIndex(c(1,5,5),c(5,5,5),TRUE);
getFullDimsIndex(c(5,5,5),c(5,5,5),TRUE);
