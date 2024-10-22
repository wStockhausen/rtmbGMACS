#--script to test objective function capabilities

#--can objfun handle a multilevel list as `data`? YES!!
#----save it to compiled model using REPORT(data)

#--can objfun handle a multilevel list as `params`? NO!!
#----`params` input to MakeADFun can be a named list, but
#----all elements must be vectors, matrices, or arrays
#----`arrays` are maximally-dimensioned (so `arrays` can't be "ragged")
#----"ragged" arrays in conventional R are implemented as lists

#--seems like one way to handle a ragged array would be to define it as a vector
#----and attach attributes that define it's "raggedness", along with functions that
#----allow one to get and set its values using the dimension indices
#------1. could implement as maximally-dimensioned arrays with NAs where appropriate
#------2. could implement as vector with index related back to ragged dimensions & vice-versa

alpha = 10;
beta  = 20;
x = 1:10;
data = list(x = x,
            y = alpha + beta*x + rnorm(length(x),0,1),
            opts=list(`1`="fizzy",`growth`="option 2"),
            ps=list(alpha=1,beta=2));
names(data$y) = paste0("x",as.character(x));
objfun<-function(params){
  #--expand the list of parameters
  getAll(data,params);

  if (data$opt[["1"]]=="fizzy") {
    z = 1;
  } else {
    z = 2;
  }
  REPORT(data) #--access via model$sumulate()$data
  REPORT(z)

  #--calculate objective function
  yp = p[ps$alpha] + p[ps$beta] * x;
  y %~% dnorm(yp,1);
  dy = y-yp;
  names(yp) = paste0("x",as.character(x));
  names(dy) = paste0("x",as.character(x));
  REPORT(y)
  REPORT(yp)  #--keeps names attributes
  REPORT(dy)  #--keeps names attributes

  return(.nll);
}

mdl = MakeADFun(objfun,list(p=c(1,2)));
rep = mdl$report();
opt = nlminb(mdl$par,mdl$fn,mdl$gr,hessian=mdl$he);

