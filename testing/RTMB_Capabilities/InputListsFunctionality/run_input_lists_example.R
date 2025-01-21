#--test functionality of input lists
require(RTMB);
require(ggplot2);

# How to:
#--1. map parameters to the same values using the "map" argument to MakeADFun
#--2. turn parameters on-and-off using the "map" argument to MakeADFun
#--3. set bounds on parameters by setting `upper` and `lower` in the optimization
#------when using the optimizer `nlminb`

like_type = "t"; #--"t" or "normal"

# Simulate artificial data from random effects model
set.seed(234)
n = 100         # Number of individuals
n_i = 4
n_j = 10
i = sample(1:n_i,size=n,replace=TRUE)   # Levels of factor 1
j = sample(1:n_j,size=n,replace=TRUE)   # Levels of factor 2
beta = rnorm(n_i,0,1)                   # Coefficients of factor 1
u = rnorm(n_j,0,1)                      # Random effects vector associated with 2

if (like_type=="normal"){
  d = rnorm(n,0,0.2);     #--devs
  y = beta[i] + u[j] + d; #--observations
} else
if (like_type=="t"){
  d = rt(n,n/2,0);
  y = beta[i] + u[j] + d; #--Data
}
dfr=tibble::tibble(x=1:n,y=y,d=d) |> tidyr::pivot_longer(c(y,d));
ggplot(dfr,aes(x=x,y=value,colour=name)) + geom_line();

#--create "inputs" (most RTMB examples call this "data") list
inputs <- list(data=list(y = y,
                         i = i,                   # 1-indexing in R
                         j = j),
               opts=list(like_type=like_type));

#--create parameters list, assigning initial values
parameters <- list(beta = rep(0,n_i),
                   u = rep(0,n_j),
                   sigma = (like_type=="normal")*1.0 +    #--std dev parameter for normal dist
                           (like_type=="t"     )*10.0,    #--dof for t
                   tau = 1);                              #--hyper std dev for u

#--source the objective function
compiler::enableJIT(0);
source("objective_function.R");
objective_function(parameters);

# Fit full model
obj <- MakeADFun(objective_function,parameters, random="u",silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj);
summary(rep);
summary(rep,"fixed");
summary(rep,"random");
summary(rep,"report");


