#--based on Hans Skaug's TMB example at https://github.com/skaug/tmb-case-studies/tree/master/map_example
require(RTMB);

# How to:
#--1. map parameters to the same values using the "map" argument to MakeADFun
#--2. turn parameters on-and-off using the "map" argument to MakeADFun
#--3. set bounds on parameters by setting `upper` and `lower` in the optimization
#------when using the optimizer `nlminb`


# Simulate artificial data from random effects model
set.seed(234)
n = 100         # Number of individuals
n_i = 4
n_j = 10
i = sample(1:n_i,size=n,replace=TRUE)   # Levels of factor 1
j = sample(1:n_j,size=n,replace=TRUE)   # Levels of factor 2
beta = rnorm(n_i,0,1)                   # Coefficients of factor 1
u = rnorm(n_j,0,1)                      # Random effects vector associated with 2

y = beta[i] + u[j] + rnorm(n,0,.2)         # Data

#--create "inputs" (most RTMB examples call this "data") list
inputs <- list(y = y,
               i = i,                   # 1-indexing in R
               j = j)

#--create parameters list, assigning initial values
parameters <- list(beta = rep(0,n_i),
                   u = rep(0,n_j),
                   sigma = 1,
                   tau = 1);

#--source the objective function
source("objective_function.R");

# Fit full model
obj <- MakeADFun(objective_function,parameters, random="u",silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
print(sdreport(obj));

# Make all beta's equal
map=list(beta=as.factor(c(1,1,1,1)));#--creating vector of factors with 1 level
obj <- MakeADFun(objective_function, parameters, random="u",silent=TRUE,map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)
print(sdreport(obj));

# Group beta's in two groups
map=list(beta=as.factor(c(1,1,2,2)))
obj <- MakeADFun(objective_function, parameters, random="u",silent=TRUE, map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)
print(sdreport(obj));

# Do not estimate beta[1] (fix at initial value 0)
map=list(beta=factor(c(NA,2,3,4)));#--NA's must be assigned when creating factor
obj <- MakeADFun(objective_function, parameters, random="u",silent=TRUE, map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)
print(sdreport(obj));

# Do not estimate any beta's (fix at initial value 0)
map=list(beta=as.factor(c(NA,NA,NA,NA)))
obj <- MakeADFun(objective_function, parameters, random="u",silent=TRUE, map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)
print(sdreport(obj));

# Remove random effect from model
map=list(u=as.factor(rep(NA,length(u))),tau=as.factor(NA))
obj <- MakeADFun(objective_function, parameters, random="u",silent=TRUE, map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)
print(sdreport(obj)); print(opt$objective);

# Bounds on fixed effects and hyper parameters (not "u")
map=list()               #--no "mirroring"
U = c(rep(10,n_i),10,10)
names(U) = c(rep("beta",n_i),"sigma","tau")
L = c(rep(-10,n_i),1e-10,1e-10)
names(L) = names(U)
obj <- MakeADFun(objective_function, parameters, random="u",silent=TRUE, map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr,lower=L,upper=U)
print(sdreport(obj)); print(opt$objective);

# Model without random effect but with parameter bounds
map=list(u=as.factor(rep(NA,length(u))),tau=as.factor(NA)) #--remove all u's and hyper-parameter tau
obj <- MakeADFun(objective_function, parameters, random="u",silent=TRUE, map=map)
print(names(U))  # To see which component of L and U should be removed
idx = which(names(U)=="tau");
opt <- nlminb(obj$par, obj$fn, obj$gr,lower=L[-6],upper=U[-6])
print(sdreport(obj)); print(opt$objective);
