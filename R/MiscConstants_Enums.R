#--"enums" (most based on glmmTMB/R/enum.R)----
##--.valid_link: based on glmmTMB/R/enum.R{.valid_link} - enum for valid link function for stat family definition----
.valid_link <- c(
  log                 = "log",       # 0
  logit               = "logit",     # 1
  probit              = "probit",    # 2
  inverse             = "inverse",   # 3
  cloglog             = "cloglog",   # 4
  identity            = "identity",  # 5
  sqrt                = "sqrt",      # 6
  lambertW            = "lambertW"   # 7
)
##--.valid_family: based on glmmTMB/R/enum.R{.valid_family} - enum for valid distributions functions----
.valid_family <- c(
  gaussian = "gaussian",                  #    0
  binomial = "binomial",                  #  100
  betabinomial ="betabinomial",           #  101
  beta         ="beta",                   #  200
  ordbeta = "ordbeta",                    #  201
  Gamma ="Gamma",                         #  300
  poisson ="poisson",                     #  400
  truncated_poisson ="truncated_poisson", #  401
  genpois ="genpois",                     #  402
  compois ="compois",                     #  403
  truncated_genpois ="truncated_genpois", #  404
  truncated_compois ="truncated_compois", #  405
  nbinom1 ="nbinom1",                     #  500
  nbinom2 ="nbinom2",                     #  501
  nbinom12 ="nbinom12",                   #  502
  truncated_nbinom1 ="truncated_nbinom1", #  550
  truncated_nbinom2 ="truncated_nbinom2", #  551
  t ="t",                                 #  600
  tweedie = "tweedie",                    #  700
  lognormal = "lognormal",                #  800
  skewnormal = "skewnormal",              #  900
  bell = "bell"                           # 1000
)
##--.valid_covstruct: based on glmmTMB/R/enum.R{.valid_covstruct} - enum for valid covariance structures----
.valid_covstruct <- c(
  diag    = "diag",     # 0,
  us      = "us",       # 1,
  cs      = "cs",       # 2,
  ar1     = "ar1",      # 3,
  ou      = "ou",       # 4,
  exp     = "exp",      # 5,
  gau     = "gau",      # 6,
  mat     = "mat",      # 7,
  toep    = "toep",     # 8,
  rr      = "rr",       # 9
  homdiag = "homdiag", # 10
  propto  = "propto",  # 11
  hetar1  = "hetar1",  # 12
  homcs   = "homcs",   # 13
  homtoep = "homtoep"  # 14
)
##--.valid_smooths: character vector for valid smooth functions----
.valid_smooths <- c("s",
                    "ti");

##--.valid_simcode: based on glmmTMB/R/enum.R{.valid_simcode} - enum for valid simulation options----
.valid_simcode <- c(
  zero = "zero",    # 0,
  fix = "fix",      # 1,
  random = "random" # 2
)

