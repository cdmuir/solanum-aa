# Write Stan models
source("r/header.R")

write_solanum_aa("aa1")
write_solanum_aa("aa2")
write_solanum_aa("aa3")
write_solanum_aa("aa4")
write_solanum_aa("aa5")

rstan:::rstudio_stanc("stan/solanum-aa1.stan")
rstan:::rstudio_stanc("stan/solanum-aa2.stan")
rstan:::rstudio_stanc("stan/solanum-aa3.stan")
rstan:::rstudio_stanc("stan/solanum-aa4.stan")
rstan:::rstudio_stanc("stan/solanum-aa5.stan")
