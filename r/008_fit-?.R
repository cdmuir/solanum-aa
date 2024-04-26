# testing out fitting AA model over posterior
source("r/header.R")

# args = commandArgs(trailingOnly = TRUE)
args = list(1)

stan_data = read_rds("data/stan_data.rds")

aa1 = cmdstan_model("stan/aa1.stan", dir = "stan/bin")

fit1 = aa1$sample(
  data = stan_data[[args[[1]]]],
  chains = 1L,
  parallel_chains = 1L,
  seed = 20240203,
  iter_warmup = 2e3,
  iter_sampling = 2e3,
  refresh = 2e2
)

fit1$summary(c("rhosq_aa_acc", "etasq_aa_acc"))
fit1$summary(c("rhosq_ppfd_aa", "etasq_ppfd_aa"))

# with phylo model
# default settings
# 2e3 iter, all hit max treedepth
# 
# with nonphylo model
# default settings
# 2e3 iter, all hit max treedepth. why?
# 
# 
# 
# fit_sim$save_object(glue("objects/fit_sim{n}.rds"))
