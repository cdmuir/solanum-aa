# testing out fitting AA model over posterior
source("r/header.R")

aa1 = cmdstan_model("stan/aa1.stan", dir = "stan/bin")

stan_data = read_rds("data/stan_data.rds")
n_draw = length(stan_data)
names(stan_data) = str_pad(seq_len(n_draw), floor(log10(n_draw)) + 1, "left", "0")

stan_data[687:4000] |>
  iwalk(\(ll, draw_id) {

    # need to fix par_table
    # focal_pars = get_par_table("aa1") |>
    #   dplyr::filter(type == "real") |>
    #   pull(parameter)
    focal_pars = c(
      "b0_aa",
      "b_aa_light_intensity_2000",
      "b_aa_light_treatment_high",
      "rhosq_aa_acc",
      "etasq_aa_acc",
      "sigma_aa_acc_id",
      "b0_log_sigma_aa" ,
      "b_log_sigma_aa_light_intensity_2000",
      "b_log_sigma_aa_light_treatment_high"
    )
    
    x = 5
    ad = 0.8
    n_divergent = Inf
    converged = FALSE
    while ((n_divergent > 10 | !converged) & (x < 11)) {
      fit_aa1 = aa1$sample(
        data = ll,
        seed = 987587829 + as.numeric(draw_id) + x,
        refresh = x * 2000 * 0.1,
        chains = 4L,
        parallel_chains = 4L,
        iter_warmup = x * 1000,
        iter_sampling = x * 1000,
        thin = x,
        adapt_delta = ad
      )
      n_divergent = nuts_params(fit_aa1) |>
        subset(Parameter == "divergent__") |>
        pull(Value) |>
        sum()
      s = fit_aa1$summary(focal_pars)
      converged = all(s$rhat < 1.01)
      x = x + 1
      # x = 2 * x
      ad = min(ad * 1.1, 0.99)
    }
    
    if (n_divergent <= 10 & converged) {
      fit_aa1$save_object(paste0("objects/fit-aa1/", draw_id, ".rds"))
    }
    
  })



# OLD. delete once no longer useful ----
# args = commandArgs(trailingOnly = TRUE)
args = list(1)

stan_data = read_rds("data/stan_data.rds")

aa1 = cmdstan_model("stan/aa1.stan", dir = "stan/bin")

system.time({
  fit1 = aa1$sample(
  data = stan_data[[args[[1]]]],
  chains = 1L,
  parallel_chains = 1L,
  seed = 20240203,
  iter_warmup = 2e3,
  iter_sampling = 2e3,
  refresh = 2e2
)
})

# fit1$summary(c("rhosq_aa_acc", "etasq_aa_acc"))
# fit1$summary(c("rhosq_ppfd_aa", "etasq_ppfd_aa"))

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

# exploring how brms does fit to emualte
# 58 seconds
aa_post = read_rds("objects/aa_post.rds") |>
  filter(.draw <= 100) |>
  split(~ .draw)

library(future)
plan(multisession, workers = 10)
fit1 = brm_multiple(
  bf(
    aa ~
      light_intensity * light_treatment +
      (light_intensity * light_treatment | acc) +
      (1 | acc_id),
    sigma ~ light_intensity * light_treatment
  ),
  data = aa_post,
  backend = "cmdstanr",
  chains = 1,
  silent = 0,
  iter = 8000, thin = 4,
  control = list(adapt_delta = 0.99)
)

summary(fit1)
