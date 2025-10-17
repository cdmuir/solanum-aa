# Test whether LMA mediates affect growth light treatment on AA (steady state)
source("r/header.R")

d1 = read_rds("objects/stan_data_df_steadystate.rds")
tr = read_rds("data/phylogeny.rds")
A = vcv(tr, corr = TRUE)

fixed_effects = c(
  "1",
  "light_treatment",
  "light_intensity",
  "mi(llma)",
  "light_treatment + light_intensity",
  "light_treatment * light_intensity",
  "light_treatment + mi(llma)",
  "light_treatment * mi(llma)",
  "light_intensity + mi(llma)",
  "light_intensity * mi(llma)",
  "light_treatment + light_intensity + mi(llma)",
  "light_treatment * light_intensity + mi(llma)",
  "light_treatment * mi(llma) + light_intensity",
  "light_treatment + light_intensity * mi(llma)",
  "light_treatment * light_intensity * mi(llma)"
)

random_effects = c(
  "(1 | gr(acc, cov = A))",
  "(1 + light_treatment | gr(acc, cov = A))",
  "(1 + light_intensity | gr(acc, cov = A))",
  "(1 + light_treatment + light_intensity | gr(acc, cov = A))",
  "(1 + light_treatment * light_intensity | gr(acc, cov = A))"
)

sigma_effects = c("light_treatment + light_intensity")

# Cross all combinations
set.seed(602275201)
model_forms = expand.grid(
  fixed = fixed_effects,
  random = random_effects,
  sigma = sigma_effects,
  stringsAsFactors = FALSE
) %>%
  mutate(seed = sample(1e9, nrow(.)),
         model = paste0("model_", row_number()))

# Build and fit each model
plan(multisession, workers = 19)

aa_models = model_forms |>
  dplyr::select(fixed, random, sigma, seed) |>
  future_pmap(function(fixed, random, sigma, seed) {
    fml = bf(as.formula(paste(
      "aa | se(se_aa, sigma = TRUE) ~", fixed, "+", random
    )), as.formula(paste("sigma ~", sigma)))
    
    crit = 0
    x = 1

    while (crit == 0 & x < 24) {
      m = brm(
        formula = fml + bf(llma |
                             mi() ~ light_treatment + (1 + light_treatment |
                                                         gr(acc, cov = A))),
        data = d1,
        data2 = list(A = A),
        backend = "cmdstanr",
        chains = 1,
        iter = x * aa_args$n_iter_init,
        thin = x,
        family = student(),
        save_pars = save_pars(all = TRUE),
        seed = seed
      )
      
      crit = summarise_draws(m) |>
        summarize(c1 = (max(rhat, na.rm = TRUE) < aa_args$max_rhat),
                  c2 = (min(ess_bulk, na.rm = TRUE) > aa_args$min_ess)) |>
        mutate(
          n_divergent = nuts_params(m) |>
            subset(Parameter == "divergent__") |>
            pull(Value) |>
            sum(),
          c3 = (n_divergent < aa_args$max_divergent),
          c4 = c1 * c2 * c3
        ) |>
        pull(c4)
      
      x = x + 1

    }
    
    return(m)
    
  },
  .progress = TRUE,
  .options = furrr_options(seed = TRUE))

# Extract loo objects
loos = future_map(
  aa_models,
  \(m) loo(
    m,
    resp = "aa",
    moment_match = FALSE,
    reloo = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

# Name the models
names(loos) = paste0("model_", seq_along(loos))

# Create table
aa_loo_table = tibble(
  model = names(loos),
  looic = sapply(loos, \(x) x$estimates["looic", "Estimate"]),
  se = sapply(loos, \(x) x$estimates["looic", "SE"])
) |>
  arrange(looic) |>
  mutate(
    delta_looic = looic - min(looic),
    best_model = delta_looic == 0
  ) |>
  left_join(model_forms, by = join_by(model))

write_rds(aa_loo_table, "objects/aa_sty_loo_table2.rds")
best_model_index = str_extract(aa_loo_table[1, "model"], "\\d+") |>
  as.integer()
write_rds(aa_models[[best_model_index]], "objects/fit_aa2_sty.rds")
