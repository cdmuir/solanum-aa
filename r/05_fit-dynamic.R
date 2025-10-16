# Fit quadratic parameters to each curve separately (dynamic A)
source("r/header.R")

rh_curves = read_rds("data/trimmed_rh_curves.rds") |>
  mutate(ci = as.numeric(as.factor(curve)))

dir = "objects/curve-fits/dynamic/"
plan(multisession, workers = 19)

rh_curves |>
  split(~ curve) |>
  future_iwalk(\(df, curve_id) {
    s = paste0(dir, curve_id, ".rds")
    if (!file.exists(s)) {
      crit = 0
      x = 1
      ad = 0.8
      
      while (crit == 0 & x < 24) {
        m = brm(
          log_Adyn ~ poly(log_gsw, 2, raw = TRUE),
          iter = x * aa_args$n_iter_init,
          thin = x,
          data = df,
          chains = 1,
          cores = 1,
          backend = "cmdstanr",
          control = list(adapt_delta = ad),
          seed = 360036340 + df$ci[1]
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
        ad = min(ad * 1.1, 0.99)
        
      }
      
      write_rds(m, s)
    } else {
      message(paste0("Curve fit for ", curve_id, " already exists. Skipping."))
    }
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))

fit_dynamic_diagnostics = list.files(dir, full.names = TRUE) |>
  future_map_dfr(\(.x) {
    m = read_rds(.x)
    s = summarise_draws(m) |>
      filter(variable != "lprior")
    tibble(
      file = .x,
      r2 = bayes_R2(m)[, "Estimate"],
      max_rhat = max(s$rhat),
      min_ess = min(s$ess_bulk),
      n_divergent = nuts_params(m) |>
        subset(Parameter == "divergent__") |>
        pull(Value) |>
        sum()
    )
  }, .progress = TRUE)

write_rds(fit_dynamic_diagnostics, "objects/fit_dynamic_diagnostics.rds")
