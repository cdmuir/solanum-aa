# Fit quadratic parameters to each curve separately
# Note: I removed ARMA because it was causing strange issues. Custom Stan model does it correctly.
source("r/header.R")

rh_curves = read_rds("data/trimmed_rh_curves.rds") |>
  mutate(ci = as.numeric(as.factor(curve)))

plan(multisession, workers = 19)
rh_curves |>
  split( ~ curve) |>
  future_iwalk(\(df, curve_id) {
    if (!file.exists(paste0("objects/curve-fits/", curve_id, ".rds"))) {
      x = 2
      ad = 0.8
      n_divergent = Inf
      while (n_divergent > 10) {
        fit_curve = brm(
          log_A ~ poly(log_gsw, 2, raw = TRUE),
          iter = x * 2000,
          thin = x,
          data = df,
          chains = 4,
          cores = 4,
          backend = "cmdstanr",
          control = list(adapt_delta = ad),
          seed = 360036340 + df$ci[1]
        )
        n_divergent = nuts_params(fit_curve) |>
          subset(Parameter == "divergent__") |>
          pull(Value) |>
          sum()
        x = x + 1
        ad = min(ad * 1.1, 0.99)
      }
      
      write_rds(fit_curve,
                paste0("objects/curve-fits/", curve_id, ".rds"))
    } else {
      message(paste0("Curve fit for ", curve_id, " already exists. Skipping."))
    }
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))
