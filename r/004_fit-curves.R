# Fit quadratic parameters to each curve separately
source("r/header.R")

rh_curves = read_rds("data/trimmed_rh_curves.rds")

rh_curves |>
  split(~ curve) |>
  iwalk(\(df, curve_id) {
    x = 1
    ad = 0.8
    n_divergent = Inf
    while (n_divergent > 10) {
      fit_curve = brm(
        log_A ~ poly(log_gsw, 2, raw = TRUE) +
          arma(
            time = NA,
            gr = NA,
            p = 1,
            q = 0,
            cov = FALSE
          ),
        iter = x * 2000,
        thin = x,
        data = df,
        chains = 4,
        cores = 4,
        backend = "cmdstanr",
        control = list(adapt_delta = ad)
      )
      n_divergent = nuts_params(fit_curve) |>
        subset(Parameter == "divergent__") |>
        pull(Value) |>
        sum()
      x = x + 1
      ad = min(ad * 1.1, 0.99)
    }
    
    write_rds(fit_curve, paste0("objects/curve-fits/", curve_id, ".rds"))
    
  })
