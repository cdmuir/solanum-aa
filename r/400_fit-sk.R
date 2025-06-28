# Fit CDWeibulll mode to each curve separately
source("r/header.R")

sk_dir = "objects/sk-curves/"
if (!dir.exists(sk_dir)) {
  dir.create(sk_dir)
}

rh_curves = read_rds("data/trimmed_rh_curves.rds") |>
  mutate(ci = as.numeric(as.factor(curve))) |>
  mutate(t_sec = elapsed - min(elapsed), .by = ci)
  
plan(multisession, workers = 2)
rh_curves |>
  split( ~ curve) |>
  magrittr::extract(seq_len(4)) |>
  future_iwalk(\(df, curve_id) {
    file = paste0(sk_dir, curve_id, ".rds")
    if (!file.exists(file)) {
      
      # set prior?
      
      x = 2
      ad = 0.8
      n_divergent = Inf
      
      while (n_divergent > 10) {
        fit_curve = brm(
          bf(
            gsw ~ gf + dg * exp(-(t_sec / tau)^lambda),
            gf ~ 1,
            dg ~ 1,
            tau ~ 1,
            lambda ~ 1,
            nl = TRUE
          ),
          iter = x * 2000,
          thin = x,
          data = df,
          chains = 1,
          cores = 1,
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
      
      write_rds(fit_curve, file)
    } else {
      message(paste0("Curve fit for ", curve_id, " already exists. Skipping."))
    }
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))



fit1 = read_rds(paste0(sk_dir, "LA0107-C_pseudohypo_150.rds")) 

ce = conditional_effects(fit1) 

ggplot(ce[[1]], aes(t_sec)) +
  geom_lineribbon(mapping = aes(y = estimate__, ymin = lower__, ymax = upper__),
                  color = "steelblue", fill = "grey") +
  geom_point(data = fit1$data, mapping = aes(y = gsw))


prior1 <- prior(normal(0.2251668, 0.01), nlpar = "gf") +
  prior(normal(-0.1514164, 0.02), nlpar = "dg", ub = 0) +
  prior(normal(600, 600), nlpar = "tau", lb = 0)

fit1 = brm(
  bf(
    gsw ~ gf + dg * exp(-(t_sec / tau)^lambda),
    gf ~ 1 + (1 | curve),
    dg ~ 1 + (1 | curve),
    tau ~ 1 + (1 | curve),
    lambda ~ 1 + (1 | curve),
    # autocor = ~ ar(time = t_sec, p = 1, cov = TRUE),
    nl = TRUE
  ),
  data = df1,
  backend = "cmdstanr",
  chains = 1,
  # prior = prior1,
  family = student
)

fit1
conditional_effects(fit1)

df2 = tibble(t_sec = seq(0, max(df1$t_sec, length.out = 1e2)))
p2 = predict(fit1, newdata = df2)
df2 = mutate(df2, gsw = p2[, "Estimate"], lower = p2[, "Q2.5"], upper = p2[, "Q97.5"])

ggplot(df1, aes(t_sec, gsw)) +
  geom_point() +
  geom_line(data = df2, aes(t_sec, gsw), color = "blue") +
  geom_ribbon(data = df2, aes(t_sec, ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Time (s)", y = "gsw")





