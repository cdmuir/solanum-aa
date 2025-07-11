# modeling ideas (move elsewhere, material below is good)
# based on eq. A2 in Buckley et al. (2023)
# X is osmolality gradient (delta pi in buckley paper)
# X_f is target osmolality
# dX / dt = a (X - X_f)  # buckley paper, results in Vico eqn I think
#
# Now suppose that a is a function of X
# a = a0 + a1 * X
# dX / dt = (a0 + a1 * X) * (X - X_f)
# this only has implicit soln according to ChatGpt.
# alternate forms usually result in sigmoidal
# 
# what if we go back to dX / dt = a (X - X_f) but assume that d X_f / dt = b (X_f - X_inf)
# we get this according to chatGPT:
solve_X_dynamics <- function(t, X0, Xf0, X_inf, a, b) {
  if (a == b) {
    # Special case: a == b
    term1 <- X_inf
    term2 <- a * (Xf0 - X_inf) * t * exp(-a * t)
    term3 <- (X0 - X_inf) * exp(a * t)
    return(term1 + term2 + term3)
  } else {
    # General case: a != b
    term1 <- X_inf
    term2 <- (a * (Xf0 - X_inf)) / (b - a) * exp((b - 2 * a) * t)
    term3 <- (X0 - X_inf + (a * (Xf0 - X_inf)) / (b - a)) * exp(a * t)
    return(term1 - term2 + term3)
  }
}

# Fit CDWeibulll mode to each curve separately
source("r/header.R")

sk_dir = "objects/sk-curves/"
if (!dir.exists(sk_dir)) {
  dir.create(sk_dir)
}

rh_curves = read_rds("data/trimmed_rh_curves.rds") |>
  mutate(ci = as.numeric(as.factor(curve))) |>
  mutate(t_sec = elapsed - min(elapsed), .by = ci)
  
form_vico = gsw ~ gf + dg * exp(-(t_sec / tau))
form_cdweibull = gsw ~ gf + dg * exp(-(t_sec / tau) ^ lambda)

bform_vico = bf(
  form_vico,
  gf ~ 1,
  logdg ~ 1,
  tau ~ 1,
  nl = TRUE
)

bform_cdweibull = bf(
  form_cdweibull,
  gf ~ 1,
  logdg ~ 1,
  tau ~ 1,
  lambda ~ 1,
  nl = TRUE
)

# nls ----
safe_nls = safely(nls)

nls_summary = rh_curves |>
  split(~ curve) |>
  future_imap_dfr(\(df, curve_id) {
    # df = filter(rh_curves, curve == "LA0107-E_pseudohypo_2000")
    init = list(gf = min(df$gsw),
                dg = diff(range(df$gsw)),
                tau = 50)
    
    fit_vico = safe_nls(form_vico, data = df, start = init)
    
    if (is.null(fit_vico$result)) {
      b_vico = tibble(
        parameter = names(init),
        value = NA,
        aic = NA,
        model = "vico"
      )
      init = c(init, lambda = 1)
    } else {
      b_vico = as_tibble(coef(fit_vico$result), rownames = "parameter") |>
        mutate(model = "vico", aic = AIC(fit_vico$result))
      init = c(as.list(coef(fit_vico$result)), lambda = 1)
    }
    
    fit_cdweibull = safe_nls(form_cdweibull, data = df, start = init)
    
    if (is.null(fit_cdweibull$result)) {
      b_cdweibull = tibble(
        parameter = names(init),
        value = NA,
        aic = NA,
        model = "cdweibull"
      )
    } else {
      b_cdweibull = as_tibble(coef(fit_cdweibull$result), rownames = "parameter") |>
        mutate(model = "cdweibull", aic = AIC(fit_cdweibull$result))
    }
    
    bind_rows(b_vico, b_cdweibull) |>
      mutate(curve = curve_id)
    
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))


nls_summary |>
  filter(curve == "LA3778-K_pseudohypo_150")

write_rds(nls_summary, "objects/nls_summary.rds")

# brms ----
plan(multisession, workers = 9)
rh_curves |>
  split( ~ curve) |>
  future_iwalk(\(df, curve_id) {
    file = paste0(sk_dir, curve_id, ".rds")
    if (!file.exists(file)) {
      
      # need to figure out somethign with priors/init
      # prior1 = prior(normal(1/3, 1/3), nlpar = "gf") +
      #   prior(normal(-1/4, 1/4), nlpar = "dg", ub = 0) +
      #   # prior(normal(-1, 1), nlpar = "lambda", ub = 0) +
      #   prior(normal(100, 100), nlpar = "tau", lb = 0)
      df = filter(rh_curves, curve == "LA0107-C_pseudohypo_2000")
      init = list(gf = min(df$gsw),
                  dg = diff(range(df$gsw)),
                  tau = 50)
      fit_vico = nls(form_vico, data = df, start = init)
      init = c(as.list(coef(fit_vico)), lambda = 1)
      fit_cdweibull = nls(form_cdweibull, data = df, start = init)
      AIC(fit_vico, fit_cdweibull)
      
      x = 2
      ad = 0.8
      n_divergent = Inf
      
      while (n_divergent > 10) {
        fit_curve = brm(
          formula = fml,
          iter = x * 2000,
          thin = x,
          data = df,
          chains = 1,
          cores = 1,
          # prior = prior1,
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

# zip
zip::zipr(
  "objects/sk-curves.zip",
  list.files(sk_dir, full.names = TRUE),
  recurse = TRUE
)

nls(gsw ~ gf + dg * exp(-(t_sec / tau)),
    data = df,
    start = list(gf = 0.3, dg = 0.2, tau = 60)
)

fit1 = read_rds(paste0(sk_dir, "LA0107-C_pseudohypo_2000.rds")) 
plot(fit1$data$t_sec, fit1$data$gsw)
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





