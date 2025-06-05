source("r/header.R")

stan_rh_curves = read_rds("data/stan_rh_curves.rds")
stan_rh_curves$acc1 = NULL

m5 = cmdstan_model("stan/solanum-aa5.stan", dir = "stan/bin")

fit5 = m5$sample(
  data = stan_rh_curves,
  seed = 508504744,
  chains = 4,
  parallel_chains = 4,
  refresh = 2e2,
  iter_warmup = 2e3,
  iter_sampling = 2e3,
  thin = 2
)

# no convergence at all even with 2e3 iterations.
fit5$save_object("objects/solanum-aa5.rds")
fit5 = read_rds("objects/solanum-aa5.rds")
s = fit5$summary()
mcmc_trace(fit5$draws("lp__"))
mcmc_trace(fit5$draws("b_aa_light_intensity_2000"))
s$variable[5000:6000]

test = cmdstan_model("stan/test.stan", dir = "stan/bin")

# FOR LA0107, 2e3 sampling iterations gets convergence for most things. Only the B_curve() parameters need some time.
# FOR LA0107 + LA4778, 4e3 sampling iterations converges well. Scaling up to all species before trying resid
# FOR ALL ACCESSIONS, 4e3 sampling iterations converges well
# try:
# adding resid covariance back in.

fit = test$sample(
  data = stan_rh_curves,
  seed = 508504744,
  chains = 4,
  parallel_chains = 4,
  refresh = 4e2,
  iter_warmup = 4e3,
  iter_sampling = 4e3,
  thin = 4
)

s = fit$summary()
s |>
  arrange(desc(rhat)) |>
  print(n = 100)

i = 900
ii = (stan_rh_curves$curve == i)

# data
df1 = tibble(
  scaled_log_gsw = stan_rh_curves$scaled_log_gsw[ii],
  A = stan_rh_curves$A[ii],
) |>
  mutate(log_A = log(A))

# predictions
df2 = fit$draws(c("Mu_curve", glue("B_curve[{i},{x}]", x= 1:3))) |>
  as_draws_df() |>
  dplyr::rename(B_curve1 = !!glue("B_curve[{i},1]"),
                B_curve2 = !!glue("B_curve[{i},2]"),
                B_curve3 = !!glue("B_curve[{i},3]")) |>
  mutate(b0 = `Mu_curve[1]` + B_curve1,
         b1 = `Mu_curve[2]` + B_curve2,
         b2 = `Mu_curve[3]` + B_curve3, .keep = "unused") |>
  crossing(scaled_log_gsw = seq(min(stan_rh_curves$scaled_log_gsw[ii]),
                                 max(stan_rh_curves$scaled_log_gsw[ii]),
                                 length.out = 100)) |>
  mutate(log_A = b0 + b1 * scaled_log_gsw + b2 * scaled_log_gsw^2) |>
  group_by(scaled_log_gsw) |>
  point_interval(log_A)

ggplot(df2, aes(scaled_log_gsw, log_A)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  geom_line() +
  geom_point(data = df1, aes(scaled_log_gsw, log_A), color = "red") +
  labs(x = "Scaled log(gsw)", y = "log(A)")
