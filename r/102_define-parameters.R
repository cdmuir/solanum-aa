# Define parameters for simulate synthetic data sets
source("r/header.R")

aa_hyperpars = read_rds("objects/aa_hyperpars.rds")

# testing out alternative way to simulate autocorrelation
# this code generate correlation matrices for MVN()
elapsed = c(seq(0, by = 10, length.out = 4), 100)
m1 = outer(elapsed, elapsed, "-")

R_c = R_w = diag(1, length(elapsed))

ltri_c = exp(-aa_hyperpars$b_autocorr_c * m1[which(lower.tri(m1))])
utri_c = exp(-aa_hyperpars$b_autocorr_c * t(m1)[which(upper.tri(m1))])
ltri_w = exp(-aa_hyperpars$b_autocorr_w * m1[which(lower.tri(m1))])
utri_w = exp(-aa_hyperpars$b_autocorr_w * t(m1)[which(upper.tri(m1))])

R_c[which(lower.tri(R_c))] = ltri_c
R_c[which(upper.tri(R_c))] = utri_c
R_w[which(lower.tri(R_w))] = ltri_w
R_w[which(upper.tri(R_w))] = utri_w

assert_true(isSymmetric(R_c))
assert_true(isSymmetric(R_w))

  t(R_c)[which(upper.tri(R_c))] =
  
R_w[which(lower.tri(R_w))] = R_w[which(upper.tri(R_w))] =
  exp(-aa_hyperpars$b_autocorr_w * m1[which(lower.tri(m1))])

exp(-aa_hyperpars$b_autocorr_c * elapsed)

# just testing that arima.sim works the way I think
arima.sim(list(ar = 0.9), n = 2, sd = 1)

aa_hyperpars = read_rds("objects/aa_hyperpars.rds")

set.seed(20231227)

rep_vector = with(aa_hyperpars, str_c("r", str_pad(
  seq_len(n_rep), floor(log10(n_rep)) + 1, "left", "0"
)))
pts_vector = with(aa_hyperpars, str_c("p", str_pad(
  seq_len(n_pts), floor(log10(n_pts)) + 1, "left", "0"
)))
leaf_type_vector = c("amphi", "pseudohypo")


aa_pars = with(
  aa_hyperpars,
  crossing(nesting(
    crossing(rep = rep_vector,
             leaf_type = leaf_type_vector),
    error_intercept = rnorm(2 * n_rep, 0, sigma_error_intercept)
  ),
  pts = pts_vector) |>
    split(~ rep + leaf_type) |>
    # this will need to be adjusted for time interval between points
    map_dfr(\(x) {
      x |>
        mutate(
          error_resid = arima.sim(list(ar = rho_error_resid), n = nrow(x), sd = sigma_error_resid),
          c_a = c_a,
          flow = flow,
          g_bw = g_bw,
          P = P,
          RH = RH,
          s = s,
          T_air = T_air,
          T_leaf = T_leaf,
          sigma_c = sigma_c,
          sigma_w = sigma_w 
        )
    })
)


df_sim = with(
  aa_hyperpars,
  crossing(
    rep = rep_vector,
    leaf_type = leaf_type_vector,
    nesting(pts = pts_vector,
            i = seq(0, 1, length.out = n_pts))
  ) |>
    full_join(aa_pars, by = join_by(rep, pts, leaf_type)) |>
    mutate(
      min_log_gsw = (leaf_type == "amphi") * log(min_gsw_amphi) +
        (leaf_type == "pseudohypo") * log(min_gsw_pseudohypo),
      max_log_gsw = (leaf_type == "amphi") * log(max_gsw_amphi) +
        (leaf_type == "pseudohypo") * log(max_gsw_pseudohypo),
      log_gsw_real = min_log_gsw + (max_log_gsw - min_log_gsw) * i,
      log_A_real = intercept + error_intercept + slope * log_gsw_real,
      g_sw = exp(log_gsw_real),
      A = exp(log_A_real),
      K = K_amphi * 1 / (leaf_type == "amphi")
    ) 
)

df_sim = df_sim |>
  li6800_simulate() |>
  # NEED TO ADD IN TIME-CORRELATED ERROR USING ARIMA
  li6800_add_error(select(df_sim, sigma_c, sigma_w)) |>
  li6800_estimate() 

ggplot(df_sim, aes(g_sw, A, color = leaf_type)) +
  facet_wrap(~rep) +
  geom_point()

ggplot(df_sim, aes(gsw_hat, A_hat, color = leaf_type)) +
  facet_wrap(~rep) +
  geom_point()

write_rds(df_sim, "synthetic-data/df_sim.rds")
