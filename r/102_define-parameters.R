# Define parameters for simulate synthetic data sets
source("r/header.R")
set.seed(20231227)

aa_hyperpars = read_rds("objects/aa_hyperpars.rds")


# just testing that arima.sim works the way I think
# map_dfr(seq_len(1e4), \(x) {
#   tibble(r = x, name = 1:3, value = arima.sim(list(ar = 0.9), n = 3, sd = 1)[1:3])
# }) |>
#   pivot_wider() |>
#   select(-r) |>
#   cor()

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
      n = nrow(x)
      x |>
        mutate(
          elapsed = seq(0, length.out = n, by = interval),
          R_c = list(make_autocorr_matrix(elapsed, b_autocorr_c)),
          R_w = list(make_autocorr_matrix(elapsed, b_autocorr_w)),
          b_autocorr_c = b_autocorr_c,
          b_autocorr_w = b_autocorr_w,
          error_CO2r = mvnfast::rmvn(1, rep(0, n), R_c[[1]]),
          # error_resid = arima.sim(list(ar = rho_error_resid), n = nrow(x), sd = sigma_error_resid),
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
