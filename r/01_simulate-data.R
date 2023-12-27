# Simulate synthetic data sets to validate statistical models
source("r/header.R")

set.seed(20231227)

# hyper parameters
aa_hyperpars = list(
  n_acc = 1, # number of accessions
  n_rep = 1e1, # number of replicates per accession per treatment
  n_pts = 1e1, # number of points per curve
  min_gsw_amphi = 0.10,
  max_gsw_amphi = 0.50,
  min_gsw_pseudohypo = 0.05,
  max_gsw_pseudohypo = 0.25,
  intercept = 0, # later on, change to mu_intercept, sigma_intercept
  slope = 1,  # later on, change to mu_slope, sigma_slope
  sigma_error_intercept = 1,
  sigma_error_resid = 0.1
)

rep_vector = with(aa_hyperpars, str_c("r", str_pad(
  seq_len(n_rep), floor(log10(n_rep)) + 1, "left", "0"
)))
pts_vector = with(aa_hyperpars, str_c("p", str_pad(
  seq_len(n_pts), floor(log10(n_pts)) + 1, "left", "0"
)))
leaf_type_vector = c("amphi", "pseudohypo")

aa_pars = with(aa_hyperpars,
               crossing(nesting(
                 crossing(rep = rep_vector,
                          leaf_type = leaf_type_vector),
                 error_intercept = rnorm(2 * n_rep, 0, sigma_error_intercept)
               ),
               pts = pts_vector) |>
                 mutate(error_resid = rnorm(2 * n_rep * n_pts, 0,  sigma_error_resid)))

df_sim = with(
  aa_hyperpars,
  crossing(
    rep = rep_vector,
    leaf_type = leaf_type_vector,
    nesting(
      pts = pts_vector,
      i = seq(0, 1, length.out = n_pts)
    )
  ) |>
    full_join(aa_pars, by = join_by(rep, pts, leaf_type)) |>
    mutate(
      min_log_gsw = (leaf_type == "amphi") * log(min_gsw_amphi) +
        (leaf_type == "pseudohypo") * log(min_gsw_pseudohypo),
      max_log_gsw = (leaf_type == "amphi") * log(max_gsw_amphi) +
        (leaf_type == "pseudohypo") * log(max_gsw_pseudohypo),
      log_gsw = min_log_gsw + (max_log_gsw - min_log_gsw) * i,
      log_A = intercept + error_intercept + slope * log_gsw + error_resid
    )
)

# ggplot(df_sim, aes(log_gsw, log_A, color = leaf_type)) +
#   facet_wrap(~rep) +
#   geom_point()

write_rds(df_sim, "synthetic-data/df_sim.rds")
