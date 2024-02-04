# Define parameters for simulate synthetic data sets
source("r/header.R")
set.seed(20231227)

aa_hyperpars = read_rds("objects/aa_hyperpars.rds")

aa_hyperpars |>
  transpose() |>
  iwalk(\(.hpar, .i) {
    
    id_vector = LETTERS[seq_len(.hpar$n_id)]
    pts_vector = with(.hpar, str_c("p", str_pad(
      seq_len(n_pts), floor(log10(n_pts)) + 1, "left", "0"
    )))
    leaf_type_vector = c("amphi", "pseudohypo")
    
    aa_pars = 
      with(.hpar,
           crossing(
             nesting(
               crossing(
                 
                 # id level coefficients (should change to acc_id levels when
                 # acc is added)
                 tibble(
                   id = id_vector,
                   b_intercept_id = rnorm(n_id, 0, sigma_intercept_id),
                   b_slope_id = rnorm(n_id, 0, sigma_slope_id)
                 ),
                 leaf_type = leaf_type_vector
                 
                 ),
               
               # leaf_type:id coefficients (should change to leaf_type:accc_id
               # when acc is added)
               b_intercept_error = rnorm(2 * n_id, 0, sigma_intercept_error)
             ),
             pts = pts_vector
           ) |>
             split( ~ id + leaf_type) |>
             map_dfr(\(x) {
               n = nrow(x)
               x |>
                 mutate(
                   elapsed = seq(0, length.out = n, by = interval),
                   R_c = list(make_autocorr_matrix(elapsed, b_autocorr_c)),
                   R_w = list(make_autocorr_matrix(elapsed, b_autocorr_w)),
                   b_autocorr_c = b_autocorr_c,
                   b_autocorr_w = b_autocorr_w,
                   c_a = c_a,
                   flow = flow,
                   g_bw = g_bw,
                   P = P,
                   RH = RH,
                   s = s,
                   T_air = T_air,
                   T_leaf = T_leaf,
                   sigma_c = sigma_c,
                   sigma_w = sigma_w,
                   error_CO2r = sigma_c * rmvn(1, rep(0, n), R_c[[1]])[1, ],
                   error_CO2s = sigma_c * rmvn(1, rep(0, n), R_c[[1]])[1, ],
                   error_H2Or = sigma_w * rmvn(1, rep(0, n), R_w[[1]])[1, ],
                   error_H2Os = sigma_w * rmvn(1, rep(0, n), R_w[[1]])[1, ]
                 )
               })
             )
    
    df_sim = with(
      .hpar,
      crossing(
        leaf_type = leaf_type_vector,
        id = id_vector,
        nesting(pts = pts_vector,
                i = seq(0, 1, length.out = n_pts))
      ) |>
        full_join(aa_pars, by = join_by(id, pts, leaf_type)) |>
        mutate(
          min_log_gsw = (leaf_type == "amphi") * log(min_gsw_amphi) +
            (leaf_type == "pseudohypo") * log(min_gsw_pseudohypo),
          max_log_gsw = (leaf_type == "amphi") * log(max_gsw_amphi) +
            (leaf_type == "pseudohypo") * log(max_gsw_pseudohypo),
          intercept = mu_intercept + b_intercept_error + b_intercept_id,
          slope = mu_slope + b_slope_id,
          log_gsw = min_log_gsw + (max_log_gsw - min_log_gsw) * i,
          A = intercept + slope * log_gsw,
          g_sw = exp(log_gsw),
          K = K_amphi * 1 / (leaf_type == "amphi")
        )
    )
    
    df_sim = df_sim |>
      li6800_simulate() |>
      # add autocorrelated error
      mutate(
        CO2_r = c_0 + error_CO2r,
        CO2_s = c_a + error_CO2s,
        H2O_r = w_0 + error_H2Or,
        H2O_s = w_a + error_H2Os
      ) |>
      li6800_estimate()
    
    # Quick plot for checking
    # ggplot(df_sim, aes(g_sw, A)) +
    #   facet_wrap(~ id) +
    #   geom_point()
    
    write_rds(df_sim, glue("synthetic-data/df_sim{n}.rds", n = str_pad(.i, 4L, "left", "0")))
    
  })
