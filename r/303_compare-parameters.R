# Compare measured to estimated parameter values
source("r/header.R")

fit_aa1 = read_rds("objects/fit_aa1.rds")

library("loo")
options(mc.cores = 10)

log_lik_aa1 = fit_aa1$draws("log_lik")
reff_aa1 = relative_eff(exp(log_lik_aa1)) 
loo_aa1 = loo(log_lik_aa1, r_eff = reff_aa1)
print(loo_aa1)

mcmc_trace(fit_aa1$draws("b_aa_light_treatment_high"))
s = fit_aa1$summary()

fit_aa1$summary(c( "rho_resid", "b0_aa",
"b_aa_light_intensity_2000",
"b_aa_light_treatment_high",
"log_sigma_aa_light_intensity_2000_acc",
"log_sigma_aa_light_intensity_2000_acc_id",
"log_sigma_aa_light_treatment_high_acc",
"log_sigma_aa_light_treatment_high_acc_id",
"log_sigma_aa_acc",
"log_sigma_aa_acc_id",
"b0_log_sigma_aa",
"b_log_sigma_aa_light_intensity_2000",
"b_log_sigma_aa_light_treatment_high"
)) |>
  select(variable, median, q5, q95)

fit_aa1$summary(c( "b_aa_light_treatment_high_acc"
)) |>
  select(variable, median, q5, q95) |>
  print(n = 33)

# NEED TO UPDATE
pars = c("A")

fit_dat = read_rds("objects/fit_dat.rds")

rh_curves = read_rds("data/prepared_rh_curves.rds")
stan_rh_curves = read_rds("data/stan_rh_curves.rds")
    
assert_true(stan_rh_curves$n == nrow(thinned_rh_curves))
# another assertion to make sure order is correct?

par_string = glue("^({.s}|w)\\[([0-9]+),([0-9]+)\\]$", .s = str_c(pars, collapse = "|"))
    
# Component weights per each curve
df_w = fit_dat$draws("w") |>
  as_draws_df() |>
  pivot_longer(matches(par_string), values_to = "w") |>
  mutate(
    curve_number = str_replace(name, par_string, "\\2"),
    comp = str_replace(name, par_string, "\\3")
  ) |>
  select(-name)

fit_dat$draws("lp__")     |> mcmc_trace()
    # |>
      full_join(
        select(df_sim, row, curve_number),
        by = join_by(curve_number),
        relationship = "many-to-many"
      )
    
    # parameter estimates
    df_par_est = fit_sim$draws(pars) |>
      as_draws_df() |>
      pivot_longer(matches(par_string)) |>
      mutate(
        par = str_replace(name, par_string, "\\1"),
        row = as.integer(str_replace(name, par_string, "\\2")),
        comp = str_replace(name, par_string, "\\3")
      ) |>
      select(-name) |>
      pivot_wider(names_from = "par") |>
      full_join(select(df_sim, row, curve_number),
                by = join_by(row),
                relationship = "many-to-many")
    
    # Join weights and parameter estimates
    df_w_par = full_join(df_w,
                        df_par_est,
                        join_by(.chain, .iteration, .draw, curve_number, comp, row)) |>
      pivot_wider(
        id_cols = c(starts_with("."), "row"),
        names_from = "comp",
        values_from = c("w", all_of(pars))
      )
    
    .e = environment()
    
    pars |>
      map(\(.p, n_comp) {
        glue("df_w_par = mutate(df_w_par, {.p}_est = ({numerator}) / ({denominator}))",
             numerator = str_c(glue("{.p}_{i} * w_{i}", i = seq_len(n_comp)), collapse = " + "),
             denominator = str_c(glue("w_{i}", i = seq_len(n_comp)), collapse = " + "))
      }, n_comp = stan_sim$n_comp) |>
      map(\(.x) parse(text = .x)) |>
      walk(eval, envir = .e)
    
    df_w_par = df_w_par |>
      select(starts_with("."), "row", ends_with("_est")) |>
      pivot_longer(cols = ends_with("_est")) |>
      group_by(name, row) |>
      point_interval(value)
    
    # Join parameter estimates with simulated values and summarize fit
    pars |>
      map_dfr(\(.p){
        .d = full_join(
          df_w_par |>
          filter(name == paste0(.p, "_est")),
        df_sim |>
          select(row, light_treatment, leaf_type, id, all_of(c(.p, hats[[.p]]))),
        by = join_by(row)
        )
        truth = .d[.p] # true value
        hat = .d[hats[.p]] # LICOR estimate (simulated)
        est = .d$value # stan estimate
        
        tibble(
          sim = glue("sim{.x}"),
          var = rep(.p, 3),
          q1 = c("true", "true", "simulated"),
          q2 = c("simulated", "estimated", "estimated"),
          r = c(
            cor(truth, hat),
            cor(truth, est),
            cor(hat, est)
          )
        )
      })
    
  }) |>
  write_rds("objects/fit_sim_summary_pars.rds")
