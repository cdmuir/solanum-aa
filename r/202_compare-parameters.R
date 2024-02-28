# Compare simulated to estimated parameter values
source("r/header.R")

dat = list.files("synthetic-data", "df_sim[0-9]{4}.rds", full.names = TRUE)
fit = list.files("objects", "fit_sim[0-9]{4}.rds", full.names = TRUE)
n_dat = str_extract(dat, "[0-9]{4}")
n_fit = str_extract(fit, "[0-9]{4}")

assert_set_equal(n_dat, n_fit)

pars = c("A", "c_0", "g_sw", "w_0")

# for each parameter, there should be a corresponding "hat" value, which is the 
# estimate that the LI-6800 would report. These are called different things in 
# the synthetic data.
hats = c(A = "A_hat", c_0 = "CO2_r", g_sw = "gsw_hat", w_0 = "H2O_r")
assert_set_equal(pars, names(hats))

.x = "0001"
# n_dat |>
  # map_dfr(\(.x) {
    df_sim = read_rds(glue("synthetic-data/df_sim{.x}.rds"))
    stan_sim = read_rds(glue("synthetic-data/stan_sim{.x}.rds"))
    
    # check order
    assert_true(all(df_sim$curve_number == stan_sim$curve))
    
    fit_sim = read_rds(glue("objects/fit_sim{.x}.rds"))
    
    par_string = glue("^({.s}|w)\\[([0-9]+),([0-9]+)\\]$", .s = str_c(pars, collapse = "|"))
    
    # Component weights per each curve
    df_w = fit_sim$draws("w") |>
      as_draws_df() |>
      pivot_longer(matches(par_string), values_to = "w") |>
      mutate(
        curve_number = str_replace(name, par_string, "\\2"),
        comp = str_replace(name, par_string, "\\3")
      ) |>
      select(-name) |>
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
    
    # Join parameter estimates with simulated values
    df_par = pars |>
      map(\(.p) {
        full_join(
          df_w_par |>
           filter(name == paste0(.p, "_est")),
          df_sim |>
            select(row, light_treatment, leaf_type, id, all_of(c(.p, hats[[.p]]))),
          by = join_by(row)
        )
      })

    # WORKING HERE
    # Summarize fit
    tibble(
      sim = glue("sim{.x}"),
      var = rep("w_0", 3),
      q1 = c("true", "true", "simulated"),
      q2 = c("simulated", "estimated", "estimated"),
      r = c(
        cor(df_w0$w_0, df_w0$H2Or_sim),
        cor(df_w0$w_0, df_w0$w0_est),
        cor(df_w0$H2Or_sim, df_w0$w0_est)
      )
    )
    
  # }) |>
  # write_rds("objects/fit_sim_summary_pars.rds")
