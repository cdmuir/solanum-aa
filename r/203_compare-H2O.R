# Compare simulated to estimated [H2O] values
source("r/header.R")

dat = list.files("synthetic-data", "df_sim[0-9]{4}.rds", full.names = TRUE)
fit = list.files("objects", "fit_sim[0-9]{4}.rds", full.names = TRUE)
n_dat = str_extract(dat, "[0-9]{4}")
n_fit = str_extract(fit, "[0-9]{4}")

assert_set_equal(n_dat, n_fit)

# w_0 and H2O_r ----
n_dat |>
  map_dfr(\(.x) {
    df_sim = read_rds(glue("synthetic-data/df_sim{.x}.rds"))
    stan_sim = read_rds(glue("synthetic-data/stan_sim{.x}.rds"))
    
    # check order
    assert_true(all(df_sim$curve_number == stan_sim$curve))
    
    fit_sim = read_rds(glue("objects/fit_sim{.x}.rds"))
    
    par_string = "^(w|w_0)\\[([0-9]+),([0-9]+)\\]$"
    
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
    
    # w_0 estimates
    df_w0_est = fit_sim$draws("w_0") |>
      as_draws_df() |>
      pivot_longer(matches(par_string), values_to = "w_0") |>
      mutate(
        row = as.integer(str_replace(name, par_string, "\\2")),
        comp = str_replace(name, par_string, "\\3")) |>
      select(-name) |>
      full_join(select(df_sim, row, curve_number),
                by = join_by(row),
                relationship = "many-to-many")
    
    # Join weights and c_0 estimates
    df_w_w0 = full_join(df_w,
                        df_w0_est,
                        join_by(.chain, .iteration, .draw, curve_number, comp, row)) |>
      pivot_wider(
        id_cols = c(starts_with("."), "row"),
        names_from = "comp",
        values_from = c("w", "w_0")
      ) |>
      mutate(w0_est = (w_1 * w_0_1 + w_2 * w_0_2) / (w_1 + w_2)) |>
      group_by(row) |>
      point_interval(w0_est)
    
    # Join w_0 estimates with simulated values
    df_w0 = full_join(
      df_sim |>
        # Simulated w_0
        select(light_treatment, leaf_type, id, w_0, H2Or_sim = H2O_r) |>
        mutate(row = row_number()),
      
      df_w_w0,
      by = join_by(row)
    )
    
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
    
  }) |>
  write_rds("objects/fit_sim_summary_w0.rds")
