# Compare simulated to estimated [CO2] values
source("r/header.R")

dat = list.files("synthetic-data", "df_sim[0-9]{4}.rds", full.names = TRUE)
fit = list.files("objects", "fit_sim[0-9]{4}.rds", full.names = TRUE)
n_dat = str_extract(dat, "[0-9]{4}")
n_fit = str_extract(fit, "[0-9]{4}")

assert_set_equal(n_dat, n_fit)

n_dat |>
  map_dfr(\(.x) {
    df_sim = read_rds(glue("synthetic-data/df_sim{.x}.rds"))
    stan_sim = read_rds(glue("synthetic-data/stan_sim{.x}.rds"))
    
    # check order
    assert_true(all(df_sim$curve_number == stan_sim$curve))
    
    fit_sim = read_rds(glue("objects/fit_sim{.x}.rds"))
    
    par_string = "^(w|c_0)\\[([0-9]+),([0-9]+)\\]$"
    
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
    
    # c_0 estimates
    df_c0_est = fit_sim$draws("c_0") |>
      as_draws_df() |>
      pivot_longer(matches(par_string), values_to = "c_0") |>
      mutate(# parameter = str_replace(name, par_string, "\\1"),
        row = as.integer(str_replace(name, par_string, "\\2")),
        comp = str_replace(name, par_string, "\\3")) |>
      select(-name) |>
      full_join(select(df_sim, row, curve_number),
                by = join_by(row),
                relationship = "many-to-many")
    
    # Join weights and c_0 estimates
    df_w_c0 = full_join(df_w,
                        df_c0_est,
                        join_by(.chain, .iteration, .draw, curve_number, comp, row)) |>
      pivot_wider(
        id_cols = c(starts_with("."), "row"),
        names_from = "comp",
        values_from = c("w", "c_0")
      ) |>
      mutate(c0_est = (w_1 * c_0_1 + w_2 * c_0_2) / (w_1 + w_2)) |>
      group_by(row) |>
      point_interval(c0_est)
    
    # Join c_0 estimates with simulated values
    df_c0 = full_join(
      df_sim |>
        # Simulated c_0
        select(light_treatment, leaf_type, id, c_0, CO2r_sim = CO2_r) |>
        mutate(row = row_number()),
      
      df_w_c0,
      by = join_by(row)
    )
    
    # Summarize fit
    tibble(
      sim = glue("sim{.x}"),
      var = rep("c_0", 3),
      q1 = c("true", "true", "simulated"),
      q2 = c("simulated", "estimated", "estimated"),
      r = c(
        cor(df_c0$c_0, df_c0$CO2r_sim),
        cor(df_c0$c_0, df_c0$c0_est),
        cor(df_c0$CO2r_sim, df_c0$c0_est)
      )
    )
    
  }) |>
  write_rds("objects/fit_sim_summary_c0.rds")
