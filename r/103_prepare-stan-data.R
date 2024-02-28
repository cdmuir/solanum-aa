# Prepare synthetic data for Stan

source("r/header.R")

list.files("synthetic-data", pattern = "df_sim[0-9]{4}.rds", full.names = TRUE) |>
  walk(\(.x) {
    n = str_extract(.x, "[0-9]{4}")
    df_sim = read_rds(.x)
    
    stan_sim = df_sim |>
      # select variables needed for Stan
      dplyr::select(
        A,
        curve,
        id,
        leaf_type,
        leaftype_x_id,
        light_intensity,
        light_treatment,
        elapsed,
        flow,
        g_bw,
        g_sw,
        K,
        P,
        RH,
        s,
        T_air,
        T_leaf,
        CO2_r,
        CO2_s,
        H2O_r,
        H2O_s
      ) |>
      compose_data()
    
    # Manual changes
    stan_sim$n_comp = 2
    stan_sim$n_pts = df_sim |>
      summarise(n = n(), .by = "curve") |>
      pull(n)

    write_rds(stan_sim, glue("synthetic-data/stan_sim{n}.rds"))
    
  })
