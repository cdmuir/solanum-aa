# Prepare synthetic data for Stan

source("r/header.R")

list.files("synthetic-data", pattern = "df_sim[0-9]{4}.rds", full.names = TRUE) |>
  walk(\(.x) {
    n = str_extract(.x, "[0-9]{4}")
    df_sim = read_rds(.x)
    df_sim_env = environment()
    
    stan_sim = list(
      n_pts = length(unique(df_sim$pts)),
      n_id = length(unique(df_sim$id)),
      n_leaf_type = length(unique(df_sim$leaf_type)),
      n_light_treatment = length(unique(df_sim$light_treatment))
    )
    
    c(
      "elapsed",
      "flow",
      "g_bw",
      "g_sw",
      "K",
      "P",
      "RH",
      "s",
      "T_air",
      "T_leaf",
      "CO2_r",
      "CO2_s",
      "H2O_r",
      "H2O_s"
    ) |>
      map(\(var) {
        glue(
          "stan_sim${var} = array(df_sim${var}, dim = c({i}, {j}, {k}, {l}))",
          i = stan_sim$n_pts,
          j = stan_sim$n_id,
          k = stan_sim$n_leaf_type,
          l = stan_sim$n_light_treatment
        )
      }) |>
      map(\(s) parse(text = s)) |>
      walk(eval, envir = df_sim_env)
    
    write_rds(stan_sim, glue("synthetic-data/stan_sim{n}.rds"))
    
  })
