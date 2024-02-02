# Prepare synthetic data for Stan

source("r/header.R")

df_sim = read_rds("synthetic-data/df_sim.rds")

stan_sim = list(
  n_pts = length(unique(df_sim$pts)),
  n_rep = length(unique(df_sim$rep)),
  n_leaf_type = length(unique(df_sim$leaf_type))
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
      "stan_sim${var} = array(df_sim${var}, dim = c({i}, {j}, {k}))",
      i = stan_sim$n_pts,
      j = stan_sim$n_rep,
      k = stan_sim$n_leaf_type
    )
  }) |>
  map(\(s) parse(text = s)) |>
  walk(eval, envir = globalenv())

write_rds(stan_sim, "objects/stan_sim.rds")