# Summarize parameter convergence and extract draws 
source("r/header.R")

dir_steadystate = "objects/curve-fits/steadystate/"
dir_dynamic = "objects/curve-fits/dynamic/"

steadystate_curve_fits = list.files(dir_steadystate)
dynamic_curve_fits = list.files(dir_dynamic)

# Summary of parameter convergence
steadystate_curve_fits |>
  map_dfr(\(.x) {
    read_rds(paste0(dir_steadystate, .x)) |>
      posterior::summarize_draws() |>
      mutate(file = .x)
  }, .progress = TRUE) |>
  write_rds("objects/sty-summary.rds")

dynamic_curve_fits |>
  map_dfr(\(.x) {
    read_rds(paste0(dir_dynamic, .x)) |>
      posterior::summarize_draws() |>
      mutate(file = .x)
  }, .progress = TRUE) |>
  write_rds("objects/dyn-summary.rds")

# Curve fit draws
steadystate_curve_fits |>
  map_dfr(\(.x) {
    read_rds(paste0(dir_steadystate, .x)) |>
      as_draws_df() |>
      mutate(file = .x)
  }, .progress = TRUE) |>
  write_rds("objects/sty-draws.rds")

dynamic_curve_fits |>
  map_dfr(\(.x) {
    read_rds(paste0(dir_dynamic, .x)) |>
      as_draws_df() |>
      mutate(file = .x)
  }, .progress = TRUE) |>
  write_rds("objects/dyn-draws.rds")
