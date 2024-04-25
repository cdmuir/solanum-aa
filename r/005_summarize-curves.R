# Summarize parameter convergence and extract draws 
source("r/header.R")

curve_fits = list.files("objects/curve-fits")

# Summary of parameter convergence
curve_fits |>
  map_dfr(\(.x) {
    read_rds(paste0("objects/curve-fits/", .x)) |>
      posterior::summarize_draws() |>
      mutate(file = .x)
  }) |>
  write_rds("objects/curve-fits-summary.rds")

# Curve fit draws
curve_fits |>
  map_dfr(\(.x) {
    read_rds(paste0("objects/curve-fits/", .x)) |>
      as_draws_df() |>
      mutate(file = .x)
  }) |>
  write_rds("objects/curve-fits-draws.rds")
