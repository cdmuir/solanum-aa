# Summarize parameter convergence and extract draws 
source("r/header.R")

# TEMPORARY: THERE MAY BE A DISCREPANCY BETWEEN THE CURVE FITS AND THE CURVES IN TRIMMED DATA.
# this can be removed once the curve fits are updated to match the trimmed data
rh_curves = read_rds("data/trimmed_rh_curves.rds")

curve_fits = list.files("objects/curve-fits")

all(paste0(unique(rh_curves$curve), ".rds") %in% curve_fits)
x = curve_fits %in% paste0(unique(rh_curves$curve), ".rds")
curve_fits = curve_fits[x]

# Summary of parameter convergence
curve_fits |>
  map_dfr(\(.x) {
    read_rds(paste0("objects/curve-fits/", .x)) |>
      posterior::summarize_draws() |>
      mutate(file = .x)
  }, .progress = TRUE) |>
  write_rds("objects/curve-fits-summary.rds")

# Curve fit draws
curve_fits |>
  map_dfr(\(.x) {
    read_rds(paste0("objects/curve-fits/", .x)) |>
      as_draws_df() |>
      mutate(file = .x)
  }, .progress = TRUE) |>
  write_rds("objects/curve-fits-draws.rds")
