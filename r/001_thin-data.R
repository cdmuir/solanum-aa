# Thin data to reduce oversampling small regions
source("r/header.R")

rh_curves = read_rds("data/rh_curves.rds")

thinned_rh_curves = rh_curves |>
  split( ~ acc_id + curve_type + light_intensity) |>
  map_dfr(thin_data, bin_width = aa_args$thinning_interval, .progress = TRUE)

# stats on thinning
thinning_summary = thinned_rh_curves |>
  summarize(
    n_total = n(),
    n_keep = sum(keep),
    .by = c("acc_id", "light_intensity", "light_treatment", "curve_type")
  ) |>
  mutate(percent_thin = 100 * (1 - n_keep / n_total))

# ggplot(thinning_summary, aes(percent_thin)) +
#   facet_grid(vars(light_intensity, light_treatment), vars(curve_type)) +
#   geom_histogram()

thinned_rh_curves |>
  filter(keep) |>
  select(-keep) |>
  write_rds("data/thinned_rh_curves.rds")

write_rds(thinning_summary, "objects/thinning_summary.rds")
