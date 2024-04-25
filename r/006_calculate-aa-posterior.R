# Calculate the posterior distribution for each AA estimate
source("r/header.R")

# Summarize RH curves to get overlap of log_gsw values ----
trimmed_rh_curves = read_rds("data/trimmed_rh_curves.rds") |>
  summarize(
    min_log_gsw = min(log_gsw),
    max_log_gsw = max(log_gsw),
    .by = c(
      "acc",
      "acc_id",
      "light_treatment",
      "light_intensity",
      "leaf_type"
    )
  )

  curve_fits_summary = read_rds("objects/curve-fits-summary.rds")

# I need to decide if any of these require re-fitting
curve_fits_summary |>
  filter(rhat > 1.05 | ess_bulk < 100)

# Draws from the curve fits ----
curve_fits_draws = read_rds("objects/curve-fits-draws.rds") |>
  rename(
    b0 = b_Intercept,
    b1 = b_polylog_gsw2rawEQTRUE1,
    b2 = b_polylog_gsw2rawEQTRUE2
  ) |>
  select(starts_with("."), matches("b[0-2]"), file) |>
  mutate(accid_leaftype_lightintensity = str_remove(file, ".rds"), .keep = "unused") |>
  separate_wider_delim(accid_leaftype_lightintensity, delim = "_", 
                       names = c("acc_id", "leaf_type", "light_intensity"))

# Join the range of log_gsw values to the curve fits and calculated AA ----
aa_post = full_join(
  curve_fits_draws,
  trimmed_rh_curves,
  by = join_by(acc_id, leaf_type, light_intensity)
) |>
  mutate(log_gsw = min_log_gsw * (leaf_type == "amphi") + max_log_gsw * (leaf_type == "pseudohypo")) |>
  select(-min_log_gsw, -max_log_gsw) |>
  pivot_wider(names_from = leaf_type, values_from = c(b0, b1, b2, log_gsw)) |>
  mutate(
    upper_int = aa_int(log_gsw_pseudohypo, b0_amphi, b0_pseudohypo, b1_amphi, b1_pseudohypo, b2_amphi, b2_pseudohypo),
    lower_int = aa_int(log_gsw_amphi, b0_amphi, b0_pseudohypo, b1_amphi, b1_pseudohypo, b2_amphi, b2_pseudohypo),
    aa = upper_int - lower_int
  )

# Check that AA estimates converged ----
aa_summary = aa_post |>
  select(starts_with("."), acc_id, aa) |>
  split(~ acc_id) |>
  map_dfr(\(.x) {
    .x |>
      select(-acc_id) |>
      as_draws_df() |>
      posterior::summarise_draws() |>
      mutate(acc_id = .x$acc_id[1])
  })

# I need to decide if any of these require re-fitting
aa_summary |>
  filter(rhat > 1.05 | ess_bulk < 100)

# Write posterior of AA estimates ----
write_rds(aa_post, "objects/aa_post.rds")
