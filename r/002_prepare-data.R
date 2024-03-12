# Prepare data for analysis
source("r/header.R")

read_rds("data/thinned_rh_curves.rds") |>
  mutate(
    RH = H2O_s / li6800_svp(Tair, Pa),
    leaf_type = case_when(
      curve_type == "1-sided RH" ~ "pseudohypo",
      curve_type == "2-sided RH" ~ "amphi"
    )
  ) |>
  separate_wider_delim(acc_id,
                       "-",
                       names = c("acc", "id"),
                       cols_remove = FALSE) |>
  # make unique ID for each leaf_type within acc_id
  unite("lightintensity_x_acc_id", acc_id, light_intensity, remove = FALSE) |>
  # make unique ID for each curve
  unite("curve", acc_id, leaf_type, light_intensity, remove = FALSE) |>
  mutate(
    log_gsw = log(gsw),
    scaled_log_gsw = (log_gsw - mean(log_gsw)) / sd(log_gsw)
  ) |>
  write_rds("data/prepared_rh_curves.rds")
