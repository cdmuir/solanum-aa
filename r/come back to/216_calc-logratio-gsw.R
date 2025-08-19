source("r/header.R")

curve_fits = list.files("objects/curve-fits", full.names = TRUE)

curve_fits_draws = read_rds("objects/curve-fits-draws.rds") |>
  rename(b0 = b_Intercept, b1 = b_polylog_gsw2rawEQTRUE1, b2 = b_polylog_gsw2rawEQTRUE2) |>
  dplyr::select(starts_with("."), matches("b[0-2]"), file) |>
  mutate(accid_leaftype_lightintensity = str_remove(file, ".rds")) |>
  separate_wider_delim(
    accid_leaftype_lightintensity,
    delim = "_",
    names = c("acc_id", "leaf_type", "light_intensity")
  )

tmp = curve_fits_draws |>
  summarize(.by = c(acc_id, light_intensity)) |>
  reframe(map2_dfr(acc_id, light_intensity, \(id, li) {
    
    # Find Amax of pseudohypo curve
    m1 = read_rds(curve_fits[str_detect(curve_fits, glue("{id}_pseudohypo_{li}"))])
    max_log_A = max(m1$data$log_A)
    approx_log_gsw = m1$data$log_gsw[which.max(m1$data$log_A)]
    
    # Find log_gsw in pseudohypo curve need to achieve max_log_A
    post_hypo = curve_fits_draws |>
      filter(acc_id == id, leaf_type == "pseudohypo", light_intensity == li) |>
      mutate(
        x1 = (-b1 + sqrt(b1^ 2 - 4 * b2 * (b0 - max_log_A))) / (2 * b2),
        x2 = (-b1 - sqrt(b1^ 2 - 4 * b2 * (b0 - max_log_A))) / (2 * b2),
        dx1 = (x1 - approx_log_gsw) ^ 2,
        dx2 = (x2 - approx_log_gsw) ^ 2,
        log_gsw = ifelse(dx1 < dx2, x1, x2)
      ) |>
      pull(log_gsw)
    
    # Find log_gsw in amphi curve need to achieve max_log_A
    m2 = read_rds(curve_fits[str_detect(curve_fits, glue("{id}_amphi_{li}"))])
    approx_log_gsw = m2$data |>
      mutate(dA = (log_A - max_log_A) ^ 2) |>
      filter(dA == min(dA)) |>
      pull(log_gsw)
    
    post_amphi = curve_fits_draws |>
      filter(acc_id == id, leaf_type == "amphi", light_intensity == li) |>
      mutate(
        x1 = (-b1 + sqrt(b1^ 2 - 4 * b2 * (b0 - max_log_A))) / (2 * b2),
        x2 = (-b1 - sqrt(b1^ 2 - 4 * b2 * (b0 - max_log_A))) / (2 * b2),
        dx1 = (x1 - approx_log_gsw) ^ 2,
        dx2 = (x2 - approx_log_gsw) ^ 2,
        log_gsw = ifelse(dx1 < dx2, x1, x2)
      ) |>
      pull(log_gsw)
    
    tibble(
      acc_id = id,
      light_intensity = li,
      gsw_log_ratio = (post_amphi - post_hypo)
    )
  }, .progress = TRUE))


tmp |>
  filter(is.na(gsw_log_ratio), light_intensity == "2000") |>
  summarize(n = n(), .by = c(acc_id, light_intensity)) |>
  arrange(n) |>
  print(n = Inf)

tmp |>
  filter(!is.na(gsw_log_ratio)) |>
  summarize(n = n(),
            gsw_log_ratio = median(gsw_log_ratio),
            .by = c(acc_id, light_intensity)) |>
  filter(n > 3500) |>
  mutate(percent_change = (exp(gsw_log_ratio) - 1) * 100) |>
  pull(percent_change) |> hist()


# Specify plant
id = "LA3124-O"

# Specify measurement light intensity
li = "2000"

# Find Amax of pseudohypo curve
m1 = read_rds(curve_fits[str_detect(curve_fits, glue("{id}_pseudohypo_{li}"))])
max_log_A = max(m1$data$log_A)
approx_log_gsw = m1$data$log_gsw[which.max(m1$data$log_A)]

# Find log_gsw in pseudohypo curve need to achieve max_log_A
post_hypo = curve_fits_draws |>
  filter(acc_id == id, leaf_type == "pseudohypo", light_intensity == li) |>
  mutate(
    x1 = (-b1 + sqrt(b1^ 2 - 4 * b2 * (b0 - max_log_A))) / (2 * b2),
    x2 = (-b1 - sqrt(b1^ 2 - 4 * b2 * (b0 - max_log_A))) / (2 * b2),
    dx1 = (x1 - approx_log_gsw) ^ 2,
    dx2 = (x2 - approx_log_gsw) ^ 2,
    log_gsw = ifelse(dx1 < dx2, x1, x2)
  ) |>
  pull(log_gsw)

# Find log_gsw in amphi curve need to achieve max_log_A
m2 = read_rds(curve_fits[str_detect(curve_fits, glue("{id}_amphi_{li}"))])
approx_log_gsw = m2$data |>
  mutate(dA = (log_A - max_log_A) ^ 2) |>
  filter(dA == min(dA)) |>
  pull(log_gsw)

post_amphi = curve_fits_draws |>
  filter(acc_id == id, leaf_type == "amphi", light_intensity == li) |>
  mutate(
    x1 = (-b1 + sqrt(b1^ 2 - 4 * b2 * (b0 - max_log_A))) / (2 * b2),
    x2 = (-b1 - sqrt(b1^ 2 - 4 * b2 * (b0 - max_log_A))) / (2 * b2),
    dx1 = (x1 - approx_log_gsw) ^ 2,
    dx2 = (x2 - approx_log_gsw) ^ 2,
    log_gsw = ifelse(dx1 < dx2, x1, x2)
  ) |>
  pull(log_gsw)


(exp(post_amphi - post_hypo) - 1) * 100

tibble(
  acc_id = id,
  light_intensity = li,
  gsw_log_ratio = (post_amphi - post_hypo),
))