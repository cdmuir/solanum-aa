# Compare simulated to estimated [CO2] values
source("r/header.R")

dat = list.files("synthetic-data", "df_sim[0-9]{4}.rds", full.names = TRUE)
fit = list.files("objects", "fit_sim[0-9]{4}.rds", full.names = TRUE)
n_dat = str_extract(dat, "[0-9]{4}")
n_fit = str_extract(fit, "[0-9]{4}")

assert_set_equal(n_dat, n_fit)

# c_0 and CO2_r ----
n_dat |>
  map_dfr(\(.x) {
    df_sim = read_rds(glue("synthetic-data/df_sim{.x}.rds"))
    fit_sim = read_rds(glue("objects/fit_sim{.x}.rds"))
    
    par_string = "^c_0\\[([0-9]+),([0-9]+),([0-9]+),([0-9]+)\\]$"
    
    df_c0 = full_join(
      # Simulated c_0
      df_sim |>
        select(light_treatment, leaf_type, id, pts, c_0, CO2r_sim = CO2_r),
      
      # Estimated c_0
      fit_sim$draws("c_0") |>
        as_draws_df() |>
        pivot_longer(starts_with("c_0"), values_to = "c_0") |>
        mutate(
          pts = str_c(
            "p",
            str_replace(name, par_string, "\\1") |>
              str_pad(2L, "left", "0")
          ),
          id = LETTERS[str_replace(name, par_string, "\\2") |>
                         as.numeric()],
          lt1 = str_replace(name, par_string, "\\3"),
          leaf_type = case_when(lt1 == 1 ~  "amphi",
                                lt1 == 2 ~ "pseudohypo"),
          lt2 = str_replace(name, par_string, "\\4"),
          light_treatment = case_when(lt2 == 1 ~  "high",
                                lt2 == 2 ~ "low")
        ) |>
        select(-name, -lt1, -lt2) |>
        summarize(c0_est = median(c_0), .by = c(pts, id, leaf_type, light_treatment)),
      by = join_by(light_treatment, leaf_type, id, pts)
    )
    
    # Summarize fit
    tibble(
      sim = glue("sim{.x}"),
      var = rep("c_0", 3),
      q1 = c("true", "true", "simulated"),
      q2 = c("simulated", "estimated", "estimated"),
      r = c(
        cor(df_c0$c_0, df_c0$CO2r_sim),
        cor(df_c0$c_0, df_c0$c0_est),
        cor(df_c0$CO2r_sim, df_c0$c0_est)
      )
    )
    
  }) |>
  write_rds("objects/fit_sim_summary_c0.rds")


# Plot
# True c_0 versus simulated CO2_r
# ggplot(df_c0, aes(c_0, CO2r_sim, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   coord_equal()

# Simulated CO2_s versus estimated c_a
# ggplot(df_c0, aes(CO2r_sim, c0_est, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point()

# True c_0 versus estimated c_0
# ggplot(df_c0, aes(c_0, c0_est, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   coord_equal()

# c_a and CO2_s ----
# This is not that informative since c_a is fixed, but code might be useful for
# debugging later
# df_ca = full_join(
#   # Simulated c_a
#   df_sim |>
#     select(leaf_type, rep, pts, c_a, CO2s_sim = CO2_s),
#   
#   # Estimated c_a
#   fit_sim$draws("c_a") |>
#     as_draws_df() |>
#     pivot_longer(starts_with("c_a"), values_to = "c_a") |>
#     mutate(
#       pts = str_c(
#         "p",
#         str_replace(name, par_string, "\\1") |>
#           str_pad(2L, "left", "0")
#       ),
#       rep = str_c(
#         "r",
#         str_replace(name, par_string, "\\2") |>
#           str_pad(2L, "left", "0")
#       ),
#       lt = str_replace(name, par_string, "\\3"),
#       leaf_type = case_when(lt == 1 ~  "amphi",
#                             lt == 2 ~ "pseudohypo")
#     ) |>
#     select(-name, -lt) |>
#     summarize(ca_est = median(c_a), .by = c(pts, rep, leaf_type)),
#   by = join_by(leaf_type, rep, pts)
# )
# 
# # True c_a versus simulated CO2_s
# ggplot(df_ca, aes(c_a, CO2s_sim, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   coord_equal()
# 
# # Simulated CO2_s versus estimated c_a
# ggplot(df_ca, aes(CO2s_sim, ca_est, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point()
# 
# # True c_a versus estimated c_a
# ggplot(df_ca, aes(c_a, ca_est, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   coord_equal()
