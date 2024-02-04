# Compare simulated to estimated g_sw values
source("r/header.R")

dat = list.files("synthetic-data", "df_sim[0-9]{4}.rds", full.names = TRUE)
fit = list.files("objects", "fit_sim[0-9]{4}.rds", full.names = TRUE)
n_dat = str_extract(dat, "[0-9]{4}")
n_fit = str_extract(fit, "[0-9]{4}")

assert_set_equal(n_dat, n_fit)

n_dat |>
  map_dfr(\(.x) {
    df_sim = read_rds(glue("synthetic-data/df_sim{.x}.rds"))
    fit_sim = read_rds(glue("objects/fit_sim{.x}.rds"))
    
    par_string = "^g_sw\\[([0-9]+),([0-9]+),([0-9]+)\\]$"
    
    df_gsw = full_join(
      # Simulated g_sw
      df_sim |>
        select(leaf_type, id, pts, g_sw, gsw_sim = gsw_hat),
      
      # Estimated g_sw
      fit_sim$draws("g_sw") |>
        as_draws_df() |>
        pivot_longer(starts_with("g_sw"), values_to = "g_sw") |>
        mutate(
          pts = str_c(
            "p",
            str_replace(name, par_string, "\\1") |>
              str_pad(2L, "left", "0")
          ),
          id = LETTERS[str_replace(name, par_string, "\\2") |>
                         as.numeric()],
          lt = str_replace(name, par_string, "\\3"),
          leaf_type = case_when(lt == 1 ~  "amphi",
                                lt == 2 ~ "pseudohypo")
        ) |>
        select(-name,-lt) |>
        summarize(
          gsw_est = median(g_sw),
          .by = c(pts, id, leaf_type)
        ),
      by = join_by(leaf_type, id, pts)
    )
    
    # Summarize fit
    tibble(
      sim = glue("sim{.x}"),
      var = rep("g_sw", 3),
      q1 = c("true", "true", "simulated"),
      q2 = c("simulated", "estimated", "estimated"),
      r = c(
        cor(df_gsw$g_sw, df_gsw$gsw_sim),
        cor(df_gsw$g_sw, df_gsw$gsw_est),
        cor(df_gsw$gsw_sim, df_gsw$gsw_est)
      )
    )
    
  }) |>
  write_rds("objects/fit_sim_summary_gsw.rds")

# Plot
# # True g_sw versus simulated g_sw
# ggplot(df_gsw, aes(g_sw, gsw_sim, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   coord_equal()
# 
# # Simulated g_sw versus estimated g_sw
# ggplot(df_gsw, aes(gsw_sim, gsw_est, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   coord_equal()
# 
# # True g_sw versus estimated g_sw
# ggplot(df_gsw, aes(g_sw, gsw_est, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   coord_equal()
