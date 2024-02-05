# Compare simulated to estimated A values
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
    
    par_string = "^A\\[([0-9]+),([0-9]+),([0-9]+),([0-9]+)\\]$"
    df_A = full_join(
      # Simulated A
      df_sim |>
        select(light_treatment, leaf_type, id, pts, A, A_sim = A_hat),
      
      # Estimated A
      fit_sim$draws("A") |>
        as_draws_df() |>
        pivot_longer(starts_with("A"), values_to = "A") |>
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
        summarize(A_est = median(A), .by = c(pts, id, leaf_type, light_treatment)),
      by = join_by(light_treatment, leaf_type, id, pts)
    )
    
    # Summarize fit
    tibble(
      sim = glue("sim{.x}"),
      var = rep("g_sw", 3),
      q1 = c("true", "true", "simulated"),
      q2 = c("simulated", "estimated", "estimated"),
      r = c(
        cor(df_A$A, df_A$A_sim),
        cor(df_A$A, df_A$A_est),
        cor(df_A$A_sim, df_A$A_est)
      )
    )
    
  }) |>
  write_rds("objects/fit_sim_summary_A.rds")

# # Plot
# # True A versus simulated A
# ggplot(df_A, aes(A, A_sim, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   coord_equal()
# 
# # Simulated A versus estimated A
# ggplot(df_A, aes(A_sim, A_est, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   coord_equal()
# 
# # True A versus estimated A
# ggplot(df_A, aes(A, A_est, color = leaf_type)) +
#   facet_wrap(~ rep) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   coord_equal()
