# Plot fitted rh curves with AA values ----
source("r/header.R")

rh_curves = read_rds("data/prepared_rh_curves.rds")
df_b = read_rds("objects/df_b.rds")
df_aa = read_rds("objects/df_aa.rds")

sd_log_gsw = sd(rh_curves$log_gsw)
mu_log_gsw = mean(rh_curves$log_gsw)

pdf("figures/rh_curves.pdf", width = 6, height = 4)

plot_list = df_b |>
  filter(str_detect(acc_id, "LA2172")) |>
  split( ~ acc_id + light_intensity + leaf_type) |>
  map_dfr(crossing, i = seq(0, 1, length.out = 10)) |>
  mutate(
    scaled_log_gsw = min_scaled_log_gsw + i * (max_scaled_log_gsw - min_scaled_log_gsw),
    log_A = b0 + b1 * scaled_log_gsw + b2 * scaled_log_gsw ^ 2
  ) |>
  group_by(acc_id, light_intensity, leaf_type, scaled_log_gsw) |>
  point_interval(log_A) |>
  split( ~ acc_id) |>
  magrittr::extract(20) |>
  map(\(.x) {
    a = first(.x$acc_id)
    
    df_aa_acc_id = df_aa |>
      filter(acc_id == a)
    aa_2000 = df_aa_acc_id |>
      filter(light_intensity == 2000) |>
      pull(aa) |>
      signif(3)
    aa_150 = df_aa_acc_id |>
      filter(light_intensity == 150) |>
      pull(aa) |>
      signif(3)
    df_points = filter(rh_curves, acc_id == a) |>
      select(gsw, A, light_intensity, leaf_type)
    
    ggplot(
      mutate(
        .x,
        A = exp(log_A),
        .lower = exp(.lower),
        .upper = exp(.upper),
        gsw = exp(scaled_log_gsw * sd_log_gsw + mu_log_gsw)
      ),
      aes(
        x = gsw,
        y = A,
        fill = leaf_type,
        color = leaf_type,
        linetype = light_intensity
      )
    ) +
      geom_ribbon(alpha = 0.2,
                  mapping = aes(ymin = .lower, ymax = .upper)) +
      geom_line() +
      geom_point(data = df_points, aes(shape = leaf_type), size = 2) +
      labs(title = paste0("Accession: ", a),
           subtitle = bquote(atop(AA[2000] == .(aa_2000),  ~ AA[150] == .(aa_150)))) +
      scale_x_log10() +
      scale_y_log10()
  })

pdf("figures/rh_curves.pdf", width = 6, height = 4)
print(plot_list)
dev.off()
