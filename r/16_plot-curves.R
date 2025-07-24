# Plot quadratic fits to each curve
source("r/header.R")

# set plotting parameters
fig_width = 6
fig_height = 5
fontsize_pt = 10
hjust = -0.1
vjust = 0.5

# load data
accession_info = read_rds("data/accession-info.rds")
plant_info = read_rds("data/plant-info.rds")
aa_summary = read_rds("objects/aa_summary.rds")
fit_curve_diagnostics = read_rds("objects/fit_curve_diagnostics.rds") |>
  mutate(
    acc_id = str_extract(file, "LA[0-9]{4}A*-[A-Z]{1}[AB]{0,1}"),
    `Curve type` = str_extract(file, "amphi|pseudohypo"),
    light_intensity = str_extract(file, "150|2000")
  )

# checks
assert_true(all(!is.na(fit_curve_diagnostics$acc_id))) # Check for NA values in acc_id
assert_true(all(!is.na(fit_curve_diagnostics$`Curve type`))) # Check for NA values in `Curve type`
assert_true(all(!is.na(fit_curve_diagnostics$light_intensity))) # Check for NA values in light_intensity

plan(multisession, workers = 9)
ags_curves = future_map(
  unique(aa_summary$acc_id),
  \(.x) {
    # Extract fit, r2, data, and predictions for each curve type and light intensity
    df1 = filter(aa_summary, acc_id == .x) |>
      crossing(curve_type = c("amphi", "pseudohypo")) |>
      mutate(
        file = glue(
          "objects/curve-fits/{acc_id}_{curve_type}_{light_intensity}.rds"
        ),
        fit = map(file, read_rds),
        r2 = map(fit, bayes_R2),
        data = map(fit, ~ .x$data),
        pred = map2(fit, data, \(.x, .y) {
          d1 = tibble(log_gsw = seq(min(.y$log_gsw), max(.y$log_gsw), length.out = 100))
          predict(.x, newdata = d1) |>
            as_tibble() |>
            bind_cols(d1) |>
            rename(log_A = Estimate)
          
        })
      )
    
    # Data points
    df_points = pmap_dfr(dplyr::select(df1, acc_id, light_intensity, curve_type, data),
                         \(acc_id, light_intensity, curve_type, data) {
                           tibble(
                             acc_id = acc_id,
                             light_intensity = light_intensity,
                             `Curve type` = curve_type,
                             log_gsw = data$log_gsw,
                             log_A = data$log_A
                           )
                         })  |>
      mutate(
        Measurement = light_intensity |>
          factor(levels = c("2000", "150")) |>
          fct_recode(low = "150", high = "2000"),
        across(starts_with("log_"), exp, .names = "{sub('^log_', '', .col)}"),
        gr = paste0(acc_id, "_", `Curve type`, "_", light_intensity)
      )
    
    # Predictions
    df_lines = pmap_dfr(dplyr::select(df1, acc_id, light_intensity, curve_type, pred),
                        \(acc_id, light_intensity, curve_type, pred) {
                          tibble(
                            acc_id = acc_id,
                            light_intensity = light_intensity,
                            `Curve type` = curve_type,
                            log_gsw = pred$log_gsw,
                            log_A = pred$log_A,
                            log_lower = pred$`Q2.5`,
                            log_upper = pred$`Q97.5`
                          )
                        })  |>
      mutate(
        Measurement = light_intensity |>
          factor(levels = c("2000", "150")) |>
          fct_recode(low = "150", high = "2000"),
        across(starts_with("log_"), exp, .names = "{sub('^log_', '', .col)}"),
        gr = paste0(acc_id, "_", `Curve type`, "_", light_intensity)
      )
    
    # Summary statistics
    df_text_r2 = df_lines |>
      arrange(acc_id, `Curve type`, Measurement, gsw) |>
      summarize(
        gsw = last(gsw),
        A = last(A),
        .by = c("acc_id", "Curve type", "light_intensity", "Measurement")
      ) |>
      left_join(fit_curve_diagnostics,
                by = join_by(acc_id, `Curve type`, light_intensity)) |>
      mutate(r2_label = glue("italic(r)^2 == {r2_formatted}", r2_formatted = format(r2, digits = 3))) |>
      mutate(gsw = max(gsw), .by = Measurement)
    
    df_text_aa = df_text_r2 |>
      summarize(
        gsw = max(gsw),
        A = mean(A),
        .by = c("acc_id", "light_intensity", "Measurement")
      ) |>
      left_join(aa_summary, by = join_by(acc_id, light_intensity)) |>
      mutate(aa_label = glue(
        "AA == {aa_formatted} %+-% {se2}",
        aa_formatted = format(median, digits = 3),
        se2 = format(2 * sd, digits = 3)
      ))
    
    # Title
    main = glue(
      "{acc_id} ({sp})",
      acc_id = .x,
      sp = filter(accession_info, accession == str_extract(.x, "LA[0-9]{4}A*"))$species
    )
    sub = glue("growth light intensity: {x}",
               x = filter(plant_info, acc_id == .x)$light_treatment)
    
    # Create the plot
    p = ggplot(df_points, aes(gsw, A, shape = Measurement)) +
      scale_shape_manual(values = c("low" = 19, "high" = 21)) +
      scale_color_manual(values = c(
        "amphi" = "tomato",
        "pseudohypo" = "steelblue"
      )) +
      geom_ribbon(
        data = df_lines,
        mapping = aes(
          ymin = lower,
          ymax = upper,
          group = gr
        ),
        fill = "grey",
        alpha = 0.5
      ) +
      geom_point(
        mapping = aes(color = `Curve type`),
        fill = "white",
        alpha = 0.5
      ) +
      geom_line(data = df_lines,
                mapping = aes(color = `Curve type`, group = gr)) +
      geom_text(
        data = df_text_r2,
        mapping = aes(label = r2_label),
        parse = TRUE,
        hjust = hjust,
        vjust = vjust,
        size = fontsize_pt / 2.845 # convert to mm
      ) +
      geom_text(
        data = df_text_aa,
        mapping = aes(label = aa_label),
        parse = TRUE,
        hjust = hjust,
        vjust = vjust,
        size = fontsize_pt / 2.845 # convert to mm
      ) +
      ggtitle(main, subtitle = sub) +
      xlab(expression(italic(g)[sw] ~ bgroup('[', mol ~ m^-2 ~ s^-1, ']'))) +
      ylab(expression(italic(A) ~ bgroup('[', paste(mu, "mol") ~ m^-2 ~ s^-1, ']')))
    
    panel_dim = get_panel_dim(p, fig_width, fig_height) |>
      map(as.numeric)
    b = ggplot_build(p)
    
    # Adjust y-coordinates to prevent text label overlap
    df_text_r2_adj = full_join(
      dplyr::select(df_text_r2, acc_id, `Curve type`, light_intensity, A, r2_label),
      dplyr::select(df_text_aa, acc_id, light_intensity, A, aa_label),
      by = join_by(acc_id, light_intensity)
    ) |>
      rowwise() |>
      mutate(
        g = list(list(
          r2 = textGrob(parse(text = r2_label), gp = gpar(fontisze = fontsize_pt)),
          aa = textGrob(parse(text = aa_label), gp = gpar(fontisze = fontsize_pt))
        )),
        blah = list(g[["r2"]]),
        r2_label_height_in = convertHeight(grobHeight(g[["r2"]]), "inches", valueOnly = TRUE),
        aa_label_height_in = convertHeight(grobHeight(g[["aa"]]), "inches", valueOnly = TRUE),
        y_range = list(b$layout$panel_params[[1]]$y.range),
        y_units_per_inch = diff(y_range) / panel_dim$panel_height,
        r2_height_data_units = r2_label_height_in * y_units_per_inch,
        aa_height_data_units = aa_label_height_in * y_units_per_inch,
        r2_y_min = A.x - 2 * (1 - vjust) * r2_height_data_units,
        r2_y_max = A.x + 2 * (1 - vjust) * r2_height_data_units,
        aa_y_min = A.y - 2 * (1 - vjust) * aa_height_data_units,
        aa_y_max = A.y + 2 * (1 - vjust) * aa_height_data_units,
        A = case_when(
          (A.x > A.y) &
            (r2_y_min < aa_y_max) ~ A.x + 1.5 * 2 * (1 - vjust) * aa_height_data_units,
          (A.x < A.y) &
            (r2_y_max > aa_y_min) ~ A.x - 1.5 * 2 * (1 - vjust) * aa_height_data_units,
          .default = A.x
        )
      ) |>
      ungroup() |>
      dplyr::select(acc_id, `Curve type`, light_intensity, A) |>
      full_join(dplyr::select(df_text_r2, -A),
                by = join_by(acc_id, `Curve type`, light_intensity))
    
    # Find the layer index for the r2 label and update its data
    r2_layer = p$layers |>
      map_lgl(\(.x) {
        if (is.null(.x$mapping$label)) {
          FALSE
        } else {
          as_name(.x$mapping$label) == "r2_label"
        }
      }) |>
      which()
    
    p$layers[[r2_layer]]$data = df_text_r2_adj
    
    b = ggplot_build(p)
    
    # Adjust x-axis limits to accommodate text labels
    xmax = b$data |>
      map_dfr(\(.x) {
        if ("label" %in% colnames(.x)) {
          .x |>
            mutate(
              g = list(textGrob(
                parse(text = label), gp = gpar(fontisze = fontsize_pt)
              )),
              label_width_in = convertWidth(grobWidth(g[[1]]), "inches", valueOnly = TRUE),
              x_range = list(b$layout$panel_params[[1]]$x.range),
              x_units_per_inch = diff(x_range[[1]]) / panel_dim$panel_width,
              x_width_data_units = label_width_in * x_units_per_inch,
              x_max = x + 2.1 * (1 - hjust) * x_width_data_units
            ) |>
            summarize(x_max = max(x_max))
        } else {
          NULL
        }
      }) |>
      pull(x_max) |>
      max()
    
    return(p + xlim(b$layout$panel_params[[1]]$x.range[1], xmax))
    
  },
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)


pdf("figures/ags-curves.pdf",
    width = fig_width,
    height = fig_height)
for (i in seq_along(ags_curves)) {
  print(ags_curves[[i]])
}

dev.off()
