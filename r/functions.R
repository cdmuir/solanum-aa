number_to_word = function(num) {
  if (num >= 0 && num <= 10) {
    words = c("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")
    return(words[num + 1])
  } else {
    return(num)
  }
}

## Add information to aa_stats
add_to_stats = function(...) {
  old = read_rds("objects/aa_stats.rds")
  new = list(...)
  if (any(names(new) %in% names(old))) {
    old[intersect(names(old), names(new))] = NULL
  }
  write_rds(c(old, new), "objects/aa_stats.rds")
}

## Function to thin data before analysis
thin_data = function(.d, bin_width, min_n = 20L) {
  
  assert_data_frame(.d)
  assert_number(bin_width, lower = 0)
  assert_int(min_n)
  
  if (nrow(.d) > min_n) {
    d1 = .d |>
      mutate(log_gsw = log(gsw)) |>
      dplyr::select(gsw, log_gsw) |>
      arrange(log_gsw)
    
    n = 1
    
    while (n < min_n) {
      ret = d1 |>
        summarize(min = min(log_gsw) - bin_width,
                  max = max(log_gsw) + bin_width) |>
        reframe(
          bin_start = seq(min, max - bin_width, by = bin_width),
          bin_end = seq(min + bin_width, max, by = bin_width)
        ) |>
        crossing(d1) |>
        filter(log_gsw >= bin_start & log_gsw < bin_end) |>
        mutate(bin_mid = bin_start + bin_width / 2,
               d_mid = abs(log_gsw - bin_mid)) |>
        reframe(
          gsw = gsw,
          log_gsw = log_gsw,
          d_mid = d_mid,
          keep = (d_mid == min(d_mid)),
          .by = "bin_mid"
        )
      
      assert_true(nrow(ret) == nrow(.d))
      ret$keep[1] <- ret$keep[nrow(ret)] <- TRUE
      
      n = sum(ret$keep)
      bin_width = bin_width / 2
    }
    
    return(full_join(.d, dplyr::select(ret, gsw, keep), by = join_by(gsw)))
    
  } else {
    return(mutate(.d, keep = TRUE))
  }
  
}

## Calculate AA from parameter estimates:
# log(A_amphi) = b0_a + b1_a * log(gsw) + b2_a * log(gsw) ^ 2
# log(A_hypo)  = b0_h + b1_h * log(gsw) + b2_h * log(gsw) ^ 2
# Indefinite integral of log(A_amphi / A_hypo) when log(A) is a quadratic 
# function of log(gsw)
aa_int = function(log_gsw, b0_a, b0_h, b1_a, b1_h, b2_a, b2_h) {
  
  log_gsw ^ 3 * (b2_a / 3 - b2_h / 3) + 
    log_gsw ^ 2 * (b1_a / 2 - b1_h / 2) + 
    log_gsw * (b0_a - b0_h)
  
}

## Function to do reverse log10 transformation
# From ChatGPT
reverselog10_trans = function() {
  trans_new(
    name = "reverselog10",
    transform = function(x)
      -log10(x),
    inverse = function(x)
      10^(-x),
    domain = c(1e-100, Inf)
  )
}

## Function to get ggplot panel dimensions
# Help from: https://stackoverflow.com/questions/71250264/figuring-out-panel-size-given-dimensions-for-the-final-plot-in-ggsave-to-show-c

get_panel_dim = function(gp, fig_width, fig_height) {
  output_width = unit(fig_width, "in")
  output_height = unit(fig_height, "in")
  
  # Convert plot to gtable
  gt = ggplotGrob(gp)
  
  # Find panel
  is_panel = grep("panel", gt$layout$name)
  panel_location = gt$layout[is_panel,]
  
  # Get widths/heights of everything except the panel
  width  = gt$widths[-panel_location$l]
  height = gt$heights[-panel_location$t]
  
  # Calculate as output size - size of plot decoration
  panel_height = output_height - sum(height)
  panel_width  = output_width - sum(width)
  
  list(
    panel_height = convertUnit(panel_height, "in"),
    panel_width = convertUnit(panel_width, "in")
  )
}

## Function to refactor treatments for figures
refactor_for_figure = function(.x) {
  
  mutate(
    .x,
    Growth = light_treatment |>
      factor(levels = c("low", "high")) |>
      fct_recode(sun = "high", shade = "low"),
    Measurement = light_intensity |>
      factor(levels = c("150", "2000")) |>
      fct_recode(low = "150", high = "2000")
  )
  
}

# Calculate saturating vapor pressure following the LI6800 manual
li6800_svp = function(T_degreeC, P_kPa) {
  1000 * 0.61365 * exp(17.502 * T_degreeC / (240.97 + T_degreeC)) / P_kPa
}
