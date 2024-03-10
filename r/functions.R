# 0. Import data ----

## Function to thin data before analysis
thin_data = function(.d, bin_width, min_n = 20L) {
  
  assert_data_frame(.d)
  assert_number(bin_width, lower = 0)
  assert_int(min_n)
  
  if (nrow(.d) > min_n) {
    d1 = .d |>
      mutate(log_gsw = log(gsw)) |>
      select(gsw, log_gsw) |>
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
    
    return(full_join(.d, select(ret, gsw, keep), by = join_by(gsw)))
    
  } else {
    return(mutate(.d, keep = TRUE))
  }
  
}

## Set of functions to identify linear portions of each A-gsw curve
get_2dsmooth = function(.fit) {
  
  if (inherits(.fit, "loess")) log_gsw = .fit$x
  if (inherits(.fit, "gam")) log_gsw = .fit$model$`log(gsw)`
  Ahat = predict(.fit)
  
  #calculate the first deriative and the new mean x value
  Aprime <- Ahat[-1] - diff(Ahat) / 2
  dAdg <- diff(Ahat) / diff(log_gsw)
  
  #calculate the 2nd deriative and the new mean x value
  Apprime <- Aprime[-1] - diff(Aprime)/2
  d2Adg2 <- diff(dAdg) / diff(Aprime)
  
  if (inherits(.fit, "loess")) return(tibble(d2Adg2 = d2Adg2[,1]))
  if (inherits(.fit, "gam")) return(tibble(d2Adg2 = d2Adg2))
  
}

simulate_smooth = function(.fit, n_sim) {
  
  assert_class(.fit, "lm")
  assert_int(n_sim)
  
  crossing(gsw = exp(.fit$model$`log(gsw)`),
           sim = paste0("sim", str_pad(
             seq_len(n_sim), ceiling(log10(n_sim + 1)), "left", "0"
           ))) %>%
    mutate(A = predict(.fit, .) + rnorm(nrow(.), 0, sigma(.fit))) |>
    split( ~ sim) %>%
    # map(loess, formula = A ~ log(gsw)) |>
    map(\(.x) {
      gam(A ~ s(log(gsw)), data = .x)
    })
  
}

get_cutoff = function(syn, probs) {
  assert_list(syn)
  assert_number(probs, lower = 0, upper = 1)
  syn |>
    map_dfr(get_2dsmooth) |>
    pull(d2Adg2) |>
    abs() |>
    quantile(probs = probs)
}

find_longest_linear = function(.d, .fit, cutoff, min_gap = 3L) {
  
  d1 = .d |>
    arrange(gsw) |>
    mutate(
      d2Adg2 = c(0, get_2dsmooth(.fit)$d2Adg2, 0),
      linear = abs(d2Adg2) < cutoff
    ) 
  
  # Remove small gaps of nonlinearity
  zero_positions = which(!d1$linear)
  diff_positions = diff(zero_positions)
  seq_starts = c(zero_positions[1], zero_positions[which(diff_positions > 1) + 1])
  seq_ends = c(zero_positions[which(diff_positions > 1)], zero_positions[length(zero_positions)])
  seq_lengths = seq_ends - seq_starts
  seq_starts[seq_lengths < min_gap]
  s = which(seq_lengths < min_gap) |>
    map(\(.x) seq_starts[.x]:seq_ends[.x]) |>
    flatten_int()
  d1$linear[s] = TRUE
  
  # Find longest linear piece
  one_positions = which(d1$linear)
  diff_positions = diff(one_positions)
  seq_starts = c(one_positions[1], one_positions[which(diff_positions > 1) + 1])
  seq_ends = c(one_positions[which(diff_positions > 1)], one_positions[length(one_positions)])
  
  seq_widths = log(d1$gsw[seq_ends]) - log(d1$gsw[seq_starts])
  seq_n = which.max(seq_widths)
  longest_linear = logical(nrow(d1))
  longest_linear[seq_starts[seq_n]:seq_ends[seq_n]] <- TRUE
  
  mutate(d1, longest_linear = longest_linear)
  
}

add_linear = function(.d, n_sim, ...) {
  
  # 1. fit linear regression
  fit1 = lm(A ~ log(gsw), data = .d)
  
  # 2. simulate N synthetic data sets from lm fit and then fit smooth function 
  # to each synthetic dataset
  syn = simulate_smooth(fit1, n_sim)
  
  # 3. Calculate second derivative cutoff based on null distribution
  cutoff = get_cutoff(syn, 0.99)
  
  # 4. Fit smooth function to data
  fit2 = gam(A ~ s(log(gsw)), data = .d)
  
  # 5. Add column on longest linear piece
  d2 = find_longest_linear(.d, fit2, cutoff, ...)
  
  d2
  
}

# 1. Simulate synthetic data ----
# n.b. Many specialized functions related to simulating and fitting LI6800 data 
# are in r/licor-functions.R

# Function to calculate the decay, per s, in autocorrelation between data points
calculate_corr_decay = function(rho, t) {
  assert_numeric(rho, lower = -1, upper = 1, any.missing = FALSE, min.len = 1L)
  assert_numeric(t, lower = 0, any.missing = FALSE, min.len = 1L)
  assert_true(length(t) %in% c(1L, length(rho)))
  # rho = exp(-b*t)
  # log(rho) = - b * t
  # -log(rho) / t = b
  -log(rho) / t
}

# Function to generate correlation matrix for simulating correlated error
make_autocorr_matrix = function(elapsed, b_autocorr) {
  # `elapsed` should be vector of time points in a series
  assert_numeric(elapsed, lower = 0, any.missing = FALSE, finite = TRUE)
  # `b_autocorr` should be a number indicating the autocorrelation decay rate
  assert_number(b_autocorr, lower = 0, finite = TRUE)
  
  m1 = outer(elapsed, elapsed, "-")
  R = diag(1, length(elapsed))
  
  ltri = exp(-b_autocorr * m1[which(lower.tri(m1))])
  utri = exp(-b_autocorr * t(m1)[which(upper.tri(m1))])
  
  R[which(lower.tri(R))] = ltri
  R[which(upper.tri(R))] = utri
 
  assert_true(isSymmetric(R))
  
  R
  
}

# 3. Fit data ----

# Calculate AA from parameter estimates:
# log(A_amphi) = b0_a + b1_a * log(gsw) + b2_a * log(gsw) ^ 2
# log(A_hypo)  = b0_h + b1_h * log(gsw) + b2_h * log(gsw) ^ 2
# Indefinite integral of log(A_amphi / A_hypo) when log(A) is a quadratic 
# function of log(gsw)
aa_int = function(log_gsw, b0_a, b0_h, b1_a, b1_h, b2_a, b2_h) {
  
  log_gsw ^ 3 * (b2_a / 3 - b2_h / 3) + 
    log_gsw ^ 2 * (b1_a / 2 - b1_h / 2) + 
    log_gsw * (b0_a - b0_h)
  
}
