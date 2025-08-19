
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

## Functions to write Stan models ----
write_solanum_aa = function(model) {
  assert_string(model)
  assert_choice(model, paste0("aa", 1:5))
  
  write_lines(
    c(
      solanum_aa_functions(),
      solanum_aa_data(),
      solanum_aa_parameters(model),
      solanum_aa_model(model),
      solanum_aa_generated(model)
    ),
    glue("stan/solanum-{model}.stan")
  )
  
  invisible()
  
}

get_par_table = function(model) {
  read_csv("data/parameters.csv", col_types = "cciddciiiii") |>
    pivot_longer(cols = starts_with("aa"),
                 names_to = "model",
                 values_to = "value") |>
    dplyr::filter(model == !!model, value == 1L)
}

solanum_aa_functions = function() {
  read_lines("stan/solanum-aa-functions.stan")
}

solanum_aa_data = function() {
  read_lines("stan/solanum-aa-data.stan")
}

solanum_aa_parameters = function(model) {
  
  par_table = get_par_table(model)
  
  # parameters
  pars = par_table |>
    mutate(
      lb1 = if_else(is.na(lb), "", glue("lower={lb}")),
      ub1 = if_else(is.na(ub), "", glue("upper={ub}")),
      sep = if_else(is.na(lb) | is.na(ub), "", ", "),
      bnds = if_else(bounds == "0", "", glue("<{lb1}{sep}{ub1}>")),
      par = glue("{type}{bnds} {parameter};")
    ) |>
    pull(par)
  
  # transformed parameters
  # only works with log_sigma -> sigma transformation
  tpars1 = par_table |>
    dplyr::filter(str_detect(parameter, "^log_sigma")) |>
    mutate(
      new_parameter = str_remove(parameter, "log_"),
      tpar = glue("{type} {new_parameter};"),
      trans = glue("{new_parameter} = exp({parameter});")
    ) |>
    dplyr::select(tpar, trans) |>
    pivot_longer(everything()) |>
    mutate(name = factor(name, levels = c("tpar", "trans"))) |>
    arrange(name, value) |>
    pull(value)
  
  # Covariance matrix from correlation matrix
  tpars2 = par_table |>
    dplyr::filter(str_detect(parameter, "^R_")) |>
    mutate(
      new_parameter = str_replace(parameter, "R_", "Sigma_"),
      type1 = str_remove(type, "corr_") |>
        str_replace("\\[([0-9]+)\\]", "[\\1,\\1]"),
      tpar = glue("{type1} {new_parameter};"),
      trans = glue("{new_parameter} = quad_form_diag({parameter}, {sigma});", 
                   sigma = str_replace(parameter, "R_", "sigma_"))
    ) |>
    dplyr::select(tpar, trans) |>
    pivot_longer(everything()) |>
    mutate(name = factor(name, levels = c("tpar", "trans"))) |>
    arrange(name, value) |>
    pull(value)
  
  tpars = c(tpars1, tpars2)
  
  c(
    glue("parameters {{
  {str_c(pars, collapse = '\n  ')}
}}"),
    glue("transformed parameters {{
  {str_c(tpars, collapse = '\n  ')}
}}")
  )
  
}

solanum_aa_model = function(model) {
  
  c(
    "model {",
    solanum_aa_priors(model),
    solanum_aa_aa(model),
    solanum_aa_rh(),
    solanum_aa_pai(model),
    "}"
  )
  
}

solanum_aa_priors = function(model) {
  
  par_table = get_par_table(model)
  
  # identify phylogentically structured parameters
  phylo_s = "^b_(.*)_acc$"
  phylo_pars = par_table |>
    dplyr::filter(str_detect(parameter, phylo_s)) |>
    pull(parameter)
  
  # priors on phylogenetic structure
  x = str_replace(phylo_pars, phylo_s, "\\1")
  phylo_priors = c(
    glue("  matrix[n_acc,n_acc] Sigma_{x}_acc;"),
    glue("  Sigma_{x}_acc = cov_GPL1(Dmat, etasq_{x}_acc, rhosq_{x}_acc, 0);")
  )
  
  # priors
  priors = par_table |>
    mutate(prior = glue("  {parameter} ~ {prior};")) |>
    pull(prior)
  
  c(
    "  // priors on phylogenetic structure",
    glue("  {str_c(phylo_priors, collapse = '\n  ')}"),
    "",
    "  // priors",
    glue("  {str_c(priors, collapse = '\n  ')}"),
    ""
  )
  
}

solanum_aa_aa = function(model) {
  
  b_2000 = switch(
    model,
    aa1 = "b_aa_light_intensity_2000",
    aa2 = "b_aa_light_intensity_2000",
    aa3 = "b_aa_light_intensity_2000 + b_aa_light_intensity_2000_acc[acc[i]]",
    aa4 = "b_aa_light_intensity_2000",
    aa5 = "b_aa_light_intensity_2000 + b_aa_light_intensity_2000_acc[acc[i]]"
  )
  
  b_high = switch(
    model,
    aa1 = "b_aa_light_treatment_high",
    aa2 = "b_aa_light_treatment_high",
    aa3 = "b_aa_light_treatment_high",
    aa4 = "b_aa_light_treatment_high + b_aa_light_treatment_high_acc[acc[i]]",
    aa5 = "b_aa_light_treatment_high + b_aa_light_treatment_high_acc[acc[i]]"
  )
  
  b_aa_2000_high = switch(
    model,
    aa1 = "",
    aa2 = "b_aa_2000_high * (light_intensity[i] == 2) * (light_treatment[i] == 2) +",
    aa3 = "",
    aa4 = "",
    aa5 = "b_aa_2000_high * (light_intensity[i] == 2) * (light_treatment[i] == 2) +"
  )
  
  read_lines("stan/solanum-aa-aa.stan") |>
    str_replace_all(
      c(
        "\\{b_2000\\}" = b_2000,
        "\\{b_high\\}" = b_high,
        "\\{b_aa_2000_high\\}" = b_aa_2000_high
      )
    ) 
  
}

solanum_aa_rh = function() {
  read_lines("stan/solanum-aa-rh.stan")
}

solanum_aa_pai = function(model) {
  
  aa_acc = switch(
    model,
    aa1 = "b0_aa + b_aa_acc;",
    aa2 = "b0_aa + b_aa_acc;",
    aa3 = "b0_aa + b_aa_acc + b_aa_light_intensity_2000 + b_aa_light_intensity_2000_acc;",
    aa4 = "b0_aa + b_aa_acc + b_aa_light_treatment_high + b_aa_light_treatment_high_acc;",
    aa5 = "b0_aa + b_aa_acc + b_aa_light_intensity_2000 + b_aa_light_intensity_2000_acc + b_aa_light_treatment_high + b_aa_light_treatment_high_acc;"
  )
  
  read_lines("stan/solanum-aa-pai.stan") |>
    str_replace("\\{aa_acc\\}", aa_acc)
  
}

solanum_aa_generated = function(model) {
  
  c(
    "generated quantities {",
    "  // calculated log-likelihood to estimate LOOIC for model comparison",
    "  vector[n_lightintensity_x_acc_id] log_lik;",
    str_replace(solanum_aa_aa(model), "target \\+", "log_lik[i] "),
    "}"
  )
  
}

## Functions for checkpointing Stan models ----
inititalize_stan = function(
    model_code,
    data,
    iter_typical = 150,
    path,
    init = NULL,
    max_treedepth = 10L
) {
  
  path = chkptstanr::create_folder(folder_name = path)
  
  ## write model code
  stan_code_path = cmdstanr::write_stan_file(
    code = model_code,
    dir = paste0(path, "/stan_model"),
    basename = "model"
  )
  
  ## compile model
  m = cmdstan_model(stan_code_path, dir = path)
  
  ## initiate fit
  sample_chunk = m$sample(
    data = data,
    refresh = 0,
    init = list(init),
    output_dir = path,
    output_basename = "model",
    chains = 1L,
    parallel_chains = 1L,
    iter_warmup = iter_typical,
    iter_sampling = 0,
    save_warmup = TRUE,
    thin = 1e0,
    max_treedepth = max_treedepth,
    adapt_engaged = TRUE,
    adapt_delta = 0.8,
    step_size = NULL
  )
  
  stan_state = chkptstanr::extract_stan_state(sample_chunk, "warmup")
  
  write_rds(data, paste0(path, "/data.rds"))
  write_rds(stan_state, paste0(path, "/cp_info/cp_info_0.rds"))
  
  invisible()
  
}

## Modified from chkptstanr::chkpt_stan()
chkpt_stan1 = function(path,
                       iter_warmup = 1000,
                       iter_sampling = 1000,
                       thin = 1,
                       iter_per_chkpt = 100,
                       chkpt_progress = TRUE,
                       max_treedepth = 10L) {
  chkpt_set_up = chkptstanr::chkpt_setup(iter_sampling, iter_warmup, iter_per_chkpt)
  
  data = readRDS(file = paste0(path, "/data.rds"))
  
  ## path to model code
  stan_code_path = paste0(path, "/stan_model/model.stan")
  
  ## compile model
  m = cmdstan_model(stan_code_path, dir = path)
  
  ## check on current status
  cp_files = list.files(paste0(path, "/cp_info"))
  
  ## get stan state
  checkpoints = as.numeric(gsub(".*info_(.+).rds.*", "\\1", cp_files))
  last_chkpt = max(checkpoints)
  stan_state = readRDS(file = paste0(path, "/cp_info/",
                                     cp_files[which.max(checkpoints)]))
  # fix NAs in stan_state$init
  stan_state$inits = stan_state$inits |> 
    map(\(.chain) {
      .chain |>
        map(\(.par) {
          .par = ifelse(is.na(.par), 0, .par)
        })
    })
  
  if (last_chkpt == chkpt_set_up$total_chkpts) {
    return(message("Checkpointing complete"))
  } else {
    cp_seq = seq(last_chkpt + 1, chkpt_set_up$total_chkpts)
  }
  
  for (i in cp_seq) {
    if (i <= chkpt_set_up$warmup_chkpts) {
      stan_phase = "warmup"
      sample_chunk = m$sample(
        data = data,
        refresh = 0,
        init = stan_state$inits,
        output_dir = path,
        output_basename = "model",
        chains = 1L,
        parallel_chains = 1L,
        iter_warmup = iter_per_chkpt,
        iter_sampling = 0,
        save_warmup = TRUE,
        thin = thin,
        max_treedepth = max_treedepth,
        adapt_engaged = TRUE,
        adapt_delta = 0.8,
        step_size = stan_state$step_size_adapt,
        inv_metric = stan_state$inv_metric
      )
    }
    else {
      stan_phase = "sample"
      sample_chunk = m$sample(
        data = data,
        refresh = 0,
        init = stan_state$inits,
        output_dir = path,
        output_basename = "model",
        chains = 1L,
        parallel_chains = 1L,
        iter_warmup = 0,
        iter_sampling = iter_per_chkpt,
        save_warmup = FALSE,
        thin = thin,
        max_treedepth = max_treedepth,
        adapt_engaged = FALSE,
        step_size = stan_state$step_size_adapt,
        inv_metric = stan_state$inv_metric
      )
      sample_chunk$save_object(paste0(path, "/cp_samples/",
                                      "samples_", i, ".rds"))
    }
    stan_state = chkptstanr::extract_stan_state(sample_chunk, stan_phase)
    saveRDS(object = stan_state,
            file = paste0(path, "/cp_info/cp_info_",
                          i, ".rds"))
    cat(chkptstanr:::progress_bar(chkpt_set_up, phase = stan_phase, i = i))
    if (i == chkpt_set_up$total_chkpts) {
      message("Checkpointing complete")
    }
  }
  
}

sort_by_number = function(.x) {
  .x[.x |>
       str_remove_all(".*samples_") |>
       str_remove(".rds") |>
       as.integer() |>
       order()]
}

combine_chkpt_draws1 = function(path) {
  
  draws = list.files(paste0(path, "/cp_samples"), full.names = TRUE) |> 
    sort_by_number() |>
    map(read_rds) |>
    map(\(.x) .x$draws())
  
  draws_dim = dim(draws[[1]])
  draws_per_cp = draws_dim[1]
  chains = draws_dim[2]
  draws_total = draws_per_cp * length(draws)
  param_names = dimnames(draws[[1]])$variable
  param_total = length(param_names)
  draws_abind = abind::abind(draws, along = 1)
  draws_array = array(
    data = 0,
    dim = c(draws_total, chains,
            param_total),
    dimnames = list(
      iteration = c(1:draws_total),
      chain = c(1:chains),
      variable = param_names
    )
  )
  for (i in seq_along(param_names)) {
    draws_array[, , i] = draws_abind[, , param_names[i]]
  }
  class(draws_array) = c("draws_array", "draws", "array")
  return(draws_array)
  
}
