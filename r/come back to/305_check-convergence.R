# Check model convergence (MIGHT NEED TO EDIT TO INCLUDE TERMS FOR RH CURVE TERMS)
# This scripts needs some adjusting after models finish running

source("r/header.R")

# args = commandArgs(trailingOnly = TRUE)
args = list("aa1", "0")

m = args[1]
chain = args[2]
chkpt_path1 = glue("checkpoints/chkpt_folder_{m}_{chain}")

# use this once model finishes running
# draws = read_rds(draws, glue("objects/draws_{m}_{chain}.rds"))
draws = combine_chkpt_draws1(path = chkpt_path1) # temporary version

focal_pars = get_par_table(m) |>
  dplyr::filter(type == "real") |>
  pull(parameter)

s = posterior::summarize_draws(draws)
s |>
  dplyr::filter(variable %in% focal_pars)

par_summary = draws$summary(focal_pars)

paste0("aa", 1:5) |>
  walk(\(model) {
    fit = read_rds(glue("objects/fit_{model}.rds"))
    
    focal_pars = get_par_table(model) |>
      dplyr::filter(type == "real") |>
      pull(parameter)
    
    par_summary = fit$summary(focal_pars)
    
    assert_true(all(par_summary$rhat < 1.01))
    assert_true(all(par_summary$ess_tail > 1e3))
    
  })
