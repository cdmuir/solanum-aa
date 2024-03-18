# Check model convergence
source("r/header.R")

paste0("aa", 1:5) |>
  walk(\(model) {
    fit = read_rds(glue("objects/fit_{model}.rds"))
    
    focal_pars = get_par_table(model) |>
      dplyr::filter(length == "1") |>
      pull(parameter)
    
    par_summary = fit$summary(focal_pars)
    
    assert_true(all(par_summary$rhat < 1.01))
    assert_true(all(par_summary$ess_tail > 1e3))
    
  })
