# Check convergence and refit until criteria are met
source("r/header.r")

c("fit_aa1", "fit_aa2", "fit_stomata") |>
  walk(\(.m) {
    s = glue("objects/{.m}.rds")
    m = read_rds(s)
    
    crit = summarise_draws(m) |>
      summarize(c1 = (max(rhat, na.rm = TRUE) < aa_args$max_rhat),
                c2 = (min(ess_bulk, na.rm = TRUE) > aa_args$min_ess)) |>
      mutate(c3 = c1 * c2) |>
      pull(c3)
    
    while (crit == 0) {
      m = update(m,
                 iter = 4 * ndraws(m),
                 seed = m$fit@metadata$metadata$seed)
      crit = summarise_draws(m) |>
        summarize(c1 = (max(rhat, na.rm = TRUE) < aa_args$max_rhat),
                  c2 = (min(ess_bulk, na.rm = TRUE) > aa_args$min_ess)) |>
        mutate(c3 = c1 * c2) |>
        pull(c3)
    }
    
    write_rds(m, s)
    
  })
