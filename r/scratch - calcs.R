# simplified back-of the envelope calculation for how much water could be save by switching to amphistomy
source("r/header.R")

rh_curves = read_rds("data/trimmed_rh_curves.rds")
tr = read_rds("data/phylogeny.rds")
A = vcv(tr, corr = TRUE)

plan(multisession, workers = 19)

set.seed(125691028)
tibble(acc = unique(rh_curves$acc), seed = sample(1e9, length(unique(rh_curves$acc)))) |>
  future_pwalk(function(acc, seed) {
    .acc = acc
    crit = 0
    x = 1
    
    while (crit == 0 & x < 24) {
      m = brm(
        formula = log_A ~ leaf_type + light_intensity * log_gsw + (light_intensity * log_gsw |
                                                                     acc_id),
        
        data = filter(rh_curves, acc == .acc),
        backend = "cmdstanr",
        chains = 1,
        iter = x * aa_args$n_iter_init,
        thin = x,
        family = student(),
        save_pars = save_pars(all = TRUE),
        seed = seed
      )
      
      crit = summarise_draws(m) |>
        summarize(c1 = (max(rhat, na.rm = TRUE) < aa_args$max_rhat),
                  c2 = (min(ess_bulk, na.rm = TRUE) > aa_args$min_ess)) |>
        mutate(
          n_divergent = nuts_params(m) |>
            subset(Parameter == "divergent__") |>
            pull(Value) |>
            sum(),
          c3 = (n_divergent < aa_args$max_divergent),
          c4 = c1 * c2 * c3
        ) |>
        pull(c4)
      
      x = x + 1
      
    }
    
    write_rds(m, glue("objects/simple-Ags-fits/{acc}.rds"))
    return(m)
  },
  .progress = TRUE,
  .options = furrr_options(seed = TRUE))
