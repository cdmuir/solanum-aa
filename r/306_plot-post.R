# Plot MCMC trace plots and compare prior to posterior
# NOT DONE - WAITING ON MORE FINAL MODELS TO FINISH RUNNING
source("r/header.R")

paste0("aa", 1:5) |>
  walk(\(model) {
    fit = read_rds(glue("objects/fit_{model}.rds"))
    
    focal_pars = get_par_table(model) |>
      dplyr::filter(length == "1")
    
    par_draws = fit$draws(focal_pars$parameter) |>
      as_draws_df()
    
    gp = focal_pars$parameter |>
      map(\(.p) {

        gp_trace = mcmc_trace(par_draws, .p)
        df_post = fit$draws(.p) |>
          as_draws_df() |>
          mutate(distribution = "posterior")
        df_prior = tibble(distribution = "prior",
                          p = seq(0.005, 0.995, length.out = 1000),)
        
        df_prior$x = focal_pars |>
          dplyr::filter(parameter == .p) |>
          pull(prior) |>
          str_replace("normal\\(", "qnorm(df_prior$p,") %>%
          parse(text = .) |>
          eval()
        
        df_prior[.p] = focal_pars |>
          dplyr::filter(parameter == .p) |>
          pull(prior) |>
          str_replace("normal\\(", "dnorm(df_prior$x,") %>%
          parse(text = .) |>
          eval()
        
        gp_dist = ggplot() +
          geom_density(
            aes_string(
              x = .p,
              color = "distribution",
              fill = "distribution"
            ),
            data = df_post,
            alpha = 0.5
          ) +
          geom_area(
            aes_string(
              x = "x",
              y = .p,
              color = "distribution",
              fill = "distribution"
            ),
            data = df_prior,
            alpha = 0.5
          ) +
          scale_color_manual(values = c(
            "posterior" = "steelblue",
            "prior" = "tomato"
          )) +
          scale_fill_manual(values = c(
            "posterior" = "steelblue",
            "prior" = "tomato"
          )) +
          labs(title = .p,
               y = "probability density",
               x = "parameter value")
        
        plot_grid(gp_trace, gp_dist, ncol = 1)
      })
    
    pdf(glue("figures/{model}_trace-dist.pdf"),
        width = 6,
        height = 4)
    gp
    dev.off()
  })