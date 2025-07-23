rm(list = ls())

library(ape)
library(bayesplot)
library(brms)
library(checkmate)
library(cmdstanr)
library(cowplot)
library(dplyr)
library(forcats)
library(furrr)
library(ggplot2)
library(glue)
library(grid)
library(loo)
library(lubridate)
library(magrittr)
library(mgcv)
library(mvnfast)
library(phangorn)
library(progressr)
library(purrr)
library(readr)
library(rlang)
library(rootSolve)
library(scales)
library(stringr)
library(tibble)
library(tidybayes)
library(tidyr)
library(TreeTools)

source("r/functions.R")

# format of acceptable IDs
id_string = "^(LA[0-9]{4}A*)-[A-Z]{1}[AB]{0,1}$"

theme_set(theme_cowplot())

# options 
aa_args = read_rds("objects/aa_args.rds") |>
  c(
    aa_outlier_threshold = 3,
    thinning_interval = 0.05,
    n_iter_init = 2e3,
    max_divergent = 10,
    max_rhat = 1.01,
    min_ess = 400
  )

write_rds(aa_args, "objects/aa_args1.rds")

# Fit Rubisco model
x_depth = c(145, 830.5)
y_depth = c(0, 700)
x_rubisco = c(558, 56.5)
y_rubisco = c(0, 100)
nishio_carbon_1993_fig5 = read_csv("data/nishio_carbon_1993_fig5.csv",
                                   show_col_types = FALSE) |>
  mutate(
    depth = raw_depth * (diff(y_depth) / diff(x_depth)) - x_depth[1],
    rubisco = raw_rubisco * (diff(y_rubisco) / diff(x_rubisco)) -
      x_rubisco[1] * (diff(y_rubisco) / diff(x_rubisco)),
    # Reversing so that 0 is at abaxial, 1 is at adaxial, to match Earles et al.
    rel_depth = rev(depth / (min(depth) + max(depth)))
  )

fit_rubisco = gam(rubisco ~ s(rel_depth), data = nishio_carbon_1993_fig5)
