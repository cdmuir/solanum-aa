rm(list = ls())

library(ape)
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
library(phangorn)
library(posterior)
library(purrr)
library(readr)
library(rlang)
library(rstan)
library(scales)
library(stringr)
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
