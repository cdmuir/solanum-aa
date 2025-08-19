# install_packages.R

# Combine both sets of packages and remove duplicates
pkgs <- unique(c(
  "ape",
  "brms",
  "checkmate",
  "cmdstanr",
  "cowplot",
  "dplyr",
  "forcats",
  "furrr",
  "ggplot2",
  "glue",
  "loo",
  "lubridate",
  "magrittr",
  "phangorn",
  "posterior",
  "purrr",
  "readr",
  "rlang",
  "rstan",
  "scales",
  "stringr",
  "tidybayes",
  "tidyr",
  "TreeTools",
  "kableExtra",
  "knitr",
  "tibble"
))

# Install any missing packages from CRAN
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Run install
invisible(lapply(pkgs, install_if_missing))

# Special handling for cmdstanr (on CRAN but often recommended via GitHub for latest)
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
}

message("All requested packages are installed.")