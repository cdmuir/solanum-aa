rm(list = ls())

library(bayesplot)
library(brms)
library(checkmate)
library(cmdstanr)
library(cowplot)
library(dplyr)
library(ggplot2)
library(glue)
library(lubridate)
library(magrittr)
library(purrr)
library(readr)
library(scales)
library(stringr)
library(tibble)
library(tidybayes)
library(tidyr)

source("r/functions.R")
source("r/licor-functions.R")

# format of acceptable IDs
id_string = "^(LA[0-9]{4}A*|nelsonii|sandwicense)-[A-Z]{1}[A]{0,1}$"

theme_set(theme_cowplot())
