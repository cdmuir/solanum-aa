# Modify functions to save checkpoints for Stan models
# this should eventually be moved to functions.R
source("r/header.R")
library(chkptstanr)

# for debuggin
model_code = model_code = read_lines("stan/solanum-aa1.stan") |> paste(collapse = "\n")
data = read_rds("data/stan_rh_curves.rds")
iter_warmup = 40
iter_sampling = 40
thin = 1
iter_per_chkpt = 10
chkpt_progress = TRUE
path = "chkpt_aa1"
init = read_rds("objects/init.rds")
max_treedepth = 10L

# chkpt_stan1 = function(
#     model_code,
#     data,
#     iter_warmup = 1000,
#     iter_sampling = 1000,
#     thin = 1,
#     iter_per_chkpt = 100,
#     chkpt_progress = TRUE,
#     path, 
#     init = NULL,
#     max_treedepth = 10L
# ) {
#   
  # assert_true(iter_warmup %% iter_per_chkpt == 0, "iter_warmup must be a multiple of iter_per_chkpt")
  # assert_true(iter_sampling %% iter_per_chkpt == 0, "iter_sampling must be a multiple of iter_per_chkpt")
  # assert_true(iter_per_chkpt < iter_warmup, "iter_per_chkpt must be less than iter_warmup")
  # assert_true(iter_per_chkpt < iter_sampling, "iter_per_chkpt must be less than iter_sampling")
  if (!dir.exists(path)) {
    create_folder(folder_name = path)
  }
  
  ## write model code
  stan_code_path = cmdstanr::write_stan_file(
    code = model_code,
    dir = paste0(path, "/stan_model"),
    basename = "model"
  )
  
  ## compile model
  m = cmdstan_model(stan_code_path, dir = path)
  
  ## calculate number of checkpoints
  n_chkpts = ceiling(iter_sampling / iter_per_chkpt)
  
  # check on current status
  cp_files = list.files(paste0(path, "/cp_info"))
  
  if (length(cp_files) == 0) {
    # initiate fit
    m$sample(
      data = data,
      refresh = 0,
      init = init,
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
# WORKING THROUGH OPTIONS HERE
step_size = NULL,
metric = NULL,
metric_file = NULL,
inv_metric = NULL,
init_buffer = NULL,
term_buffer = NULL,
window = NULL,
fixed_param = FALSE,
show_messages = TRUE,
show_exceptions = TRUE,
diagnostics = c("divergences", "treedepth", "ebfmi")
    )
    ## write model to file
    ## compile model
  } else {
    # continue fitting
  }
}
