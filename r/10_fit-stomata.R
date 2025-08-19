# Fit stomata model
source("r/header.R")

fit_aa = read_rds("objects/fit_aa1.rds")
d1 = fit_aa$data

tr = read_rds("data/phylogeny.rds")
A = vcv(tr, corr = TRUE)

plant_info = read_rds("data/plant-info.rds") |>
  dplyr::select(accession, replicate, light_treatment)

df_stomata = read_rds("data/stomata.rds") |>
  left_join(plant_info, by = join_by(accession, replicate)) |>
  rename(acc = accession) |>
  filter(acc %in% unique(d1$acc),
         if_all(contains("stomatal_density"), ~ !is.na(.x)),
         if_all(contains("guard_cell_length"), ~ !is.na(.x))) |>
  mutate(
    lower_sd = log(lower_stomatal_density_mm2),
    upper_sd = log(upper_stomatal_density_mm2),
    lower_gcl = log(lower_guard_cell_length_um),
    upper_gcl = log(upper_guard_cell_length_um)
  )

fixed_effects = c("1", "light_treatment")

random_effects = c("(1 | gr(acc, cov = A))", "(1 + light_treatment | gr(acc, cov = A))")

sigma_effects = c("1", "light_treatment")

# Cross all combinations
set.seed(704392913)
model_forms = expand.grid(
  fixed = fixed_effects,
  random = random_effects,
  sigma = sigma_effects,
  stringsAsFactors = FALSE
) %>%
  mutate(seed = sample(1e9, nrow(.)))

# Build and fit each model
plan(multisession, workers = 19)

stomata_models = future_pmap(model_forms, function(fixed, random, sigma, seed) {
  fml = bf(as.formula(
    paste(
      "mvbind(lower_sd, upper_sd, lower_gcl, upper_gcl) ~",
      fixed,
      "+",
      random
    )
  ), as.formula(paste("sigma ~", sigma)))
  brm(
    formula = fml + set_rescor(TRUE),
    data = df_stomata,
    data2 = list(A = A),
    backend = "cmdstanr",
    chains = 1,
    family = student(),
    save_pars = save_pars(all = TRUE),
    seed = seed
  )
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

# Extract loo objects
loos = future_map(
  stomata_models,
  \(m) loo(
    m,
    # resp = "aa",
    moment_match = FALSE,
    reloo = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

# Name the models
names(loos) = paste0("model_", seq_along(loos))

# Create table
stomata_loo_table = tibble(
  model = names(loos),
  looic = sapply(loos, \(x) x$estimates["looic", "Estimate"]),
  se = sapply(loos, \(x) x$estimates["looic", "SE"])
) |>
  arrange(looic) |>
  mutate(delta_looic = looic - min(looic), best_model = delta_looic == 0)

write_rds(stomata_loo_table, "objects/stomata_loo_table.rds")
best_model_index = str_extract(stomata_loo_table[1, "model"], "\\d+") |>
  as.integer()
write_rds(stomata_models[[best_model_index]], "objects/fit_stomata.rds")
