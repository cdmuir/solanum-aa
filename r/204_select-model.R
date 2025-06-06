# Fit nonphylogenetic model on posterior estimate of AA to determine model form  for phylogenetic model ----
source("r/header.R")

d1 = read_rds("objects/stan_data_df.rds")
# based on loo_compare, no evidence for amphi_first being significant, so removed

fixed_effects = c(
  "1",
  "light_treatment",
  "light_intensity",
  "light_treatment + light_intensity",
  "light_treatment * light_intensity"
)

random_effects = c(
  "(1 | acc)",
  "(1 + light_treatment | acc)",
  "(1 + light_intensity | acc)",
  "(1 + light_treatment + light_intensity | acc)",
  "(1 + light_treatment * light_intensity | acc)"
)

# Cross all combinations
set.seed(602275201)
model_forms = expand.grid(fixed = fixed_effects,
                          random = random_effects,
                          stringsAsFactors = FALSE) %>%
  mutate(seed = sample(1e9, nrow(.)))

# Build and fit each model
plan(multisession, workers = 19)

aa_models = future_pmap(model_forms, function(fixed, random, seed) {
  fml = as.formula(paste("aa ~", fixed, "+", random))
  brm(
    formula = fml,
    data = d1,
    backend = "cmdstanr",
    chains = 1,
    family = student(),
    save_pars = save_pars(all = TRUE),
    seed = seed
  ) |>
    add_criterion("loo", moment_match = TRUE, reloo = TRUE)
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

# Extract loo objects
loos = lapply(aa_models, function(m) m$criteria$loo)

# Name the models
names(loos) = paste0("model_", seq_along(loos))

# Create table
aa_loo_table = tibble(
  model = names(loos),
  looic = sapply(loos, \(x) x$estimates["looic", "Estimate"]),
  se = sapply(loos, \(x) x$estimates["looic", "SE"])
) |>
  arrange(looic) |>
  mutate(
    delta_looic = looic - min(looic),
    best_model = delta_looic == 0
  )

write_rds(aa_loo_table, "objects/aa_loo_table.rds")
