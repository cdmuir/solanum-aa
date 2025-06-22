# Fit phylogenetic model on posterior estimate of AA 
source("r/header.R")

d1 = read_rds("objects/stan_data_df.rds")
tr = read_rds("data/phylogeny.rds")
A = vcv(tr, corr = TRUE)

# based on loo_compare, no evidence for amphi_first being significant, so removed

fixed_effects = c(
  "1",
  "light_treatment",
  "light_intensity",
  "light_treatment + light_intensity",
  "light_treatment * light_intensity"
)

random_effects = c(
  "(1 | gr(acc, cov = A))",
  "(1 + light_treatment | gr(acc, cov = A))",
  "(1 + light_intensity | gr(acc, cov = A))",
  "(1 + light_treatment + light_intensity | gr(acc, cov = A))",
  "(1 + light_treatment * light_intensity | gr(acc, cov = A))"
)

sigma_effects = c(
  "1",
  "light_treatment",
  "light_intensity",
  "light_treatment + light_intensity"
)

# Cross all combinations
set.seed(622232590)
model_forms = expand.grid(fixed = fixed_effects,
                          random = random_effects,
                          sigma = sigma_effects,
                          stringsAsFactors = FALSE) %>%
  mutate(seed = sample(1e9, nrow(.)))

# Build and fit each model
plan(multisession, workers = 19)

aa_models = future_pmap(model_forms, function(fixed, random, sigma, seed) {
  fml = bf(as.formula(paste("aa ~", fixed, "+", random)), as.formula(paste("sigma ~", sigma)))
  brm(
    formula = fml,
    data = d1,
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
  aa_models,
  \(m) loo(
    m,
    resp = "aa",
    moment_match = FALSE,
    reloo = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

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

write_rds(aa_loo_table, "objects/aa_loo_table1.rds")
best_model_index = str_extract(aa_loo_table[1, "model"], "\\d+") |>
  as.integer()
write_rds(aa_models[[best_model_index]], "objects/fit_aa1.rds")
