# Fit nonphylogenetic model on posterior estimate of AA to determine model form  for phylogenetic model ----
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
set.seed(602275201)
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
best_model_index = str_extract(aa_loo_table[1, "model"], "\\d+") |>
  as.integer()
write_rds(aa_models[[best_model_index]], "objects/fit_aa.rds")


df_new = crossing(
  acc = unique(d1$acc),
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000")
) |>
  mutate(row = row_number())

df_pred1 = posterior_epred(aa_models[[25]], newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row, names_to = "draw", values_to = "aa") |>
  full_join(df_new, by = "row") |>
  group_by(acc, light_treatment, light_intensity) |>
  point_interval(aa) |>
  group_by(light_treatment, light_intensity) |>
  arrange(aa) |>
  mutate(x = order(aa))
  
ggplot(df_pred1, aes(x, aa, ymin = .lower, ymax = .upper, group = acc)) +
  facet_grid(light_intensity ~ light_treatment) +
  geom_pointinterval() +
  geom_hline(yintercept = 0, linetype = "dashed") 

df_pred2 = posterior_epred(aa_models[[25]], newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row, names_to = "draw", values_to = "aa") |>
  full_join(df_new, by = "row") |>
  dplyr::select(-row) |>
  pivot_wider(names_from = light_treatment, values_from = aa) |>
  mutate(d_aa = high - low) |>
  group_by(acc, light_intensity) |>
  point_interval(d_aa) |>   
  group_by(light_intensity) |>
  arrange(d_aa) |>
  mutate(x = order(d_aa))

ggplot(df_pred2, aes(x, d_aa, ymin = .lower, ymax = .upper, group = acc)) +
  facet_grid(light_intensity ~ .) +
  geom_pointinterval() +
  geom_hline(yintercept = 0, linetype = "dashed") 
