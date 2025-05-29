# Check for accessions that are outliers
source("r/header.R")

aa_post = read_rds("objects/aa_post.rds")

aa_summary1 = aa_post |>
  group_by(acc, acc_id, light_treatment,light_intensity) |>
  point_interval(aa)

fit1 = brm(
  bf(
    aa ~ light_treatment * light_intensity + (light_treatment * light_intensity |
                                                acc) + (1 | acc_id),
    sigma ~ light_treatment * light_intensity
  ),
  family = student,
  data = aa_summary1,
  backend = "cmdstanr",
  iter = 2000,
  warmup = 2000,
  chains = 1,
  seed = 125409099
)


# Extract posterior draws of group-level effects
draws_re = fit1 |>
  as_draws_df() |>
  # random intercepts by group
  select(glue("r_acc_id[{acc_id},Intercept]", acc_id = unique(fit1$data$acc_id)))  

sigma_acc_id <- as_draws_df(fit1)$`sd_acc_id__Intercept`

# Get group names
group_names <- names(draws_re)

# For each posterior draw, simulate a value from the group-level distribution
simulated_re <- matrix(
  rnorm(n = length(sigma_acc_id) * length(group_names), mean = 0, sd = rep(sigma_acc_id, each = length(group_names))),
  ncol = length(group_names)
)

colnames(simulated_re) <- group_names

# Get posterior mean of each actual group's effect
observed_means <- colMeans(draws_re)

# For each group, calculate tail probability (Bayesian p-value)
tail_probs <- map_dbl(seq_along(group_names), function(i) {
  group <- group_names[i]
  sims <- simulated_re[, i]
  obs <- observed_means[i]
  
  # One-tailed p-value for being in tail of predictive distribution
  mean(abs(sims) > abs(obs))
})

# Combine into a tibble
outlier_summary <- tibble(
  group = str_extract(group_names, "(?<=\\[).+?(?=\\])"),
  observed = observed_means,
  p_extreme = tail_probs
) |>
  arrange(p_extreme)

outlier_summary

ranef(fit1)$acc_id[,,"Intercept"][,"Estimate"] |>
  qqnorm()
ranef(fit1)$acc_id[,,"Intercept"][,"Estimate"] |>
  qqline()

# scratch
aa_summary2 = aa_summary1 |>
  group_by(acc, light_treatment,light_intensity) |>
  summarize(aa = mean(aa))

ggplot(aa_summary1, aes(light_treatment, aa, color = light_intensity)) +
  geom_jitter(width = 0.1) +
  geom_boxplot()
fit1 = lm(aa ~ acc * light_treatment * light_intensity, data = aa_summary)
summary(aov(fit1))
plot(fit1)
hist(resid(fit1))
sigma(fit1)
