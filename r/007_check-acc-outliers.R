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
  chains = 1
)

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
