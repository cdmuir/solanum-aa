# Filter data to linear portions of each A-gsw curve
source("r/header.R")

# Example data and arguments
rh_curves = read_rds("data/rh_curves.rds") |>
  filter(acc_id == "LA0716-G", assumed_K == 0.5)

ggplot(rh_curves, aes(gsw, A)) +
  geom_point() +
  scale_x_log10()

d1 = rh_curves |>
  split(~ curve_type + light_intensity) |>
  map_dfr(add_linear, n_sim = 1e3)

ggplot(d1, aes(gsw, A, color = linear, shape = curve_type)) +
  geom_point() +
  scale_x_log10()

tmp = rh_curves |>
  filter(curve_type == "2-sided RH", light_intensity == "2000") 

lm(A ~ log(gsw), data = tmp)
plot(fit)

plot(tmp$A)
# STUFF BELOW NEEDS TO BE MOVED ELSEWHERE?
# Get error autocorrelation
# 
# 
get_ac1 = function(x) {
  cor(x[1:(length(x) - 1)], x[2:length(x)])
}

rh_curves |>
  filter(acc_id == "LA2172-AA") |>
  reframe(resid = resid(lm(A ~ log(gsw))),
          .by = c("acc", "acc_id", "light_treatment", "light_intensity", "curve_type")) |>
  summarize(r = get_ac1(resid),
            .by = c("acc", "acc_id", "light_treatment", "light_intensity", "curve_type"))

ggplot(filter(rh_curves, acc_id == "LA2172-AA"), aes(gsw, A, color = curve_type)) +
  geom_point() +
  scale_x_log10()
