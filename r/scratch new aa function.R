# FIGURE OUT FUNCTION TO PULL OUT REGION FOR FITTING
source("r/header.R")

rh_curves = read_rds("data/rh_hi_curves.rds")

df1 = filter(rh_curves, acc_id == "nelsonii-S", light_intensity == "2000")

ggplot(df1, aes(gsw, A)) +
  geom_point()

# FUNCTION STARTS HERE
# function(df1, min_n = 20, ...) {
ranges = df1 |>
  summarize(
    min_gsw = min(gsw, na.rm = TRUE),
    max_gsw = max(gsw, na.rm = TRUE),
    .by = "curve_type"
  )

df2 = filter(df1, gsw >= max(ranges$min_gsw) &
               gsw <= min(ranges$max_gsw))

# Check for minimum number of points
summarize(df2, n = n(), .by = curve_type)

ggplot(df2, aes(gsw, A)) +
  geom_point()
