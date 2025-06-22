# Effect of LMA on AA. Work in progress.
source("r/header.r")

plant_info = read_rds("data/plant-info.rds")

plant_info |>
  filter(lma_gm2 < 20) 

plant_info |>
  arrange(desc(lma_gm2)) |>
  select(acc_id, lma_gm2) |>
  print(n = 100)
  
aa_summary = read_rds("objects/aa_summary.rds")

df1 = plant_info |>
  select(acc_id, light_treatment, lma_gm2) |>
  filter(!is.na(lma_gm2)) |>
  full_join(aa_summary, by = join_by(acc_id)) |>
  filter(!is.na(lma_gm2), !is.na(median)) 

df1 |>
  filter(lma_gm2 < 20)
  pull(lma_gm2) |> hist()
  mutate(acc = str_extract(acc_id, "(LA[0-9]{4}A*|nelsonii|sandwicense)")) |>
  summarize(
    aa = median(median),
    lma = median(lma_gm2),
    .by = c(acc, light_treatment, light_intensity)
  )

ggplot(df1, aes(x = lma, y = aa, color = light_treatment, shape = light_intensity)) +
  geom_point() +
  geom_smooth(method = "lm")
