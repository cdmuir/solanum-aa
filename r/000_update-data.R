# Copy up-to-data from adaptive-amphistomy data repo
source("r/header.R")

c("processed-data/plant-info",
  "processed-data/rh_curves",
  "filtered-data/stomata") |>
  walk(\(.x) {
    file.copy(
      glue("../../data/adaptive-amphistomy/{.x}.rds"),
      glue("data/{.y}.rds", .y = str_remove(.x, "(process|filter)ed-data/"))
    )
  })
