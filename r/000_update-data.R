# Copy up-to-data from adaptive-amphistomy data repo
source("r/header.R")

c("processed-data/accession-climate",
  "processed-data/accession-gedi",
  "processed-data/plant-info",
  "filtered-data/rh_curves",
  "filtered-data/stomata",
  "objects/df_germ_summary",
  "objects/df_growth_summary") |>
  walk(\(.x) {
    file.copy(
      glue("../../data/adaptive-amphistomy/{.x}.rds"),
      glue("data/{.y}.rds", .y = str_remove(.x, "^.*/")),
      overwrite = TRUE
    )
  })

c("objects/aa_args", "objects/aa_stats") |>
  walk(\(.x) {
    file.copy(
      glue("../../data/adaptive-amphistomy/{.x}.rds"),
      glue("{.x}.rds"),
      overwrite = TRUE
    )
  })
