# Copy up-to-data from adaptive-amphistomy data repo
source("r/header.R")

# For now, data on UH Google Drive.
# Change path to  "../../data/adaptive-amphistomy" after data collection is complete

c("plant-info", "stomata") |>
  walk(\(.x) {
    file.copy(
      glue(
        "/Users/cdmuir/Library/CloudStorage/GoogleDrive-cdmuir@hawaii.edu/Shared drives/muir-lab/adaptive-amphistomy/processed-data/{.x}.rds"
      ),
      glue("data/{.x}.rds")
    )
  })
