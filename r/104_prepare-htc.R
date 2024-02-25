# Prepare file fit Stan models on htc
source("r/header.R")

file.copy("stan/solanum-aa2.stan", "htc/solanum-aa2.stan", overwrite = TRUE)
list.files("synthetic-data", pattern = "stan_sim[0-9]{4}.rds", full.names = TRUE) |>
  walk(\(.x) {
    file.copy(.x, str_replace(.x, "synthetic-data", "htc"), overwrite = TRUE)
    .n = str_extract(.x, "[0-9]+")
    
    read_lines("htc-templates/template.R") |>
      map_chr(glue, n = .n) |>
      write_lines(glue("htc/fit_{n}.R", n = .n))
    
    read_lines("htc-templates/template.sh") |>
      map_chr(glue, n = .n) |>
      write_lines(glue("htc/fit_{n}.sh", n = .n))
    
    read_lines("htc-templates/template.sub") |>
      map_chr(glue, n = .n) |>
      write_lines(glue("htc/fit_{n}.sub", n = .n))
    
  }) 
