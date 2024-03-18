# Prepare file fit Stan models on htc
source("r/header.R")

# Synthetic data (not working right now) ----
# file.copy("stan/solanum-aa3.stan", "htc/solanum-aa3.stan", overwrite = TRUE)
# list.files("synthetic-data", pattern = "stan_sim[0-9]{4}.rds", full.names = TRUE) |>
#   walk(\(.x) {
#     file.copy(.x, str_replace(.x, "synthetic-data", "htc"), overwrite = TRUE)
#     .e = environment()
#     .n = str_extract(.x, "[0-9]+")
#     
#     read_lines("htc-templates/template.R") |>
#       map_chr(glue, n = .n, .envir = .e) |>
#       write_lines(glue("htc/fit_{n}.R", n = .n, .envir = .e))
#     
#     read_lines("htc-templates/template.sh") |>
#       map_chr(glue, n = .n, .envir = .e) |>
#       write_lines(glue("htc/fit_{n}.sh", n = .n, .envir = .e))
#     
#     read_lines("htc-templates/template.sub") |>
#       map_chr(glue, n = .n, .envir = .e) |>
#       write_lines(glue("htc/fit_{n}.sub", n = .n, .envir = .e))
#     
#   }) 

list.files("stan", pattern = "solanum-aa[0-9]+", full.names = TRUE) |>
  str_extract("aa[0-9]+") |>
  walk(\(model) {
    c("R", "sh", "sub") |>
      map(\(extension, m) {
        .e = environment()
        paste0("htc-templates/fit_aa.", extension) |>
          read_lines() |>
          map_chr(glue, m = m, .envir = .e) |>
          write_lines(glue(
            "htc/fit_{m}.{ext}",
            m = m,
            ext = extension,
            .envir = .e
          ))
      }, m = model)
    
  })
