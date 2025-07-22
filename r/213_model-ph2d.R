source("r/header.R")
source("r/photo_2d_pm.R")


plan(multisession, workers = 10)

handlers(global = TRUE)
handlers("cli")

# U / 2 %in% c(169, 53, 17) corresponds to approximately 10, 100, 1000 mm^2
# Stomatal density (number per mm^2)
# c(10, 100, 1000) |>
#   # convert to per um^2
#   multiply_by(1e-6) |>
#   D2U() |>
#   divide_by(2)

# Palisade porosity based on Lundgren et al. (2019), fig 5I

input = crossing(
  n_x = 53,
  n_z = c(101, 201, 401),
  stomatal_distribution = c("amphi", "hypo"),
  I_0 = round(seq(0.00005, 0.002, length.out = 10), 5),
  phi_pal = 0.1
) %>%
  # J_max prop to I_0
  mutate(J_max = I_0) %>%
  split(~ seq_len(nrow(.)))

run_photo_2d_pm = function(input, fit_rubisco) {
  
  p = progressor(along = input)
  out = future_map(input, \(.x, fit_rubisco) {
    
    r = with(.x, list(I_0 = I_0, n_x = n_x, n_z = n_z, phi_pal = phi_pal))
    ph2d = photo_2d_pm(.x$stomatal_distribution, fit_rubisco = fit_rubisco, replace = r)
    ph2d$parms$stomatal_distribution = .x$stomatal_distribution
    ph2d$parms$A_n = calc_2d_pm_An(ph2d)
    p()
    ph2d
    
  }, fit_rubisco = fit_rubisco)
  
  out
  
}

ph2d = run_photo_2d_pm(input, fit_rubisco)

write_rds(ph2d, "objects/ph2d.rds")
