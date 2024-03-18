# Prepare actual data for Stan
source("r/header.R")

phy = read_rds("data/phylogeny.rds")

rh_curves = read_rds("data/prepared_rh_curves.rds") |>
  dplyr::select(
    A,
    scaled_log_gsw,
    S,
    acc,
    curve,
    acc_id,
    leaf_type,
    light_intensity,
    lightintensity_x_acc_id,
    light_treatment
  ) 

accession_climate = read_rds("data/accession-climate.rds") |>
  dplyr::select(acc1 = accession, ppfd_mol_m2) |>
  dplyr::filter(acc1 %in% unique(rh_curves$acc)) |>
  mutate(scaled_ppfd_mol_m2 = (ppfd_mol_m2 - mean(ppfd_mol_m2)) / sd(ppfd_mol_m2))

assert_true(setequal(unique(rh_curves$acc), accession_climate$acc1))
assert_true(setequal(phy$tip.label, accession_climate$acc1))

# make divergence matrix and index in same order as acc
Dmat = cophenetic(phy)
i = as.numeric(as.factor(colnames(Dmat)))
Dmat1 = Dmat[i,i] / max(Dmat)

# Compose data
stan_rh_curves = rh_curves |>
  compose_data()

# rh_curves$light_intensity[1:10]
# stan_rh_curves$light_intensity[1:10]
# rh_curves$light_treatment[1:10]
# stan_rh_curves$light_treatment[1:10]

  # commented variables are needed for full IRGA model if I go back to that
    # elapsed,
    # flow = Flow,
    # g_bw = gbw,
    # K,
    # P = Pa,
    # RH,
    # s = S,
    # T_air = Tair,
    # T_leaf = Tleaf,
    # CO2_r,
    # CO2_s,
    # H2O_r,
    # H2O_s
  
# drop variables that will be indexed by higher-level groups
stan_rh_curves$acc = NULL
stan_rh_curves$acc_id = NULL
stan_rh_curves$light_intensity = NULL
stan_rh_curves$light_treatment = NULL
stan_rh_curves$S = NULL

# Manual changes
stan_rh_curves = c(
  stan_rh_curves,
  
  ## Mapping of lightintensity_x_acc_id to curve
  rh_curves |>
    mutate(
      # character to numeric
      across(c("lightintensity_x_acc_id", "acc", "acc_id", "curve", "light_intensity"),
             \(.x) as.numeric(as.factor(.x))),
      # factor to numeric
      across(c("light_treatment"), as.numeric),
    ) |>
    summarize(
      curve = first(curve),
      acc = first(acc),
      acc_id = first(acc_id),
      light_intensity = first(light_intensity),
      light_treatment = first(light_treatment),
      .by = c(lightintensity_x_acc_id, leaf_type)
    ) |>
    pivot_wider(names_from = "leaf_type", values_from = "curve") |>
    arrange(lightintensity_x_acc_id) |>
    dplyr::select(amphi, pseudohypo, acc, acc_id, light_intensity, light_treatment) |>
    as.list(),
  
  ## min and max scaled_log_gsw, S by curve
  rh_curves |>
    mutate(across(c("curve"), \(.x) as.numeric(as.factor(.x)))) |>
    summarize(
      min_scaled_log_gsw = min(scaled_log_gsw),
      max_scaled_log_gsw = max(scaled_log_gsw),
      S = first(S),
      .by = c(curve)
    ) |>
    arrange(curve) |>
    dplyr::select(min_scaled_log_gsw, max_scaled_log_gsw, S) |>
    as.list()
  
)

stan_rh_curves$n_pts = rh_curves |>
  summarise(n = n(), .by = "curve") |>
  pull(n)

# add SPLASH data
stan_accession_climate = compose_data(accession_climate)
stan_accession_climate$n = NULL
stan_rh_curves = c(stan_rh_curves, stan_accession_climate)

# add phylogeny
stan_rh_curves$Dmat = Dmat1

assert_false(any(duplicated(names(stan_rh_curves))))

write_rds(stan_rh_curves, "data/stan_rh_curves.rds")
