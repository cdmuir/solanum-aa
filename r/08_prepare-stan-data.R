# Prepare actual data for Stan
source("r/header.R")

plant_info = read_rds("data/plant-info.rds") |>
  mutate(llma = log(lma_gm2)) |>
  dplyr::select(accession, acc_id, light_treatment, llma)

phy = read_rds("data/phylogeny.rds")

aa_summary_steadystate = read_rds("objects/aa_summary_sty.rds") |>
  dplyr::select(acc_id, light_intensity, aa = median, se_aa = sd) |>
  mutate(scaled_aa = (aa - mean(aa)) / sd(aa)) |>
  left_join(plant_info, by = join_by(acc_id))

aa_summary_dynamic = read_rds("objects/aa_summary_dyn.rds") |>
  dplyr::select(acc_id, light_intensity, aa = median, se_aa = sd) |>
  mutate(scaled_aa = (aa - mean(aa)) / sd(aa)) |>
  left_join(plant_info, by = join_by(acc_id))

accession_gedi = read_rds("data/accession-gedi.rds") |>
  dplyr::select(accession, pai) |>
  dplyr::filter(accession %in% unique(aa_summary_steadystate$accession)) |>
  mutate(scaled_pai = (pai - mean(pai)) / sd(pai))

assert_true(setequal(
  unique(aa_summary_steadystate$accession),
  accession_gedi$accession
))
assert_true(setequal(
  unique(aa_summary_dynamic$accession),
  accession_gedi$accession
))
assert_true(setequal(phy$tip.label, accession_gedi$accession))

# make divergence matrix and index in same order as acc
Dmat = cophenetic(phy)
i = as.numeric(as.factor(colnames(Dmat)))
Dmat1 = Dmat[i, i] / max(Dmat)

# Compose data
stan_data_steadystate = aa_summary_steadystate |>
  rename(acc = accession) %T>%
  write_rds("objects/stan_data_df_sty.rds") |>
  compose_data()

stan_data_dynamic = aa_summary_dynamic |>
  rename(acc = accession) %T>%
  write_rds("objects/stan_data_df_dyn.rds") |>
  compose_data()

# add GEDI data
stan_accession_gedi = compose_data(accession_gedi)
stan_accession_gedi$n = NULL
stan_accession_gedi$accession = NULL
stan_accession_gedi$n_accession = NULL
stan_data = c(stan_data_steadystate, stan_accession_gedi)

# add phylogeny
stan_data$Dmat = Dmat1

assert_false(any(duplicated(names(stan_data))))

write_rds(stan_data, "data/stan_data.rds")
