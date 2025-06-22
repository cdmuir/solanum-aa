# Prepare actual data for Stan
source("r/header.R")

plant_info = read_rds("data/plant-info.rds") |>
  mutate(llma = log(lma_gm2)) |>
  dplyr::select(accession, acc_id, light_treatment, llma) 

trimmed_amphi_first = read_rds("data/trimmed_amphi_first.rds") |>
  summarize(
    amphi_first = unique(amphi_first),
    .by = c("acc_id")
  )

phy = read_rds("data/phylogeny.rds")

aa_summary = read_rds("objects/aa_summary.rds") |>
  dplyr::select(acc_id, light_intensity, aa = median) |>
  mutate(scaled_aa = (aa - mean(aa)) / sd(aa)) |>
  left_join(plant_info, by = join_by(acc_id)) |>
  left_join(trimmed_amphi_first, by = join_by(acc_id)) 

accession_gedi = read_rds("data/accession-gedi.rds") |>
  dplyr::select(accession, pai) |>
  dplyr::filter(accession %in% unique(aa_summary$accession)) |>
  mutate(scaled_pai = (pai - mean(pai)) / sd(pai))

assert_true(setequal(unique(aa_summary$accession), accession_gedi$accession))
assert_true(setequal(phy$tip.label, accession_gedi$accession))
    
# make divergence matrix and index in same order as acc
Dmat = cophenetic(phy)
i = as.numeric(as.factor(colnames(Dmat)))
Dmat1 = Dmat[i, i] / max(Dmat)

# Compose data
stan_data = aa_summary |>
  rename(acc = accession) %T>%
  write_rds("objects/stan_data_df.rds") |>
  compose_data()

# add GEDI data
stan_accession_gedi = compose_data(accession_gedi)
stan_accession_gedi$n = NULL
stan_accession_gedi$accession = NULL
stan_accession_gedi$n_accession = NULL
stan_data = c(stan_data, stan_accession_gedi)

# add phylogeny
stan_data$Dmat = Dmat1

assert_false(any(duplicated(names(stan_data))))

write_rds(stan_data, "data/stan_data.rds")
