# Prepare actual data for Stan
source("r/header.R")

# args = commandArgs(trailingOnly = TRUE)
args = list(1)

phy = read_rds("data/phylogeny.rds")

aa_post = read_rds("objects/aa_post.rds") |>
  mutate(scaled_aa = (aa - mean(aa)) / sd(aa))
n_draw = max(aa_post$.draw)


stan_data = seq_len(n_draw) |>
  map(\(i) {
    df_data = aa_post |>
      filter(.draw == i) |>
      select(acc, acc_id, light_intensity, light_treatment, scaled_aa)
    
    accession_climate = read_rds("data/accession-climate.rds") |>
      dplyr::select(acc1 = accession, ppfd_mol_m2) |>
      dplyr::filter(acc1 %in% unique(df_data$acc)) |>
      mutate(scaled_ppfd_mol_m2 = (ppfd_mol_m2 - mean(ppfd_mol_m2)) / sd(ppfd_mol_m2))
    
    assert_true(setequal(unique(df_data$acc), accession_climate$acc1))
    assert_true(setequal(phy$tip.label, accession_climate$acc1))
    
    # make divergence matrix and index in same order as acc
    Dmat = cophenetic(phy)
    i = as.numeric(as.factor(colnames(Dmat)))
    Dmat1 = Dmat[i, i] / max(Dmat)
    
    # Compose data
    stan_data1 = df_data |>
      compose_data()
    
    # add SPLASH data
    stan_accession_climate = compose_data(accession_climate)
    stan_accession_climate$n = NULL
    stan_data1 = c(stan_data1, stan_accession_climate)
    
    # add phylogeny
    stan_data1$Dmat = Dmat1
    
    assert_false(any(duplicated(names(stan_data1))))
    
    stan_data1
    
  })

write_rds(stan_data, "data/stan_data.rds")
