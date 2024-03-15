# Prepare data for analysis
source("r/header.R")

read_rds("data/thinned_rh_curves.rds") |>
  mutate(
    RH = H2O_s / li6800_svp(Tair, Pa),
    leaf_type = case_when(
      curve_type == "1-sided RH" ~ "pseudohypo",
      curve_type == "2-sided RH" ~ "amphi"
    )
  ) |>
  separate_wider_delim(acc_id,
                       "-",
                       names = c("acc", "id"),
                       cols_remove = FALSE) |>
  # filter out Hawaiian, cultivated Solanum, and LA3909 (not enough reps)
  dplyr::filter(str_detect(acc, "^LA[0-9]{4}"), acc != "LA3475", acc != "LA3909") |>
  # make unique ID for each leaf_type within acc_id
  unite("lightintensity_x_acc_id", acc_id, light_intensity, remove = FALSE) |>
  # make unique ID for each curve
  unite("curve", acc_id, leaf_type, light_intensity, remove = FALSE) |>
  mutate(
    log_gsw = log(gsw),
    scaled_log_gsw = (log_gsw - mean(log_gsw)) / sd(log_gsw)
  ) |>
  write_rds("data/prepared_rh_curves.rds")

# phylogeny
tr = read.tree("data/Pease_etal_TomatoPhylo_RAxMLConcatTree_alltaxa_FigS2A.nwk") |>
  drop.tip(c("LA3475", "SL2.50"))

# Use LA3909 position for LA1044
tr$tip.label[tr$tip.label == "LA3909"] = "LA1044"

# graft LA0750 as sister to LA 0716
el = tr$edge.length[tr$edge[,2] == which(tr$tip.label == "LA0716")] / 2
tr1 = AddTip(tr, "LA0716", "LA0750", el)

write_rds(tr1, "data/phylogeny.rds")
