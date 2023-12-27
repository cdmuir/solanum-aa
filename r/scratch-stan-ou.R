# Notes from Statistical Rethinking v2 on OU model as Gaussian Process in Stan
library(ape)
library(rethinking)

data(Primates301)
data(Primates301_nex)

plot(ladderize(Primates301_nex), type = "fan", font = 1, no.margin = TRUE,
     label.offset = 1, cex = 0.5)

d = Primates301
d$name = as.character(d$name)
dstan = d[complete.cases(d$group_size, d$body, d$brain), ]
spp_obs = dstan$name

dat_list = list(
  N_spp = nrow(dstan),
  M = standardize(log(dstan$body)),
  B = standardize(log(dstan$brain)),
  G = standardize(log(dstan$group_size)),
  Imat = diag(nrow(dstan))
)

m14.9 = ulam(
  alist(
    B ~ multi_normal(mu, SIGMA),
    mu <- a + bM * M + bG * G,
    matrix[N_spp, N_spp]: SIGMA <- Imat * sigma_sq,
    a ~ normal(0, 1),
    c(bM, bG) ~ normal(0, 0.5),
    sigma_sq ~ exponential(1)
  ), data = dat_list, chains = 4, cores = 4, cmdstan = TRUE
) #13.8 seconds

precis(m14.9)

# Brownian motion
tree_trimmed = keep.tip(Primates301_nex, spp_obs)
Rbm = corBrownian(phy = tree_trimmed)
V = vcv(Rbm)
Dmat = cophenetic(tree_trimmed)
plot(Dmat, V, xlab = "phylogenetic distance", ylab = "covariance")

dat_list$V = V[spp_obs, spp_obs]
dat_list$R = dat_list$V / max(V)

m14.10 = ulam(
  alist(
    B ~ multi_normal(mu, SIGMA),
    mu <- a + bM * M + bG * G,
    matrix[N_spp, N_spp]: SIGMA <- R * sigma_sq,
    a ~ normal(0, 1),
    c(bM, bG) ~ normal(0, 0.5),
    sigma_sq ~ exponential(1)
  ), data = dat_list, chains = 4, cores = 4, cmdstan = TRUE
) #15.8 seconds

precis(m14.10)

# OU
dat_list$Dmat = Dmat[spp_obs, spp_obs] / max(Dmat)

m14.11 = ulam(
  alist(
    B ~ multi_normal(mu, SIGMA),
    mu <- a + bM * M + bG * G,
    matrix[N_spp, N_spp]: SIGMA <- cov_GPL1(Dmat, etasq, rhosq, 0.0),
    a ~ normal(0, 1),
    c(bM, bG) ~ normal(0, 0.5),
    etasq ~ half_normal(1, 0.25),
    rhosq ~ half_normal(3, 0.25)
  ), data = dat_list, chains = 4, cores = 4, cmdstan = TRUE
) #27.7 seconds

precis(m14.11)

stancode(m14.11)
