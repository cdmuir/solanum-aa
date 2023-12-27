# 60SITEOT
# The upshot is that the way Menura does beta diffusion seems wrong or hard to interpret, maybe because of transformation?
# I did figure out how to use sde package to simulate beta diffusion
# But does not seem like way forward. Rather, model SD on each surface as OU process

library(ape)
library(menura)
library(sde)

# Number of tips on the tree
ntips <- 128

# Iterations set to 50 to limit run time.  Most simulations require 200000 or more
iters <- 50

# Create a random tree
tr <-  compute.brlen(stree(n=ntips, type="balanced"))
tr$edge.length = tr$edge.length * 10

pars = list(
  mu = 0.5,
  alpha = 1,
  epsilon = 1
)

f_beta = function(x, l, pars) {
  
  suppressMessages({
    d = parse(text = glue::glue("- {a} * (x - {u})", a = pars$alpha, u = pars$mu))
    s = parse(text = glue::glue("{e} * x * (1 - x)", e = pars$epsilon))
    sim = sde::sde.sim(T = l, X0 = x, drift = d, sigma = s)
  })
  sim[length(sim)]
  
}

f_beta(0.5, 10, pars)

t.tipdata <- rTraitCont(tr, f_beta, ancestor = FALSE, root.value = pars$mu, pars = pars)

hist(t.tipdata)

set.seed(1)
model.1 <- fit_model(tr = tr, tipdata = t.tipdata, rt_value = pars$mu, thin = 100,
                     iters = 200000, model = "Beta", alpha = NULL, mu = NULL, 
                     sigma = NULL, N = 240, init_method = "sim", 
                     update_method = "subtree")

# Look at the MCMC trace of parameters
summary(model.1)
model.1

## Use coda package to analyse the mcmc chain
require(coda)
plot(model.1$mcmctrace)
summary(model.1$mcmctrace)

# theta[1] = alpha
# theta[2] = mu
# theta[3] = sigma (or epsilon)
m2 = list(
  # Mparam <-  c("d", "s", "drift", "diffusion", "dx_diffusion")
  d = function(t, x, theta) {
    theta[1] * (theta[2] - x)
  },
  s = function(t, x, theta) {
    theta[3] * x * (1 - x)
  },
  s_x = function(t, x, theta) {
    0
  },
  drift = quote(alpha * (mu - x)),
  diffusion = quote(sigma * x * (1 - x)),
  dx_diffusion = quote(0)
)

# this didn't work to estimate mu because it generated negative values
model.2 <- fit_model(tr = tr, tipdata = t.tipdata, rt_value = pars$mu, thin = 1,
                     iters = 4000, model = m2, alpha = 1, mu = NULL, 
                     sigma = 1, N = 240, init_method = "sim", 
                     update_method = "subtree")

summary(model.2)
plot(model.2$mcmctrace)
