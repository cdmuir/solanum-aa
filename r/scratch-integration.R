source("r/header.R")

ll = list(n1 = 100,
          n2 = 50,
          d1 = 0, d2 = 0L)
ll$x1 = seq(0, 1, length.out = ll$n1)
ll$x2 = seq(0.5, 1.5, length.out = ll$n2)

ll$y1 = 11 + 1 * ll$x1 + rnorm(ll$n1, 0, 0.1)
ll$y2 = 10 + 1 * ll$x2 + rnorm(ll$n2, 0, 0.1)

assert_true(all(ll$y1 > 0))
assert_true(all(ll$y2 > 0))

m = cmdstan_model("stan/test-integration.stan", dir = "stan/bin")

fit = m$sample(
  data = ll,
  chains = 2L,
  parallel_chains = 2L,
  iter_warmup = 1e3,
  iter_sampling = 1e3,
  refresh = 2e1
)

fit$summary("A")
