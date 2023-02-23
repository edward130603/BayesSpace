devtools::clean_dll()
devtools::load_all()

# set.seed(1)
# t <- spatialEnhance(
#   readRDS("__tools/data.RDS"),
#   q = 4,
#   platform = "ST",
#   d = 7,
#   model = "t",
#   gamma = 2,
#   jitter_prior = 0.3,
#   jitter_scale = 3.5,
#   nrep = 1000,
#   burn.in = 100,
#   save.chain = TRUE,
#   thread_num = 4
# )
#
# library(RcppClock)

check <- function(xs) {
  chain1 <- mcmcChain(xs[[1]])
  all(sapply(xs[-1], function(x) identical(chain1, mcmcChain(x))))
}

seed <- 1
file <- "__tools/clustered.rds"
q <- 15
d <- 16
platform <- "Visium"
gamma <- 4
jitter_scale <- 5
jitter_prior <- 0.3
nrep <- 1000
verbose <- TRUE

t <- microbenchmark::microbenchmark(
  t1 = {
    set.seed(seed)
    spatialEnhance(
      readRDS(file),
      q = q,
      platform = platform,
      d = d,
      model = "t",
      gamma = gamma,
      jitter_prior = jitter_prior,
      jitter_scale = jitter_scale,
      nrep = nrep,
      burn.in = 0.1 * nrep,
      save.chain = TRUE,
      thread_num = 1,
      verbose = verbose
    )
  },
  t2 = {
    set.seed(seed)
    spatialEnhance(
      readRDS(file),
      q = q,
      platform = platform,
      d = d,
      model = "t",
      gamma = gamma,
      jitter_prior = jitter_prior,
      jitter_scale = jitter_scale,
      nrep = nrep,
      burn.in = 0.1 * nrep,
      save.chain = TRUE,
      thread_num = 2,
      verbose = verbose
    )
  },
  t4 = {
    set.seed(seed)
    spatialEnhance(
      readRDS(file),
      q = q,
      platform = platform,
      d = d,
      model = "t",
      gamma = gamma,
      jitter_prior = jitter_prior,
      jitter_scale = jitter_scale,
      nrep = nrep,
      burn.in = 0.1 * nrep,
      save.chain = TRUE,
      thread_num = 4,
      verbose = verbose
    )
  },
  t6 = {
    set.seed(seed)
    spatialEnhance(
      readRDS(file),
      q = q,
      platform = platform,
      d = d,
      model = "t",
      gamma = gamma,
      jitter_prior = jitter_prior,
      jitter_scale = jitter_scale,
      nrep = nrep,
      burn.in = 0.1 * nrep,
      save.chain = TRUE,
      thread_num = 6,
      verbose = verbose
    )
  },
  t8 = {
    set.seed(seed)
    spatialEnhance(
      readRDS(file),
      q = q,
      platform = platform,
      d = d,
      model = "t",
      gamma = gamma,
      jitter_prior = jitter_prior,
      jitter_scale = jitter_scale,
      nrep = nrep,
      burn.in = 0.1 * nrep,
      save.chain = TRUE,
      thread_num = 8,
      verbose = verbose
    )
  },
  check = check,
  times = 1L
)
