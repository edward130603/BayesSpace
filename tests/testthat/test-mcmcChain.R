test_that("Index labels are correct", {
  name <- "x"
  names <- .make_index_names(name, 2)
  expect_true(names[[2]] == "x[2]")
  
  names <- .make_index_names(name, 2, 3) 
  expect_true(names[[2]] == "x[1,2]")
  
  names <- .make_index_names(name, 2, 3, 0) 
  expect_true(names[[2]] == "x[2,1]")
})

test_that("Reading and writing works as expected", {
  out <- list()
  n_iter <- 10
  
  set.seed(149)
  z <- matrix(rep(rnorm(5), n_iter), ncol=5)
  mu <- matrix(rep(rnorm(20), n_iter), ncol=20)
  
  colnames(z) <- .make_index_names("z", 5)
  colnames(mu) <- .make_index_names("mu", 5, 4)
  
  out$z <- z
  out$mu <- mu
  
  h5.fname <- .write_chain(out)
  chain <- .read_chain(h5.fname, params=c("z", "mu"))
  
  expect_true(ncol(chain) == ncol(z) + ncol(mu))
  expect_true(nrow(chain) == n_iter)
  expect_true(all(chain[, 1:5] == z))
  expect_true(all(chain[, 6:25] == mu))
  expect_true(colnames(chain)[2] == "z[2]")
  expect_true(colnames(chain)[7] == "mu[1,2]")
  expect_true(colnames(chain)[10] == "mu[2,1]")
})

test_that("List of matrices is flattened", {
  xs <- list()
  n_iter <- 10
  n_row <- 3
  n_col <- 5
  
  set.seed(149)
  xs[[1]] <- matrix(rnorm(n_row * n_col), ncol=n_col)
  xs <- rep(xs, n_iter)
  
  x <- .flatten_matrix_list(xs, "x", n_row, n_col)
  
  expect_equal(nrow(x), n_iter)
  expect_equal(ncol(x), n_row * n_col)
  expect_equal(as.numeric(x[1, 7]), xs[[1]][2, 2])
  expect_equal(colnames(x)[7], "x[2,2]")
})

test_that("cleaning and saving works", {
  sce <- exampleSCE()
  n_rep <- 100
  n_PCs <- ncol(reducedDim(sce, "PCA"))
  n_spots <- ncol(sce)
  q <- 4
  
  sce <- spatialCluster(sce, 4, model="normal", nrep=n_rep, save.chain=T)
  chain <- mcmcChain(sce)
  
  # lambda + mu + ploglik + z
  expect_equal(ncol(chain), n_PCs * n_PCs + q * n_PCs + 1 + n_spots)
  expect_equal(nrow(chain), n_rep)
  
  sce <- spatialCluster(sce, 4, model="t", nrep=n_rep, save.chain=T)
  chain <- mcmcChain(sce)
  
  # lambda + mu + ploglik + weights + z
  expect_equal(ncol(chain), n_PCs * n_PCs + q * n_PCs + 1 + n_spots + n_spots)
  expect_equal(nrow(chain), n_rep)
})