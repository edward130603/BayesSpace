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
