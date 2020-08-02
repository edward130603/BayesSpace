test_that("neighbor finding works", {
  xs <- 1:5
  ys <- 5:1
  positions <- cbind(xs, ys)
  
  df_j <- find_neighbors(positions, 1, "manhattan")
  expect_true(all(vapply(df_j, length, integer(1)) == 0))
  
  df_j <- find_neighbors(positions, 4, "manhattan")
  expect_true(all(df_j[[1]] == c(2, 3) - 1))
  expect_true(all(df_j[[3]] == c(1, 2, 4, 5) - 1))
})
