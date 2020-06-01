context("C++")
test_that("Catch unit tests pass", {
    set.seed(149)
    expect_cpp_tests_pass("BayesSpace")
})
