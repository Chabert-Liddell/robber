
A <- matrix(rnorm(10), 5, 2)
test_that("warning in matrix argument", {
  expect_warning(robustness_emp(A))
})
test_that("error in sequence argument", {
  expect_error(robustness_emp(A, deg_seq = "block"))
})

A <- matrix(c(1, 0, 0,
              0, 1, 0,
              0, 0, 1,
              1, 1, 1), nrow = 4, ncol = 3, byrow = TRUE)
test_that("returning the right dimension", {
  expect_equal(length(robustness_emp(A, nb_iter = 1)$fun), 5)
  expect_equal(length(robustness_emp(A, ext_seq = "increasing")$fun), 5)
  expect_equal(length(robustness_emp(A, ext_seq = "decreasing")$fun), 5)
  expect_equal(length(robustness_emp(A, nb_iter = 1)$auc), 1)
})
test_that("robustness value check", {
  expect_equal(robustness_emp(A, ext_seq = "decreasing")$fun,
               c(1, 1, 2/3, 1/3, 0))
  expect_equal(robustness_emp(A, ext_seq = "increasing")$fun,
               c(1, 1, 1, 1, 0))
  expect_equal(robustness_emp(A, ext_seq = "decreasing")$auc,
               3/4)
  expect_equal(robustness_emp(A, ext_seq = "increasing")$auc,
               1)
  expect_lte(robustness_emp(A)$auc,
               1)
  expect_gte(robustness_emp(A)$auc,
             3/4)
})
