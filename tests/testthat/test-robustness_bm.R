con <- matrix(c(.5,.3,.3,.1), 2, 2)
pi  <- c(.25,.75)
rho <- c(1/3, 2/3)
nr <- 50
nc <- 30
test_that("robustness lbm error", {
  expect_error(robustness_lbm(con, c(.3, .8), rho, nr, nc))
  expect_error(robustness_lbm(con, pi, c(-2, 2.5), nr, nc))
  expect_error(robustness_lbm(3*con, pi, rho, nr, nc))
  expect_error(robustness_lbm(con, 1, rho, nr, nc))
  expect_error(robustness_lbm(con, pi, 1, nr, nc))
  expect_error(robustness_lbm(matrix(c(.5,.3,.3,.1), 4, 2), pi, rho, nr, nc))
  expect_error(robustness_lbm(matrix(c(.5,.3,.3,.1), 2, 4), pi, rho, nr, nc))
})


test_that("robustness lbm result", {
  expect_equal(length(robustness_lbm(con, pi, rho, nr, nc)$fun), nr+1)
  expect_equal(
    length(robustness_lbm(con, pi, rho, nr, nc, ext_seq = "increasing")$fun),
    nr+1)
  expect_equal(
    length(robustness_lbm(con, pi, rho, nr, nc, ext_seq = "decreasing")$fun),
    nr+1)
  expect_equal(robustness_lbm(con, pi, rho, nr, nc)$auc,
               auc_robustness_lbm(con, pi, rho, nr, nc))
  expect_lte(robustness_lbm(con, pi, rho, nr, nc)$auc,
             robustness_lbm(con, pi, rho, nr, nc, ext_seq = "increasing")$auc)
  expect_gte(robustness_lbm(con, pi, rho, nr, nc)$auc,
             robustness_lbm(con, pi, rho, nr, nc, ext_seq = "decreasing")$auc)
  expect_lte(robustness_lbm(con, pi, rho, nr, nc)$auc,
             robustness_lbm(con, pi, rho, nr+1, nc)$auc)
  expect_lte(robustness_lbm(con, pi, rho, nr, nc)$auc,
             robustness_lbm(1.1*con, pi, rho, nr, nc)$auc)
})
