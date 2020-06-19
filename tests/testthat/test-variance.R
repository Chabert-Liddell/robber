con <- matrix(c(.5,.3,.3,.1), 2, 2)
pi  <- c(.25,.75)
rho <- c(1/3, 2/3)
nr <- 50
nc <- 30
var_auc_unif_lbm(con, pi, rho, nr, nc)

test_that("variance of the lbm", {
  expect_message(var_auc_unif_lbm(con, pi, rho, 110, nc))
  expect_error(var_auc_unif_lbm(con, pi, c(-2, 2.5), nr, nc))
  expect_error(var_auc_unif_lbm(3*con, pi, rho, nr, nc))
  expect_error(var_auc_unif_lbm(con, 1, rho, nr, nc))
  expect_error(var_auc_unif_lbm(con, pi, 1, nr, nc))
  expect_error(var_auc_unif_lbm(matrix(c(.5,.3,.3,.1), 4, 2), pi, rho, nr, nc))
  expect_error(var_auc_unif_lbm(matrix(c(.5,.3,.3,.1), 2, 4), pi, rho, nr, nc))
})
