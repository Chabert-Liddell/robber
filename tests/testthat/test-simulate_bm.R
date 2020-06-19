con <- matrix(c(.5,.3,.3,.1), 2, 2)
pi  <- c(.25,.75)
rho <- c(1/3, 2/3)
nr <- 50
nc <- 30

test_that("test simulation lbm gnp", {
  expect_error(simulate_lbm(con, c(.3, .8), rho, nr, nc))
  expect_error(simulate_lbm(con, pi, c(-2, 2.5), nr, nc))
  expect_error(simulate_lbm(3*con, pi, rho, nr, nc))
  expect_error(simulate_lbm(con, 1, rho, nr, nc))
  expect_error(simulate_lbm(con, pi, 1, nr, nc))
  expect_error(simulate_lbm(matrix(c(.5,.3,.3,.1), 4, 2), pi, rho, nr, nc))
  expect_error(simulate_lbm(matrix(c(.5,.3,.3,.1), 2, 4), pi, rho, nr, nc))
  expect_equal(dim(simulate_lbm(matrix(c(.5,.3,.3,.1), 2, 2), pi, rho, nr, nc)$A),
               c(nr, nc))
  expect_equal(length(simulate_lbm(matrix(c(.5,.3,.3,.1), 2, 2), pi, rho, nr, nc)$Z),
               nr)
  expect_equal(length(simulate_lbm(matrix(c(.5,.3,.3,.1), 2, 2), pi, rho, nr, nc)$W),
               nc)
})


con <- matrix(c(50, 120, 60, 80), 2, 2)
pi  <- c(10, 40)
rho <- c(10, 20)
nr <- 50
nc <- 30
simulate_lbm(con, pi, rho, nr, nc, method = "gnm")
test_that("test simulation lbm gnm", {
  expect_error(simulate_lbm(con, c(11, 40), rho, nr, nc, method = "gnm"))
  expect_error(simulate_lbm(con, pi, c(10, 18), nr, nc, method = "gnm"))
  expect_equal(sum(simulate_lbm(con,
                                pi, rho, nr, nc, method = "gnm")$A),
               sum(con))
  expect_equal(as.vector(table(simulate_lbm(con,
                                  pi, rho, nr, nc, method = "gnm")$Z)),
               pi)
  expect_equal(as.vector(table(simulate_lbm(con,
                                  pi, rho, nr, nc, method = "gnm")$W)),
               rho)
})
