#' @title Compute the Conditional variance of the AUC of the LBM Robustness
#'
#'
#'
#' @inheritParams auc_robustness_lbm
#'
#' @return A numeric, the variance
#' @export
#' @importFrom purrr map2_dbl
#'
#' @examples
#' con <- matrix(c(.5,.3,.3,.1), 2, 2)
#' pi  <- c(.25,.75)
#' rho <- c(1/3, 2/3)
#' nr <- 50
#' nc <- 30
#' var_auc_unif_lbm(con, pi, rho, nr, nc)
var_auc_unif_lbm <- function(con, pi, rho, nr, nc) {
  if (nr > 100) message("Computation may take time for a network of this size.")
  if (any(con >1 | con <0)) {
    stop("Connectivity parameters should be between 0 and 1.")
  }
  if(nrow(con) != length(pi)) stop("Invalid number of row blocks.")
  if(ncol(con) != length(rho)) stop("Invalid number of col blocks.")
  if(! all.equal(sum(pi), 1)) stop("Vector pi must sum to 1.")
  if(! all.equal(sum(rho), 1)) stop("Vector rho must sum to 1.")
  K <- nrow(con)
  Q <- ncol(con)
  etaq <- pi %*% (1 - con)
  etaqq <- matrix(purrr::map2_dbl(
    rep(seq(Q), times = Q),
    rep(seq(Q), each = Q),
    function(.x, .y)  pi %*% ((1-con[,.x])*(1-con[,.y]))),
    nrow = Q, ncol = Q)
  result <- 0
  for (q in seq(Q)) {
    result <- result + (1/nr) * etaq[q]**(0:nr) * rho[q]
  }
  result <- - sum(result)**2
  for (m in seq(0, nr)) {
    cb <- t(outer(
      X = seq.int(0,nr),
      Y = seq.int(0,nr),
      FUN = function(nr, m, mp, l) {
        choose(m, m+mp-l) * choose(nr-m, l-m)/choose(nr, mp)
        }, nr = nr, m = m))
    l_mp <- outer(X = seq.int(0,nr), Y = seq.int(0,nr), `-`)
    l_mp[l_mp<0] <- 0
    for (q in seq(Q)) {
      result <- result + (1/(nr*nr*nc))  * sum(
        cb * etaq[q]**seq.int(0,nr) * rho[q])
      for (qp in seq(Q)) {
        result <- result + (1/(nr*nr))  * (1-1/nc) *
          sum( cb * rho[q] * rho[qp] * etaq[q]**(l_mp) *
                 etaq[qp]**(pmax(0, seq.int(0,nr)-m)) *
                 etaqq[q, qp]**(pmax(0, m - l_mp)))
      }
    }
  }
  return (as.vector(result))
}


#' @title Compute the Conditional variance of the LBM Robustness term by term
#'
#'
#'
#' @inheritParams auc_robustness_lbm
#'
#' @return A vector, the variance after m extinctions
#' @export
#' @importFrom purrr map2_dbl
#'
#' @examples
#' con <- matrix(c(.5,.3,.3,.1), 2, 2)
#' pi  <- c(.25,.75)
#' rho <- c(1/3, 2/3)
#' nr <- 50
#' nc <- 30
#' var_fun_unif_lbm(con, pi, rho, nr, nc)
var_fun_unif_lbm <- function(con, pi, rho, nr, nc) {
  if (nr > 100) message("Computation may take time for a network of this size.")
  if (any(con >1 | con <0)) {
    stop("Connectivity parameters should be between 0 and 1.")
  }
  if(nrow(con) != length(pi)) stop("Invalid number of row blocks.")
  if(ncol(con) != length(rho)) stop("Invalid number of col blocks.")
  if(! all.equal(sum(pi), 1)) stop("Vector pi must sum to 1.")
  if(! all.equal(sum(rho), 1)) stop("Vector rho must sum to 1.")
  K <- nrow(con)
  Q <- ncol(con)
  etaq <- pi %*% (1 - con)
  etaqq <- matrix(purrr::map2_dbl(
    rep(seq(Q), times = Q),
    rep(seq(Q), each = Q),
    function(.x, .y)  pi %*% ((1-con[,.x])*(1-con[,.y]))),
    nrow = Q, ncol = Q)
  result <- 0
  for (q in seq(Q)) {
    result <- result + etaq[q]**(0:nr) * rho[q]
  }
  result <- -(result**2)
  cb <- t(outer(
    X = seq.int(0,nr),
    Y = seq.int(0,nr),
    FUN = function(nr, mp, l) { # return l times mp matrix
      choose(mp, mp+mp-l) * choose(nr-mp, l-mp)/choose(nr, mp)
    }, nr = nr))
  for (m in seq(1, nr)) {
    # l_mp <- outer(X = seq.int(0,nr), Y = seq.int(0,nr), `-`)
    # l_mp[l_mp<0] <- 0
    for (q in seq(Q)) {
      result[m+1] <- result[m+1] + (1/nc)  * sum(
        cb[,m+1] * etaq[q]**seq.int(0,nr) * rho[q])
      for (qp in seq(Q)) {
        result[m+1] <- result[m+1] +  (1-1/nc) *
          sum( cb[,m+1] * rho[q] * rho[qp] * etaq[q]**(pmax(0, seq.int(0,nr)-m)) *
                 etaq[qp]**(pmax(0, seq.int(0,nr)-m)) *
                 etaqq[q, qp]**(pmax(0, 2*m - seq.int(0,nr))))
      }
    }
  }
  result[1] <- 0
  return (rev(result))
}
