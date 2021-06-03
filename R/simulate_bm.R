
#' @title Simulate a bipartite interaction network using block model
#'
#' @param con A matrix, the connectivity between blocks. If
#' \code{method = "gnp"} then each entry is the probability of interactions
#' given 2 blocks. If \code{method = "gnm"} then each entry is the number of
#' interactions between 2 blocks.
#' @param pi A vector of the same length as \code{nrow(con)},
#' the block mixture for the row species. If \code{method = "gnp"}, then
#' \code{pi} is a probability vector, if \code{method = "gnm"}, then
#' \code{pi} is the number of species per block, must sum to \code{nr}.
#' @param rho A vector of the same length as \code{ncol(con)},
#' the block mixture for the column species. If \code{method = "gnm"}, then
#' \code{rho} is a probability vector, if \code{method = "gnm"}, then
#' \code{rho} is the number of species per block, must sum to \code{nc}.
#' @param nr The number of row Species
#' @param nc The number of column Species
#' @param method One of \code{"gnp"} (the default) where the blocks size and
#' number of interactions is random and \code{"gnm"} where the blocks size and
#' number of interactions are fixed.
#'
#' @return A list of 3 elements:
#' * \code{A} The incident matrix of size \code{nr}x\code{nc}
#' * \code{Z} A vector of length \code{nr}, the block assignment of the row
#' species
#' * \code{W} A vector of length \code{nc}, the block assignment of the column
#' species
#' @export
#'
#' @examples
#' # For a random number of interactions and blocks sizes
#' con <- matrix(c(.5,.3,.3,.1), 2, 2)
#' pi  <- c(.25,.75)
#' rho <- c(1/3, 2/3)
#' nr <- 50
#' nc <- 30
#' simulate_lbm(con, pi, rho, nr, nc, method = "gnp")
#'
#' # For a fixed number of ineractions and blocks sizes
#' con <- matrix(c(50, 120, 60, 80), 2, 2)
#' pi  <- c(10, 40)
#' rho <- c(10, 20)
#' nr <- 50
#' nc <- 30
#' simulate_lbm(con, pi, rho, nr, nc, method = "gnm")
simulate_lbm <- function(con, pi, rho, nr, nc, method = "gnp") {
  # if (!requireNamespace("pkg", quietly = TRUE)) {
  #   stop("Package \"igraph\" needed for this function to work. Please install it.",
  #        call. = FALSE)
  # }
  if(nrow(con) != length(pi)) stop("Invalid number of row blocks.")
  if(ncol(con) != length(rho)) stop("Invalid number of col blocks.")
  A <- matrix(0, nr, nc)
  qr <- length(pi)
  qc <- length(rho)
  if (method == "gnp") {
    if(! all.equal(sum(pi), 1)) stop("Vector pi must sum to 1.")
    if(! all.equal(sum(rho), 1)) stop("Vector rho must sum to 1.")
    if (any(con >1 | con <0)) {
      stop("Connectivity parameters should be between 0 and 1.")
    }
    W <- sample(x = seq(qc), size = nc, replace = TRUE, prob = rho)
    Z <- sample(x = seq(qr), size = nr, replace = TRUE, prob = pi)
    for (k in seq(qr)) {
      for (q in seq(qc)) {
        A[Z == k, W == q] <- igraph::as_incidence_matrix(
          igraph::sample_bipartite(sum(Z == k), sum(W == q),
                                   type = "gnp", p = con[k,q]))
      }
    }
    return(list(A = A, W = W, Z = Z))
  }
  if (method == "gnm") {
    if (! all.equal(sum(pi),  nr)) stop("Must have sum(pi) = nr.")
    if (! all.equal(sum(rho), nc)) stop("Must have sum(rho) = nc.")
    W <-  rep(seq(qc), times = rho)
    Z <-  rep(seq(qr), times = pi)
    for (k in seq(qr)) {
      for (q in seq(qc)) {
        A[Z == k, W == q] <- igraph::as_incidence_matrix(
          igraph::sample_bipartite(pi[k], rho[q],
                                   type = "gnm", m = con[k,q]))
      }
    }
    return(list(A = A, W = W, Z = Z))
  }
}
