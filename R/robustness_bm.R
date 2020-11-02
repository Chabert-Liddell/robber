#' @title Robustness for LBM
#'
#' @description Compute the robustness for a given set of Latent Block Model
#' parameters.
#'
#' @inheritParams auc_robustness_lbm
#' @param ext_seq A string, the rule for the primary extinctions sequences,
#' one of:
#' * \code{"uniform"}, the default for uniform extinctions sequences,
#' * \code{"decreasing"} and \code{"increasing"} for primary extinctions sequences
#'  by increasing and decreasing row blocks connection probability,
#' * \code{"natural"} primary extinctions sequences with the block ordering given
#'   in the function parameter.
#' * \code{"blocks"} primary extinctions sequences for all blocks permutation.
#' @return A list of size 3 except for \code{ext_seq = "blocks"}:
#' * \code{$fun} the robustness function, a vector of size \code{nr +1}
#' * \code{$auc} the area under the curve of the robustness function
#' * \code{$block} a vector of size \code{length(pi)}, the block ordering for
#'   primary extinctions sequence by blocks. \code{NULL} if \code{ext_seq = "uniform"}.
#'
#'   If \code{ext_seq = "blocks"}, then a list of length QR!,
#'   where QR is the length of the parameter pi. Each element of the list is a
#'   list of size 3 as above.
#' @export
#'
#' @examples
#' con <- matrix(c(.5,.3,.3,.1), 2, 2)
#' pi  <- c(.25,.75)
#' rho <- c(1/3, 2/3)
#' nr <- 50
#' nc <- 30
#' my_rob <- robustness_lbm(con, pi, rho, nr, nc)
#' my_rob$fun
#' my_rob$auc
#'
#' #' con <- matrix(c(.5,.3,.3,.1), 2, 2)
#' pi  <- c(.25,.75)
#' rho <- c(1/3, 2/3)
#' nr <- 50
#' nc <- 30
#' my_rob <- robustness_lbm(con, pi, rho, nr, nc, ext_seq = "blocks")
#' length(my_rob)
#' my_rob[[1]]$fun
#' my_rob[[1]]$auc
#' my_rob[[2]]$fun
#' my_rob[[2]]$auc
robustness_lbm <- function( con, pi, rho, nr, nc, ext_seq = "uniform") {
  if (any(con >1 | con <0)) {
    stop("Connectivity parameters should be between 0 and 1.")
  }
  if(nrow(con) != length(pi)) stop("Invalid number of row blocks.")
  if(ncol(con) != length(rho)) stop("Invalid number of col blocks.")
  if(! all.equal(sum(pi), 1)) stop("Vector pi must sum to 1.")
  if(! all.equal(sum(rho), 1)) stop("Vector rho must sum to 1.")
  if (! ext_seq %in% c("uniform", "increasing", "decreasing", "natural", "blocks")) {
    stop(paste0("Unknown extinction sequence: ", ext_seq, "."))
  }
  result <- rep(0,nr+1)
  k_seq <-
    switch (
      ext_seq,
      "uniform" = {
        delta <- pi%*%(1-con)
        for (q in seq_along(rho)) {
          result <- result + exp(seq.int(nr, 0) * log(delta[q]) ) * rho[q]
        }
        return(list(fun = 1-result,
                    auc = sum(1-result)/nr,
                    block = NULL))
      },
      "increasing" = {
        k_seq <- order(con %*% rho, decreasing = TRUE)
        pi <- pi[k_seq]
        con <- con[k_seq, , drop=FALSE]
        result <- rob_block_lbm(con, pi, rho, nr, nc)
        result$block <- k_seq
        return( result )
      },
      "decreasing" = {
        k_seq <- order(con %*% rho, decreasing = FALSE)
        pi <- pi[k_seq]
        con <- con[k_seq, , drop=FALSE]
        result <- rob_block_lbm(con, pi, rho, nr, nc)
        result$block <- k_seq
        return( result )
      },
      "natural" = {
        k_seq <- seq(length(pi))
        result <- rob_block_lbm(con, pi, rho, nr, nc)
        result$block <- k_seq
        return( result )
      },
      "blocks" = {
        K <- length(pi)
        perm <- gtools::permutations(n = K, r = K)
        res <- vector("list", factorial(K))
        for (k in seq(factorial(K))) {
          k_seq <- perm[k,]
          pi_seq <- pi[k_seq]
          con_seq <- con[k_seq, , drop=FALSE]
          result <- rob_block_lbm(con_seq, pi_seq, rho, nr, nc)
          result$block <- k_seq
          res[[k]] <- result
        }
        return(res)
      }
    )
}

#' @title Compute the AUC of the LBM robustness function
#'
#' This function provides much faster computation than \code{robustness_lbm} for
#' the case where \code{ext_seq = "uniform"}, when just the AUC statistic is
#' needed and the robustness function is not needed. This is particularly
#' useful if \code{nr} gets large.
#' @param con A matrix, the connectivity parameter
#' @param pi A vector of length \code{nrow(con)}, the proportion of row blocks
#' @param rho A vector of length \code{ncol(con)}, the proportion of
#' column blocks
#' @param nr An integer, the number of row (primary) species
#' @param nc An integer, the number of column (secondary) species
#'
#' @return The AUC of the LBM Robustness function for uniform primary
#' extinctions sequences.
#' @export
#'
#' @examples
#' con <- matrix(c(.5,.3,.3,.1), 2, 2)
#' pi  <- c(.25,.75)
#' rho <- c(1/3, 2/3)
#' nr <- 50
#' nc <- 30
#' auc_robustness_lbm(con, pi, rho, nr, nc)
auc_robustness_lbm <- function(con, pi, rho, nr, nc) {
  if (any(con >1 | con <0)) {
    stop("Connectivity parameters should be between 0 and 1.")
  }
  if(nrow(con) != length(pi)) stop("Invalid number of row blocks.")
  if(ncol(con) != length(rho)) stop("Invalid number of col blocks.")
  if(! all.equal(sum(pi), 1)) stop("Vector pi must sum to 1.")
  if(! all.equal(sum(rho), 1)) stop("Vector rho must sum to 1.")
  as.vector(1 - (1/nr)*(((1-pi %*% con) - (1-pi %*% con)**(nr+1))/
                          (pi%*%con))%*%rho)
}


rob_block_lbm <- function( con, pi, rho, nr, nc) {
  if (length(pi) == 1) return(robustness_lbm(con, pi, rho, nr, nc, ext_seq = "uniform"))
  X <- partitions::compositions(nr, length(pi))
  pi <- rev(pi)
  con <- con[rev(seq(nrow(con))),,drop = FALSE]
  pmult <- apply(X, 2, function(x) dmultinom(x, prob = pi))
  rob <- rep(1,nr+1)
  storage.mode(X) <- "numeric"
  cumX <- Rfast::colCumSums(X)
  for (q in seq_along(rho)) {
    for (m in seq(0, nr)) {
      rob[m+1] <- rob[m+1] - rho[q] * sum(pmult *
                                            Rfast::rowprods(
                                              vapply(seq_along(pi),
                                                     function(i) (1-con[i,q])**pmax(0,pmin(X[i,], cumX[i,] - m)),
                                                     FUN.VALUE = rep(0, ncol(X)))))
    }
  }
  # rob
  # result <- rep(0,nr+1)
  # K <- length(pi)
  # order_prob <- matrix(0, nr+1, K)
  # order_prob[nr+1, K] <- 0
  # for (k in seq_len(K-1)) {
  #   for (m in seq.int(0, nr-1)) {
  #     order_prob[m+1, k] <-
  #       sum(stats::dbinom(x = seq.int(0, nr-m-1), size = nr,
  #                         prob = 1-sum(pi[1:k])) -
  #             stats::dbinom(x = seq.int(0, nr-m-1), size = nr,
  #                           prob = 1 - sum(pi[1:k]) + pi[k]))
  #   }
  # }
  # order_prob[,K] <- 1 - rowSums(order_prob)
  # tmp_rob <- matrix(0, nr+1, length(rho))
  # for(m in seq.int(0, nr)) {
  #   for(k in seq_len(K)) {
  #     tmp_rob[m+1,] <- tmp_rob[m+1,] +
  #       (( (pi[1:k]/sum(pi[1:k])) %*% (1-con[1:k, , drop=FALSE]))**(m) *
  #          order_prob[m+1, k])
  #   }
  # }
  # result[1:(nr+1)] <- as.vector(tmp_rob %*% rho)[(nr+1):1]
  return(list(fun = rob,
              auc = sum(rob)/nr))
}
