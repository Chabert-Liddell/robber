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
#' @param method A string, the method used to compute the robustness by block.
#' One of:
#'   * \code{"exact"} (the default), average on all possible networks
#'   * \code{"no_empty_block"} restricts the networks space to the ones with no
#'   empty block of primary species
#'   * \code{"expected"} restricts the networks space to the ones that have
#'   approximately the same number of species per block than the expected one.
#' Do not do anything for \code{ext_seq = "uniform"}.
#' @param approx_threshold A numeric, the maximum size of the possible block
#' partition allowed for exact robustness by block calculation. Higher threshold
#' gives more precise results at the cost of computation times and possibly memory
#' problem. Do not do anything for \code{ext_seq = "uniform"}. When the threshold
#' is exceeded the robustness is computed by a Monte Carlo approximation with
#' \code{approx_threshold} number of samples.
#' @param net A network, if given, the function will fit a LBM to obtain the
#' parameters of the network and then compute the robustness.
#' @param ... Option to be passed to get_\code{lbm_param} function.
#' @return A list and a robber type object:
#' * \code{$fun} the robustness function, a vector of size \code{nr +1}
#' * \code{$auc} the area under the curve of the robustness function
#' * \code{$block} a vector of size \code{length(pi)}, the block ordering for
#'   primary extinctions sequence by blocks. \code{NULL} if \code{ext_seq = "uniform"}.
#' * \code{$model}, \code{$method}, \code{$ext_seq}, \code{$param}.

#' @export
#'
#' @examples
#' con <- matrix(c(.5,.3,.3,.1), 2, 2)
#' pi  <- c(.25,.75)
#' rho <- c(1/3, 2/3)
#' nr <- 50
#' nc <- 30
#' my_rob <- robustness_lbm(con, pi, rho, nr, nc, ext_seq = "natural")
#' my_rob$fun
#' my_rob$auc
#' # A easier alternative way, if you don't know the parameters of the network:
#' data(hostparasite)
#' (robustness_lbm(net = hostparasite, ncores = 1L))
robustness_lbm <- function(con = NULL, pi = NULL, rho = NULL,
                           nr = NULL, nc = NULL,
                           ext_seq = "uniform", method = "exact",
                           approx_threshold = 1e4, net = NULL, ...) {
  if(! is.null(net)) {
    param <- get_lbm_param(net, ... = ...)
    con <- param$con
    pi <- param$pi
    rho <- param$rho
    nr <- param$nr
    nc <- param$nc
  }
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
        result <-
          list(
            model = "lbm",
             ext_seq = ext_seq,
             fun = 1-result,
             auc = sum(1-result)/nr,
             block = NULL,
             method = method,
             param = list(con = con, pi = pi, rho = rho, nr = nr, nc = nc),
             sd = sqrt(do.call(var_fun_unif_lbm,
                          list(con = con, pi = pi, rho = rho, nr = nr, nc = nc))))
        class(result) <- "robber"
        return(result)
      },
      "increasing" = {
        k_seq <- order(con %*% rho, decreasing = TRUE)
        pi <- pi[k_seq]
        con <- con[k_seq, , drop=FALSE]
        result <- rob_block_lbm(con, pi, rho, nr, nc, method, approx_threshold)
        result$block <- k_seq
        result$method <- method
        result$param <- list(con = con, pi = pi, rho = rho,
                                     nr = nr, nc = nc)
        class(result) <- "robber"
        return( result )
      },
      "decreasing" = {
        k_seq <- order(con %*% rho, decreasing = FALSE)
        pi <- pi[k_seq]
        con <- con[k_seq, , drop=FALSE]
        result <- rob_block_lbm(con, pi, rho, nr, nc, method, approx_threshold)
        result$block <- k_seq
        result$method <- method
        result$param <- list(con = con, pi = pi, rho = rho,
                             nr = nr, nc = nc)
        class(result) <- "robber"
        return( result )
      },
      "natural" = {
        k_seq <- seq(length(pi))
        result <- rob_block_lbm(con, pi, rho, nr, nc, method, approx_threshold)
        result$block <- k_seq
        result$method <- method
        result$param <- list(con = con, pi = pi, rho = rho,
                             nr = nr, nc = nc)
        class(result) <- "robber"
        return( result )
      },
      # "blocks" = {
      #   K <- length(pi)
      #   perm <- gtools::permutations(n = K, r = K)
      #   res <- vector("list", factorial(K))
      #   for (k in seq(factorial(K))) {
      #     k_seq <- perm[k,]
      #     pi_seq <- pi[k_seq]
      #     con_seq <- con[k_seq, , drop=FALSE]
      #     result <- rob_block_lbm(con_seq, pi_seq, rho, nr, nc, method, approx_threshold)
      #     result$block <- k_seq
      #     res[[k]] <- result
      #   }
      #   return(res)
      # }
    )
}

#' @title Fast computation of the LBM robustness statistic for uniform extinctions
#'
#' @description This function provides much faster computation than \code{robustness_lbm} for
#' the case where \code{ext_seq = "uniform"}, when just the AUC statistic is
#' needed and the robustness function is not needed. This is particularly
#' useful if \code{nr} gets large.
#'
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


rob_block_lbm <- function( con, pi, rho, nr, nc, method = "exact", approx_threshold = 1e4) {
  if (length(pi) == 1) return(robustness_lbm(con, pi, rho, nr, nc, ext_seq = "uniform"))
  if (choose(nr + nrow(con) -1, nr) > approx_threshold) {
    approx <- TRUE
    message(paste0("Calculation by Monte Carlo because ",
                   choose(nr + nrow(con) -1, nr),
                   " is greater than approx threshold = ",
                   approx_threshold,
                   ". \n If you don't want this change, set approx_threshold to a higher value."))
  } else {
    approx <- FALSE
  }
  rob <- rep(1,nr+1)
  pi <- rev(pi)
  con <- con[rev(seq(nrow(con))),,drop = FALSE]
  if (method == "expected") {
    X <- nr*pi
    cumX <- cumsum(X)
    for (m in seq(0, nr)) {
      rob[m+1] <- rob[m+1] -
        as.vector(rho %*% exp(colSums(pmax(0,pmin(X, cumX - m))*log(1-con))))
    }
    return(list(fun = rob,
                auc = sum(rob)/nr))
  } else {
    if (! approx) {
      X <- partitions::compositions(nr, length(pi))
    } else {
      X <- stats::rmultinom(n = approx_threshold, size = nr, prob = pi)
    }
  }
  if (method == "no_empty_block") {
    X <- X[,colSums(X==0) == 0]
  }
  if (approx) {
    pmult <- 1/ncol(X)
  } else {
    pmult <- apply(X, 2, function(x) stats::dmultinom(x, prob = pi))
  }
  if (method == "no_empty_block" & ! approx) {
    pmult <- pmult/sum(pmult)
  }
  cumX <- apply(X, 2, cumsum)
  for (q in seq_along(rho)) {
    for (m in seq(0, nr)) {
      rob[m+1] <- rob[m+1] - rho[q] *
        sum(pmult * exp ( rowSums( log(
          vapply(seq_along(pi),
                 function(i) (1-con[i,q])**pmax(0,pmin(X[i,], cumX[i,] - m)),
                 FUN.VALUE = rep(0, ncol(X))))
          )))
    }
  }
  return(list(model = "lbm",
              fun = rob,
              auc = sum(rob)/nr))
}

