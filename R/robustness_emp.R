#' @title Empirical Robustness of a Network
#'
#' @description Compute the robustness of an ecological network by averaging
#' over a great number of randomly generated primary extinctions sequences
#'
#' @param A A binary incident matrix
#' @param ext_seq A string, the rule for the primary extinctions sequences,
#' one of \code{"uniform"}, the default for uniform extinctions sequences,
#' \code{"decreasing"} and \code{"increasing"} for primary extinctions sequences
#'  by increasing and decreasing degree order sequence on the row species
#' @param nb_iter An integer, the number of generated sequences. Default to
#' \code{300}.
#' @param net_type A string, the type of network. For now, only "bipartite" is
#' available.
#' @param method A string used when \code{ext_seq = "decreasing"} or
#' \code{ext_seq = "increasing"}:
#'  * \code{"linear"} (default) set an extinction
#' probability for each row species that is linear in its degree. It is a
#' shortcut for \code{method = "power"} and \code{power = 1}.
#'  * \code{"ordered"} strictly follow the degree order, i.e. the most
#'   connected species will always goes last for \code{ext_seq = "increasing"}
#'   or first for \code{ext_seq = "decreasing"}.
#' @param power An integer (default to 1). Only used when
#' \code{method = "power"},
#' the power to which the degree are elevated to compute the primary extinctions
#'  sequences for \code{"increasing"} and \code{"decreasing"} ext_seq. For 1, it
#'  is equivalent to \code{method = "linear"}, for 0, it is the same as
#'  \code{ext_seq = "uniform"}. The largest the number the closest it will be to
#'   \code{method = "ordered"}.
#' @param block A vector of row species memberships for \code{method = "block"}
#' @return A list (robber object) with the following fields:
#' * \code{$model}, \code{$ext_seq}, \code{$method}, \code{power}
#' * \code{$fun} is the robustness function, a vector of size \code{(nrow(A) +1)}
#' * \code{$sd} Standard deviation of the secondary extinctions
#' * \code{$auc} the area under the curve of the robustness function
#' @export
#'
#' @examples
#' A <- matrix(c(1, 0, 0,
#'              0, 1, 0,
#'              0, 0, 1,
#'              1, 1, 1), nrow = 4, ncol = 3, byrow = TRUE)
#' my_rob <- robustness_emp(A, ext_seq = "increasing")
#' my_rob$fun
#' my_rob$auc
robustness_emp <- function(A, ext_seq = "uniform",
                           nb_iter = 300, net_type = "bipartite",
                           method = "ordered",
                           power = 1,
                           block = NULL) {
  if (any(! A %in% c(0, 1))) {
    warning(paste0("Only 0, 1 entry accepted.
                   Will transform all non zero entries into 1."))
    A[A!=0] <- 1
  }
  if (! ext_seq %in% c("uniform", "increasing", "decreasing")) {
    stop(paste0("Unknown extinction sequence: ", ext_seq, "."))
  }
  if (any(is.na(A))) warning("NA will be considered as 0.")
  nr <- nrow(A)
  nc <- ncol(A)
  if (ext_seq == "uniform" | method %in% c("linear", "power")) {
    prob_ext <- switch (ext_seq,
      "uniform" = rep(1, nr)/nr,
      "increasing" = pmax(1, rowSums(A, na.rm = TRUE))/sum(pmax(1, rowSums(A, na.rm = TRUE))),
      "decreasing" = pmax(1, rowSums(A, na.rm = TRUE))/sum(pmax(1, rowSums(A, na.rm = TRUE)))
    )
    sec_ext <- matrix(
      vapply(X = seq_len(nb_iter),
             FUN = function(i) {
               prim_ext <- sample(x = seq_len(nr), size = nr,
                                  prob = prob_ext**power,replace = FALSE)
               if(ext_seq == "increasing") prim_ext <- rev(prim_ext)
               c(sum(colSums(A, na.rm = TRUE) != 0),
                 vapply(X = seq.int(nrow(A)),
                        FUN = function(m) {
                          sum(colSums(A[prim_ext[-seq_len(m)], ,drop = FALSE],
                                      na.rm = TRUE) != 0)
                        },
                        FUN.VALUE = integer(length = 1)))
             },
             FUN.VALUE = integer(length = nr + 1)),
      nrow = nr+1, ncol = nb_iter)
  }
  if (ext_seq == "increasing" & method == "ordered") {
    if (length(unique(rownames(A))) != nrow(A)) {
      rownames(A) <- seq(nr)
    }
    ord <- rowSums(A)
    sec_ext <- matrix(
      vapply(X = seq_len(nb_iter),
             FUN = function(i) {
               prim_ext <- names(sort(ord[sample(x = seq(nr))]))
               c(sum(colSums(A, na.rm = TRUE) != 0),
                 vapply(X = seq.int(nrow(A)),
                        FUN = function(m) {
                          sum(colSums(A[prim_ext[-seq_len(m)], ,drop = FALSE],
                                      na.rm = TRUE) != 0)
                        },
                        FUN.VALUE = integer(length = 1)))
             },
             FUN.VALUE = integer(length = nr + 1)),
      nrow = nr+1, ncol = nb_iter)
  }
  if (ext_seq == "decreasing" & method == "ordered") {
    if (length(unique(rownames(A))) != nrow(A)) {
      rownames(A) <- seq(nr)
    }
    ord <- rowSums(A)
    sec_ext <- matrix(
      vapply(X = seq_len(nb_iter),
             FUN = function(i) {
               prim_ext <- names(sort(ord[sample(x = seq(nr))], decreasing = TRUE))
               c(sum(colSums(A, na.rm = TRUE) != 0),
                 vapply(X = seq.int(nrow(A)),
                        FUN = function(m) {
                          sum(colSums(A[prim_ext[-seq_len(m)], ,drop = FALSE],
                                      na.rm = TRUE) != 0)
                        },
                        FUN.VALUE = integer(length = 1)))
             },
             FUN.VALUE = integer(length = nr + 1)),
      nrow = nr+1, ncol = nb_iter)
  }
  if (method == "block") {
   #browser()
    if (is.null(block)) {
      stop(paste0("You must provide block memberships to use method \"block\"."))
    }
    Z <- block
    connec <- vapply(
      X = seq_len(max(Z)),
      FUN = function(k) {
        mean(rowMeans(A[Z == k,]))
      },
      FUN.VALUE = double(length = 1)
    )
    if (length(unique(rownames(A))) != nrow(A)) {
      rownames(A) <- seq(nr)
    }
    ord <- switch(ext_seq,
           "increasing" = order(connec),
           "decreasing" = rev(order(connec)),
           seq(max(Z))
    )
    Z <- ord[Z]
    eff <- c(0,table(Z))
    A <- A[order(Z),]
    sec_ext <- matrix(
      vapply(X = seq_len(nb_iter),
             FUN = function(i) {
               prim_ext <-unlist(
                 sapply(seq(max(Z)),
                        function(k) eff[k]  + sample(sum(Z == k)), simplify = TRUE))
               c(sum(colSums(A, na.rm = TRUE) != 0),
                 vapply(X = seq.int(nrow(A)),
                        FUN = function(m) {
                          sum(colSums(A[prim_ext[-seq_len(m)], ,drop = FALSE],
                                      na.rm = TRUE) != 0)
                        },
                        FUN.VALUE = integer(length = 1)))
             },
             FUN.VALUE = integer(length = nr + 1)),
      nrow = nr+1, ncol = nb_iter)
  }
  robustness <- list(model = "empirical",
                     method = method,
                     ext_seq = ext_seq,
                     power = power,
                     fun = rowMeans(sec_ext)/nc,
                     sd = vapply(seq_len(nr+1),
                                 function(i) stats::sd(sec_ext[i,]/nc),
                                 FUN.VALUE = double(length = 1)),
                     auc = sum(rowMeans(sec_ext)/nc)/nr)
  class(robustness) <- "robber"

  return (robustness)
}


#' Plot function pour robber class
#'
#' @param x An object of class robber
#' @param type Type of the plot, just function
#' @param add Whether it should be added to a previous ggplot
#' @param ... \code{lty}, \code{col}, ... to be passed to \code{geom_step()} in
#' \code{ggplot2}
#'
#' @return A \code{ggplot2} object
#' @export
#'
#' @examples
#' data("hostparasite", package = "robber")
#' rob <- robustness_emp(hostparasite)
#' plot(rob)
plot.robber <- function(x, type = "function", add = FALSE, ...) {
  r <- s <- y <-  NULL
  if (x$model == "empirical") {
    color <- "blue"
    lt <- "dashed"
  } else {
    color = "black"
    lt <- "solid"
  }
  if (is.null(x$sd)) x$sd <- 0
  if(! add) {
    requireNamespace("pammtools", quietly = TRUE)
    p <- data.frame(r = seq(0, 1, length.out = length(x$fun)),
               y = x$fun,
               s = x$sd) %>%
      ggplot2::ggplot(ggplot2::aes(x = r, y = y)) +
      pammtools::geom_stepribbon(
        ggplot2::aes(ymin = pmax(0, y - 2*s), ymax = pmin(y + 2*s, 1)),
        alpha = .2) +
      ggplot2::geom_step(col = color, linetype = lt, ... = ...) +
      ggplot2::xlab("Primary extinctions") +
      ggplot2::ylab("Robustness") #+
    # ggplot2::coord_equal() +
    # ggplot2::ggtitle(paste(x$model, x$method,
    #                        x$ext_seq, "robustness: ", round(x$auc,3)))
    p
  }
  if(add) {
    p <- ggplot2::geom_step(
      data = data.frame(x = seq(0, 1, length.out = length(x$fun)),
                                y = x$fun),
              mapping = ggplot2::aes(x = x, y = y), ... = ...)
  }
  return(p)
}
