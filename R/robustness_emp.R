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
#' \code{300}
#' @param net A string, the type of network. Only "bipartite".
#'
#' @return A list of size 2:
#' * \code{$fun} is the robustness function, a vector of size \code{(nrow(A) +1)}
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
                           nb_iter = 300, net = "bipartite") {
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
  if (ext_seq == "uniform") {
    sec_ext <- matrix(
      vapply(X = seq_len(nb_iter),
             FUN = function(i) {
               prim_ext <- sample(x = seq_len(nr), size = nr, replace = FALSE)
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
  if (ext_seq == "increasing") {
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
  if (ext_seq == "decreasing") {
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
  robustness <- list(fun = rowMeans(sec_ext)/nc,
                     auc = sum(rowMeans(sec_ext)/nc)/nr)
  return (robustness)
}
