#' @title Estimate the parameters of a LBM
#'
#' @param A A matrix, the incident matrix of the network
#' @param model_size A vector of size 2, the size of the model parameters.
#' If \code{NULL}, the default, model will be chosen as the one with the
#' highest ICL criterion among all fitted models during the inference.
#' @param ... Additional arguments pass to the inference function of
#' package \code{GREMLIN} if there are missing value and \code{blockmodels} if
#' none.
#' @return A list of the LBM parameters. The list is set
#' @export
#'
#' @examples
#' # When the size of the model is unknown:
#' data("seeddispersal")
#' param <- get_lbm_param(seeddispersal, ncores = 1L)
#' do.call(auc_robustness_lbm, param)
#'
#' # For a fixed number of blocks (useful for comparison)
#' param <- get_lbm_param(seeddispersal, model_size = c(1, 3), ncores = 1L)
#' do.call(auc_robustness_lbm, param)
#'
#' # For data with missing observations
#' \donttest{
#' A <- seeddispersal
#' A[sample(seq_len(nrow(A)), 5),
#'   sample(seq_len(ncol(A)), 10)] <- NA
#' param <- get_lbm_param(A, nbCores = 1L)
#' do.call(auc_robustness_lbm, param)}
get_lbm_param <- function(A, model_size = NULL, ...) {
  my_net <- GREMLINS::defineNetwork(A, typeInter = "inc",
                                   rowFG = "row", colFG = "col")
  if (any(is.na(A))) {
    if (is.null(model_size)) {
      my_lbm <- GREMLINS::multipartiteBM(list(my_net),
                                        v_distrib = "bernoulli",
                                        namesFG = c("row", "col"),
                                        verbose = FALSE, ...)
    } else {
      my_lbm <- GREMLINS::multipartiteBMFixedModel(list(my_net),
                                                  v_distrib = "bernoulli",
                                                  namesFG = c("row", "col"),
                                                  v_K = model_size,
                                                  verbose = FALSE, ...)
    }
    return(list(con = my_lbm$fittedModel[[1]]$paramEstim$list_theta$rowcol,
                pi = my_lbm$fittedModel[[1]]$paramEstim$list_pi$row,
                rho = my_lbm$fittedModel[[1]]$paramEstim$list_pi$col,
                nr = nrow(A),
                nc = ncol(A)))
  } else {
    if (is.null(model_size)) {
      my_lbm <-
        blockmodels::BM_bernoulli("LBM", A, verbosity = 0, plotting = "", ...)
    } else {
      my_lbm <-
        blockmodels::BM_bernoulli("LBM", A, verbosity = 0, plotting = "",
                                  exploration_direction = model_size, ...)
    }
    my_lbm$estimate()
    ifelse (is.null(model_size),
            model_size <- which.max(my_lbm$ICL),
            model_size <- which(sapply(
              seq(1, length(my_lbm$memberships)),
              function(i) {
                ncol(my_lbm$memberships[[i]]$Z1) == model_size[1] &
                  ncol(my_lbm$memberships[[i]]$Z2) == model_size[2]
                }) == TRUE))
    if (purrr::is_empty(model_size)) return(NULL)
    return(list(con = my_lbm$model_parameters[[model_size]]$pi,
                pi = colMeans(my_lbm$memberships[[model_size]]$Z1),
                rho = colMeans(my_lbm$memberships[[model_size]]$Z2),
                nr = nrow(A),
                nc = ncol(A)))
  }
}
