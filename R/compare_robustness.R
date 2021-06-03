normalize_lbm_con <- function(dens, con, pi, rho, nr = NULL, nc = NULL) {
  con <- dens * con / as.vector(pi %*% con %*% rho)
  if(any(con>=1)) warning("Connectivity > 1 set to 1.")
  pmin(con, 1-1e-6)
}


#' @title Compare the robustness for different LBM parameters
#'
#' @description This function is designed to be used with a list of parameters
#'   as the one given by the function \code{get_lbm_param()}. The default
#'   parameters are designed such that the return value for a uniform extinction
#'   sequence is bounded above by 0.5.
#'
#' @param list_param A list of list of LBM parameters, each list must contain at
#'   list \code{con}, \code{pi} and \code{rho} entries, such as the one returned
#'   by the function \code{get_lbm_param()}. The size of the network is
#'   optional.
#' @param dens The density (connectance) used to compare the network. The
#'   \code{com} parameters of \code{list_param} will be normalized so that the
#'   average probability of interaction in the network is equal to \code{dens}.
#'   Default to  0.0156.
#' @param new_nr The number of rows of the normalized networks. Default to 100.
#' @param new_nc The number of columns of the normalized networks.
#'   Default to 100.
#' @param ext_seq The distribution of the primary extinction sequence used to
#'  compare the networks. One of \code{c("uniform", "decreasing", "increasing")}.
#'  Default to \code{"uniform"}.
#'
#' @return A list of the sames length as \code{list_param}, the robustness (AUC)
#'   after normalization for all set of parameters.
#' @export
#'
#' @examples
#' par1 <- list(con = matrix(c(.5,.3,.3,.1), 2, 2),
#'              pi  = c(.25,.75),
#'              rho = c(1/3, 2/3))
#' par2 <- list(con = matrix(c(.4,.1,.1,.4), 2, 2),
#'              pi  = c(.25,.75),
#'              rho = c(1/2, 1/2))
#´ compare_robustness(list_param = list(par1, par2))
#' compare_robustness(list_param = list(par1, par2), ext_seq = "decreasing")
#' compare_robustness(list_param = list(par1, par2), ext_seq = "increasing")
compare_robustness <- function(list_param,
                               dens = 0.0156,
                               new_nr = 100,
                               new_nc = 100,
                               ext_seq = "uniform") {
  purrr::map(.x = list_param,
             .f = function(par) {
               par$con <- do.call(normalize_lbm_con, c(par, dens = dens))
               par$nr <- new_nr
               par$nc <- new_nc
               if (ext_seq == "uniform") {
                 do.call(auc_robustness_lbm, par)
               } else {
                 do.call(robustness_lbm, c(par, ext_seq = ext_seq))$auc
               }
             })
}
