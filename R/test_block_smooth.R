data("web_of_life")
net <- web_of_life[[60]]$net

smooth_emp <- function(A,  nb_iter = 300)  {
  nr <- nrow(A)
  nc <- ncol(A)
  sec_ext <- matrix(
  vapply(X = seq_len(nb_iter),
         FUN = function(i) {
           prim_ext <- sample(x = seq_len(nr), size = nr,
                              prob = rowSums(A)/sum(A),replace = FALSE)
       #    prim_ext <- rev(prim_ext)
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
  return(list(fun = rowMeans(sec_ext)/nc,
              auc = sum(rowMeans(sec_ext)/nc)/nr))
}


net <- web_of_life[[60]]$net
smooth_emp(net)
param <- get_lbm_param(net)
rob_dec <- do.call(robustness_lbm, c(param, ext_seq = "decreasing"))
rob_inc <- do.call(robustness_lbm, c(param, ext_seq = "increasing"))
rob_dec$fun/rob_dec$fun[1]


smooth_inc <- purrr::map_dbl(web_of_life, ~ smooth_emp(.x$net)$auc)
smooth_dec <- purrr::map_dbl(web_of_life, ~ smooth_emp(.x$net)$auc)
####
####
####
param <- readRDS("vignettes/wol_lbm_parameters.rds")
index <- map_dbl(param, ~ choose(.x$nr + nrow(.x$con) -1, .x$nr))
index <- index < 1e7
indexnoer <- map_lgl(param, ~ nrow(.x$con) > 1)

indexbis <- index & indexnoer
rob_bm_dec <- purrr::map_dbl(seq_along(web_of_life),
                             ~ ifelse(index[.x], do.call(robustness_lbm, c(param[[.x]][1:5], ext_seq = "decreasing"))$auc,
                                      0))
rob_bm_inc <- purrr::map_dbl(seq_along(web_of_life),
                             ~ ifelse(index[.x], do.call(robustness_lbm, c(param[[.x]][1:5], ext_seq = "increasing"))$auc,
                                      0))

data.frame(x = smooth_inc[index], y = rob_bm_inc[index]) %>%
  ggplot2::ggplot(ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point() +
  geom_abline(slope = 1, intercept = 0)
cor(x = smooth_inc[index], y = rob_bm_inc[index], method = "spearman")


data.frame(x = smooth_dec[indexbis], y = rob_bm_dec[indexbis]) %>%
  ggplot2::ggplot(ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point() +
  geom_abline(slope = 1, intercept = 0)
cor(x = smooth_dec[indexbis], y = rob_bm_dec[indexbis], method = "spearman")

tb_rob <- readRDS("vignettes/tb_robustness_wol.rds")



rob_bm_inc - tb_rob$BM_inc



data.frame(x = tb_rob$emp_inc[index], y = rob_bm_inc[index]) %>%
  ggplot2::ggplot(ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point() +
  geom_abline(slope = 1, intercept = 0)
cor(x = tb_rob$emp_inc[index], y = rob_bm_inc[index], method = "spearman")

data.frame(x = tb_rob$BM_inc[index], y = smooth_inc[index]) %>%
  ggplot2::ggplot(ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point() +
  geom_abline(slope = 1, intercept = 0)

cor(data.frame(emp = tb_rob$emp_inc[index], smooth = smooth_inc[index],
    bm = rob_bm_inc[index], bm.old = tb_rob$BM_inc[index]),
    method = "spearman", use = "pairwise.complete.obs")
cor(data.frame(emp = tb_rob$emp_dec[index], smooth = smooth_dec[index],
               bm = rob_bm_dec[index], bm.old = tb_rob$BM_dec[index]),
    method = "spearman", use = "pairwise.complete.obs")


### Essaie pour pas de block vide
###
x <- as.matrix(partitions::compositions(5, 3))
dmult <- apply(x, 2, function(x) dmultinom(x, size = 5, c(2,3,5)/10))

y <- x[,colSums(x==0) == 0]

dmult2 <- apply(y, 2, function(y) dmultinom(y, size = 5, c(2,3,5)/10))


rob_z_dec <- purrr::map_dbl(seq_along(web_of_life),
                             ~ ifelse(index[.x], do.call(robustness_lbm, c(param[[.x]][1:5], ext_seq = "decreasing"))$auc,
                                      0))
rob_z_inc <- purrr::map_dbl(seq_along(web_of_life),
                             ~ ifelse(index[.x], do.call(robustness_lbm, c(param[[.x]][1:5], ext_seq = "increasing"))$auc,
                                      0))


data.frame(x = tb_rob$emp_inc[index], y = rob_z_inc[index]) %>%
  ggplot2::ggplot(ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  ggplot2::geom_point(alpha = .5) +
  ggplot2::scale_color_gradient( low = "grey50", high = "red") +
  ggplot2::xlim(c(0.5,1)) +
  ggplot2::ylim(c(0.5,1)) +
  ggplot2::xlab(label = "AUC(R) by Monte Carlo") +
  ggplot2::ylab(label = "AUC(R) by LBM Expectation") +
  ggplot2::coord_fixed() +
  ggplot2::theme_minimal(base_size = 15, base_rect_size = .25)
cor(x = tb_rob$emp_inc[index], y = rob_bm_inc[index], method = "spearman")

data.frame(x = rob_bm_inc[index], y = rob_z_inc[index]) %>%
  ggplot2::ggplot(ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  ggplot2::coord_fixed()

cor(data.frame(emp = tb_rob$emp_inc[index], smooth = smooth_inc[index],
               bm = rob_bm_inc[index], bm.old = tb_rob$BM_inc[index],
               z = rob_z_inc[index]),
    method = "spearman", use = "pairwise.complete.obs")
cor(data.frame(emp = tb_rob$emp_dec[index], smooth = smooth_dec[index],
               bm = rob_bm_dec[index], bm.old = tb_rob$BM_dec[index],
               z = rob_z_dec[index]),
    method = "spearman", use = "pairwise.complete.obs")
