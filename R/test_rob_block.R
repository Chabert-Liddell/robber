# con <-  matrix(c(.1, .5), nrow=2)
# pi <- c(2/3, 1/3)
# rho <- 1
# nc <- 10
# nr <- 100
#
# fn1  <- function() {
#   X <- t(as.matrix(expand.grid(
#     replicate(length(pi), seq(0,nr), simplify = FALSE))))
#   X <- X[, colSums(X) == nr]
#   pmult <- apply(X, 2, function(x) dmultinom(x, prob = pi))
#   rob <- rep(1,nr+1)
#   for (q in seq_along(rho)) {
#     for(m in seq(0,nr)) {
#       rob[m+1] <- rob[m+1] - rho[q] *
#         sum(pmult *
#               apply(X, 2,
#                     function(x) prod((1-con[,q])**pmax(0,pmin(x , cumsum(x) - m)))))
#     }
#   }
#   rob
# }
# createX <- function() {
#   X <- t(as.matrix(expand.grid(
#     replicate(length(pi), seq(0,nr), simplify = FALSE))))
#   X <- X[, colSums(X) == nr]
#   pmult <- apply(X, 2, function(x) dmultinom(x, prob = pi))
#   cumX <- apply(X, 2, cumsum)
# }
#
# fn2 <- function() {
#   X <- t(as.matrix(expand.grid(
#     replicate(length(pi), seq(0,nr), simplify = FALSE))))
#   X <- X[, colSums(X) == nr]
#   pmult <- apply(X, 2, function(x) dmultinom(x, prob = pi))
#   rob <- rep(1,nr+1)
#   cumX <- apply(X, 2, cumsum)
#   for (q in seq_along(rho)) {
#     for (m in seq(0, nr)) {
#       rob[m+1] <- rob[m+1] - rho[q] * sum(pmult * exp(rowSums(log(
#         vapply(seq_along(pi),
#                function(i) (1-con[i,q]) ** pmax(0,pmin(X[i,], cumX[i,] - m)),
#                FUN.VALUE = rep(0, ncol(X)))))))
#     }
#   }
#   rob
# }
#
# fn3 <- function() {
#   X <- t(as.matrix(expand.grid(
#     replicate(length(pi), as.double(seq(0,nr)), simplify = FALSE))))
#   X <- X[, colSums(X) == nr]
#   pmult <- apply(X, 2, function(x) dmultinom(x, prob = pi))
#   rob <- rep(1,nr+1)
#   cumX <- Rfast::colCumSums(X)
#   for (q in seq_along(rho)) {
#     for (m in seq(0, nr)) {
#       rob[m+1] <- rob[m+1] - rho[q] * sum(pmult *
#       Rfast::rowprods(
#         vapply(seq_along(pi),
#                function(i) (1-con[i,q])**pmax(0,pmin(X[i,], cumX[i,] - m)),
#                FUN.VALUE = rep(0, ncol(X)))))
#     }
#   }
#   rob
# }
#
# fn4 <- function() {
#   X <- partitions::compositions(nr, length(pi))
#   pmult <- apply(X, 2, function(x) dmultinom(x, prob = pi))
#   rob <- rep(1,nr+1)
#   storage.mode(X) <- "numeric"
#   cumX <- Rfast::colCumSums(X)
#   for (q in seq_along(rho)) {
#     for (m in seq(0, nr)) {
#       rob[m+1] <- rob[m+1] - rho[q] * sum(pmult *
#                                             Rfast::rowprods(
#                                               vapply(seq_along(pi),
#                                                      function(i) (1-con[i,q])**pmax(0,pmin(X[i,], cumX[i,] - m)),
#                                                      FUN.VALUE = rep(0, ncol(X)))))
#     }
#   }
#   rob
# }
#
# microbenchmark::microbenchmark(fn1(), fn2(), fn3(), fn4(),
#                                check = "equal", times = 10)
#
# microbenchmark::microbenchmark(
#   (as.matrix(expand.grid(
#     replicate(length(pi), seq(0,nr), simplify = FALSE)))),
#   t(as.matrix(expand.grid(
#     replicate(length(pi), as.double(seq(0,nr)), simplify = FALSE)))),
#   partitions::compositions(nr, length(pi)), times = 10)
#
# microbenchmark::microbenchmark(
#   apply(X, 2, cumsum),
#   apply(t(X), 1, cumsum),
#   Rfast::colCumSums(X),
#   check = "equal"
# )
#
#
# microbenchmark::microbenchmark(
#   apply(X, 2, function(x) dmultinom(x, prob = pi)),
#   mc2d::dmultinomial(t(X), prob = pi), check = "equal", times = 100)
#
# fn1()
# fn2()
# )
# (1-con[1,q]) ** pmax(0,pmin(X[1,], cumX[1,] - m))
# (1-con[2,q]) ** pmax(0,pmin(X[2,], cumX[2,] - m))
#
#
# X <- t(as.matrix(expand.grid(0:4, 0:4)))
# X <- X[, colSums(X) == 4]
# #X <- rbind(X, 4:4 - colSums(X));
# dimnames(X) <- list(letters[1:2], NULL)
# X
#
# pmult <- apply(X, 2, function(x) dmultinom(x, prob = c(1/3,2/3)))
#
# # m = 0
# 1 - sum(pmult * ((1-con[1,])**X[1,] * (1-con[2,])**X[2,]))
# # m = 1
# X1 <-  X
# X1[1, X[2,] == 0] <- X[1, X[2,] == 0] - 1
# X1[2, X[2,] != 0] <- X[2, X[2,] != 0] - 1
# 1 - sum(pmult * ((1-con[1,])**X1[1,] * (1-con[2,])**X1[2,]))
# # m = 2
# X2 <-  X1
# X2[1, X1[2,] == 0] <- X1[1, X1[2,] == 0] - 1
# X2[2, X1[2,] != 0] <- X1[2, X1[2,] != 0] - 1
# 1 - sum(pmult * ((1-con[1,])**X2[1,] * (1-con[2,])**X2[2,]))
# # m = 3
# X3 <-  X2
# X3[1, X2[2,] == 0] <- X2[1, X2[2,] == 0] - 1
# X3[2, X2[2,] != 0] <- X2[2, X2[2,] != 0] - 1
# 1 - sum(pmult * ((1-con[1,])**X3[1,] * (1-con[2,])**X3[2,]))
#
# #m = 4
# X4 <-  X3
# X4[1, X3[2,] == 0] <- X3[1, X3[2,] == 0] - 1
# X4[2, X3[2,] != 0] <- X3[2, X3[2,] != 0] - 1
# 1 - sum(pmult * ((1-con[1,])**X4[1,] * (1-con[2,])**X4[2,]))
#
#
# zpl <- pbinom(q = seq(0,nr), size = nr, prob = 2/3)
# zpl <- cbind(zpl, pbinom(q = seq(0,nr), size = nr, prob = 1))
#
# Y <- t(as.matrix(expand.grid(0:3, 0:3)))
# Y <- Y[, colSums(Y) == 3]
# apply(Y[,Y[2,] == 0], 2, function(x) dmultinom(x, prob = c(1/3,2/3)))
# zpl[1,1] * apply(Y, 2, function(x) dmultinom(x, prob = c(1/3,2/3))) +
#   (1 - zpl[1, 1])
#
#
#
# # m = 3
# 1 - zpl[4, 1] * (1-con[1,]) - (1-zpl[4,1]) * (1-con[2,])
#
# 1 - zpl[3, 1] * (1-con[1,]) -
#   (1 - zpl[3,1]) *
#   (dbinom(3, 4, 2/3) * (1-con[1,]) + dbinom(4, 4, 2/3) * (1-con[2,]))/
#   sum(dbinom(3:4, 4, 2/3))
#
# # m = 3
# 1 - zpl[3, 1] * (1-con[1,])**2 -
#   (1 - zpl[3,1]) *
#   (dbinom(3, 4, 2/3) * (1-con[1,]) * (1-con[2,]) + dbinom(4, 4, 2/3) * (1-con[2,])**2)/
#   sum(dbinom(3:4, 4, 2/3))
# n <- 50
# rob <- rep(1, n+1)
# for (z1 in seq(0, n)) {
#   z2 <- n - z1
#   z <- c(z1, z2)
# #  for (z2 in seq(0, 4-z1)) {
#   for (m in seq(0,n)) {
#     rob[m+1] <- rob[m+1] - dmultinom(c(z1, z2), prob = rev(pi)) *
#       prod((1-rev(con))**pmax(0,pmin(z , cumsum(z) - m)))
#   }
# }
#
#
# (2/3)**4
# 1-(1/3)**4
# .5**(nr:0)*pbinom(q = seq(0,nr), size = nr, prob = 1/3) + .1**(nr:1) * (1-pbinom(q = seq(nr), size = nr, prob = 1/3) )
# library(multicool)
# gc <- multicool::genComp(3, 3, addZeros = TRUE)
# purrr::map(multicool::genComp(100, 3, addZeros = TRUE),
#            ~ allPerm(initMC(.x))) %>% unlist() %>% matrix(ncol = 3, byrow = TRUE)
#
#
# pi[1]**4
# choose(4, 1) * pi[1]**3 * pi[2]
# choose(4, 2) * pi[1]**2 * pi[2]**2
# choose(4, 3) * pi[1]**1 * pi[2]**3
# choose(4, 4) * pi[1]**0 * pi[2]**4
#
#
# choose(nr, seq(0, nr)) * pi[1]**seq(nr, 0) * pi[2]**seq(0, nr)
# cumsum(choose(nr, seq(0, nr)) * pi[1]**seq(nr, 0) * pi[2]**seq(0, nr))
#
#
# ext_prob <- pbinom(q = seq(0,nr), size = nr, prob = 2/3)
# ext_prob2 <- pbinom(q = seq(0,nr), size = nr, prob = 1/3)
#
# (1-ext_prob[1]) * 0.5 + ext_prob[1]*.1
# (1-ext_prob[2])*.5 + ext_prob[2] * (pi[1]*con[1] + pi[2] * con[2])
#
#
# 1- (1-ext_prob2[1]) * (1-con[1]) - ext_prob2[1]*(1-con[2]) #inc
# 1 - ext_prob[1] * (1-con[1]) - (1-ext_prob[1])*(1-con[2]) #dec
#
# robber::robustness_lbm(con, pi, rho, nr, nc, ext_seq = "increasing")
#
#
#
# (1-pi[2]**nr) * con[1] +pi[2]**nr * con[2]
# (1-pi[1]**nr) * con[2] + pi[1]**nr *con[1]
#
# diff(ext_prob2)
# 1 - (1-ext_prob2[2]) * (1-0.5)**2 - ext_prob2[1]*(1-.1)**2  -
#   diff(ext_prob2)[1] * (1-con[1]) * (1-con[2]) #inc
#
# 1 - ext_prob[2] * (1-0.5)**2 - (1-ext_prob[1])*(1-.1)**2  -
#   diff(ext_prob)[1] * (1-con[1]) * (1-con[2]) #dec
#
# 1 - dbinom(4, nr, 1/3) * (1-con[1])**2 -
#   (1-pbinom(1, nr, 2/3)) * (1 - con[2])**2 -
#   dbinom(1, nr, 2/3) * (1-con[1]) * (1-con[2])
#
# ## Essaie a 3
#
# con <-  matrix(c(.5,  .3, .1, .3, .1, .1), nrow = 3)
# pi <- c(1/6, 2/6, 3/6)
# rho <- c(.5,.5)
# nc <- 10
# nr <- 200
#
# rob <- rep(1, nr+1)
#
#
# X <- t(as.matrix(expand.grid(
#   replicate(length(pi), seq(0,nr), simplify = FALSE))))
# X <- X[, colSums(X) == nr]
#
# pmult <- apply(X, 2, function(x) dmultinom(x, prob = pi))
#
# rob <- rep(1,nr+1)
#
# for (q in seq_along(rho)) {
#   for(m in seq(0,nr)) {
#     rob[m+1] <- rob[m+1] -
#       rho[q] * sum(pmult *
#                      apply(X, 2,
#                            function(x) prod((1-con[,q])**pmax(0,pmin(x , cumsum(x) - m)))))
#   }
# }
#
#
#
# robber::robustness_lbm(con, pi, rho, nr, nc, ext_seq = "increasing")
#
# # m = n-1
# # 3 cas il reste 1 du block 1, 1 du block 2, 1 du block 3
# 1 - pbinom(3, nr, 5/6) * (1-con[1]) -
#   (dbinom(4, nr, 5/6) - dbinom(4, nr, 3/6)) * (1 - con[2])-
#   dbinom(4, 4, 3/6) * (1 - con[3])
#
#
# # m = 1
# # n-m = 3
#
# 1 -
#   pbinom(1, 4, 3/6) *(1 - (pi[1:2]%*%con[1:2])/sum(pi[1:2]) )**3 -
#   (1 - pbinom(1, 4, 3/6))
#
#
# X <- t(as.matrix(expand.grid(0:nr, 0:nr, 0:nr)))
# X <- X[, colSums(X) == 4]
# pmult <- apply(X, 2, function(x) dmultinom(x, prob = pi))
#
# rob <- rep(1,nr+1)
#
# for (q in seq_along(rho)) {
#   for(m in seq(0,nr)) {
#     rob[m+1] <- rob[m+1] - rho[q] * sum(pmult *
#       apply(X, 2,
#             function(x) prod((1-con[,q])**pmax(0,pmin(x , cumsum(x) - m)))))
#   }
# }
# X1 <-  X
# X1[1, X[2,] == 0] <- X[1, X[2,] == 0] - 1
# X1[2, X[2,] != 0] <- X[2, X[2,] != 0] - 1
# 1 - sum(pmult * ((1-con[1,])**X1[1,] * (1-con[2,])**X1[2,]))
#
#
# #sum(apply(z, 1, function(x) 5 %in% unlist(sapply(1:4, function(i) combn(x, i, sum)))))
# X <- X[, colSums(X) <= 4]
# X <- rbind(X, 4:4 - colSums(X));
# dimnames(X) <- list(letters[1:3], NULL)
# X
# round(apply(X, 2, function(x) dmultinom(x, prob = c(3/6,2/6,1/6))), 3)
#
#
#
# pbinom(0:4, 4, 3/6)
# pbinom(0:4, 4, 5/6)
# pbinom(0:4, 4, 6/6)
# pbinom(0:4, 4, 1/6)[1] * dbinom(1, 4, 2/5)
# for (k in seq(3)) {
#
# }
# pmult[X[1,] <= 1] * (1-con)**
#
# X[1,] -1
#
#
# combn(4, 3)
#
# (1-con[1])
#
# ###########################
# block_lbm <- function( con, pi, rho, nr, nc) {
#   result <- rep(0,nr+1)
#   K <- length(pi)
#   order_prob <- matrix(0, nr+1, K)
#   order_prob[nr+1, K] <- 0
#   for (k in seq_len(K-1)) {
#     for (m in seq.int(0, nr-1)) {
#       order_prob[m+1, k] <-
#         sum(stats::dbinom(x = seq.int(0, nr-m-1), size = nr,
#                           prob = 1-sum(pi[1:k])) -
#               stats::dbinom(x = seq.int(0, nr-m-1), size = nr,
#                             prob = 1 - sum(pi[1:k]) + pi[k]))
#     }
#   }
#   order_prob[,K] <- 1 - rowSums(order_prob)
#   tmp_rob <- matrix(0, nr+1, length(rho))
#   for(m in seq.int(0, nr)) {
#     for(k in seq_len(K)) {
#       tmp_rob[m+1,] <- tmp_rob[m+1,] +
#         (( (pi[1:k]) %*% (1-con[1:k, , drop=FALSE])/sum(pi[1:k]))**(m) *
#            order_prob[nr - m+1, k])
#     }
#   }
#   result[1:(nr+1)] <- as.vector(tmp_rob %*% rho)[(nr+1):1]
#   return(list(fun = 1 - result,
#               auc = sum(1 - result)/nr))
# }
#
# block_lbm(con, pi, rho, nr, nc)
