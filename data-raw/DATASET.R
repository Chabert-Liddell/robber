## code to prepare `DATASET` dataset goes here

# my_files <- dir(path = "~/Documents/ecological_network/weboflife/bipartite/")
# my_files <- my_files[1:(length(my_files)-2)]
# web_of_life <- lapply(
#   X = seq_along(my_files),
#   FUN = function(i) {
#     df <- read.csv(file = paste0("~/Documents/ecological_network/weboflife/bipartite/",
#                                  my_files[i]), header = TRUE, row.names = 1)
#     A <- as.matrix(df)
#     if (stringr::str_sub(my_files[i], 1, 4) == "A_HP") {
#       A <- A[,2:ncol(A)]
#     }
#     A[A!=0] <- 1
#     return(list(
#       net = A,
#       nr = nrow(A),
#       nc = ncol(A),
#       dens = mean(A),
#       id = stringr::str_sub(my_files[i], 1, -5))
#     )
#   }
# )
#
#
# usethis::use_data(web_of_life, overwrite = TRUE)
#


df <- read.csv(file = "~/Documents/ecological_network/weboflife/bipartite/A_HP_044.csv",
      header = TRUE, row.names = 1)
A <- as.matrix(df)
A <- A[,2:ncol(A)]
A[A!=0] <- 1
hostparasite <- A

df <- read.csv(file = "~/Documents/ecological_network/weboflife/bipartite/M_SD_013.csv",
               header = TRUE, row.names = 1)
A <- as.matrix(df)
A[A!=0] <- 1
seeddispersal <- A

df <- read.csv(file = "~/Documents/ecological_network/weboflife/bipartite/M_PL_039.csv",
               header = TRUE, row.names = 1)
A <- as.matrix(df)
A[A!=0] <- 1
pollination <- A

usethis::use_data(hostparasite, seeddispersal, pollination, overwrite = TRUE)
