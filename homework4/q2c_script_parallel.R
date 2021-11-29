library(parallel)
library(doParallel)

source("funcs.R")

# question 2
x_prob_parallel <- function(i, j, z, par, n.core) {

  i.split <- split(i,
                   ceiling(seq_along(i) / ceiling(length(i) / n.core) ))
  j.split <- split(j,
                   ceiling(seq_along(i) / ceiling(length(i) / n.core) ))

  result <- foreach(b = 1:n.core, .combine = 'c') %dopar% {
    source("funcs.R")
    x_prob_vec(i.split[[b]], j.split[[b]], z, par)
  }

  return(result)

}


n.core <- 12

cl <- parallel::makeCluster(n.core)
doParallel::registerDoParallel(cl)

### (c)
dim.z <- c(512, 512)
# dim.z <- c(16, 16)
par <- c(0, 1, 1)
# par <- c(0.25, 0.5, 0.5)

result.c <- find.perfect.sample(dim.z, par, n.core=n.core, seed = 17381128, m.start = 25)
saveRDS(result.c, file = "h4_q2c_result/c_perfect_sample_512_parallel.rds")

parallel::stopCluster(cl)

