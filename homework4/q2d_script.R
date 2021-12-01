source("funcs.R")
library(parallel)
library(doParallel)

n.core <- 12
sample.size <- 5000

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  dd <- 128
  par <- c(0, 0.5, 0.5)
  m.by <- 10000
} else if (length(args) == 5) {
  dd <- as.numeric(args[1])
  par <- c(as.numeric(args[2]), as.numeric(args[3]), as.numeric(args[4]))
  m.by <- as.numeric(args[5])
} else {
  stop("please provide 0 or 4 arguments: id, layer, q, B")
}

cl <- parallel::makeCluster(n.core)
doParallel::registerDoParallel(cl)

result <- perfect.gibbs.multi.grid.sampler(c(dd, dd), par, n.core = n.core,
                                           m.by = m.by, sample.size = sample.size)

saveRDS(result, file = paste("rds_folder/r", dd, paste(par, collapse = "-"), "r.rds", sep = "-"))

parallel::stopCluster(cl)
