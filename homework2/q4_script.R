args = commandArgs(trailingOnly=TRUE)


lsb = load("data/LSB.rda")

source("parallel_setup.R")

fun3 <- function(LSBs, id, layer, B = 50, q = 0, FUN = count_neighbor, core = 6) {


  # cl <- parallel::makeCluster(core)
  # doParallel::registerDoParallel(cl)

  dat.original <- LSBs[[id]]$lsb[, , layer]
  dat.original.neighbor <- FUN(dat.original)

  result <- list()
  result[[1]] <- compute_MPLE_ratio(dat.original.neighbor, q)

  if(q == 0) {
    tmp <- foreach(b=1:B, .combine = 'list') %dopar% {
      source("parallel_setup.R")
      smp = matrix(rbinom(nrow(dat.original)*ncol(dat.original), 1, 0.5), nrow = nrow(dat.original))
      dat.neighbor <- FUN(smp)

      compute_MPLE_ratio(dat.neighbor, q)
      # result[[b+1]] <- compute_MPLE_ratio(dat.neighbor, q)
    }
  } else if (q == 1) {
    opt.r <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
                   z.neighbor = dat.original.neighbor, beta_zero = TRUE,
                   method = "L-BFGS-B", lower = 0.0001)
    p = exp(opt.r$par)/(1 + exp(opt.r$par))
    tmp <- foreach(b = 1:B, .combine = 'list') %dopar% {
      source("parallel_setup.R")
      smp = matrix(rbinom(nrow(dat.original)*ncol(dat.original), 1, p), nrow = nrow(dat.original))
      dat.neighbor <- FUN(smp)

      compute_MPLE_ratio(dat.neighbor, q)
      # result[[b+1]] <- compute_MPLE_ratio(dat.neighbor, q)
    }
  } else {
    n <- nrow(dat.original)
    m <- ncol(dat.original)
    opt.r <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
                   z.neighbor = dat.original.neighbor, alpha_zero = TRUE,
                   method = "L-BFGS-B", lower = 0.0001)
    beta <- opt.r$par
    mask <- matrix(1,n,m)
    neigh <- getNeighbors(mask, c(2,2,0,0))
    block <- getBlocks(mask, 2)
    k <- 2
    tmp <- foreach(b = 1:B, .combine = 'list') %dopar% {
      source("parallel_setup.R")
      result <- swNoData(beta = beta,k = k,neigh = neigh, block = block)
      z <- matrix(max.col(result$z)[1:nrow(neigh)], nrow=nrow(mask))
      smp = z - 1
      dat.neighbor <- FUN(smp)

      compute_MPLE_ratio(dat.neighbor, q)
      # result[[b+1]] <- compute_MPLE_ratio(dat.neighbor, q)
    }
  }
  result <- c(result, tmp)

  # parallel::stopCluster(cl)

  return(result)
}

cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)

if (length(args) == 0) {
  id <- 7
  layer <- 1
  q <- 0
  B <- 2
} else if (length(args) == 4) {
  id <- args[1]
  layer <- args[2]
  q <- args[3]
  B <- args[4]
} else {
  stop("please provide 0 or 4 arguments: id, layer, q, B")
}

system.time({
  ress <- fun3(LSBs, id = id, layer = layer, B = B, q = q, FUN = count_neighbor_dl)
})

saveRDS(ress, paste("results/result", id, layer, q, B, "lsb.rds", sep='-'))

parallel::stopCluster(cl)





