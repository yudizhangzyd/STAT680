library(bayesImageS)
# lsb = load("../STAT680/local/LSB.rda")
lsb = load("./LSB.rda")

find_neighbor <- function(z, i, j) {
  dim.z <- dim(z)

  if(i > dim.z[1] | j > dim.z[2]) return(NULL)

  count <- 0

  neigh.idx <- list()
  if(i - 1 >= 1) {
    count <- count + 1
    neigh.idx[[count]] <- c(i-1, j)
  }

  if(i + 1 <= dim.z[1]) {
    count <- count + 1
    neigh.idx[[count]] <- c(i + 1, j)
  }

  if(j - 1 >= 1) {
    count <- count + 1
    neigh.idx[[count]] <- c(i, j - 1)
  }

  if(j + 1 <= dim.z[2]) {
    count <- count + 1
    neigh.idx[[count]] <- c(i, j + 1)
  }

  return(neigh.idx)

}

compute_log_cond_prob <- function(z, i, j, alpha = 0, beta = 0) {
  xi <- z[i,j]
  xi_neighbor <- sapply(find_neighbor(z, i, j), function(nn) {
    z[nn[1], nn[2]]
  })

  beta_coef <- sum(xi == xi_neighbor)
  num <- xi * alpha + beta_coef * beta
  dem <- log(exp(num) + exp( (1-xi)*alpha + (length(xi_neighbor) - beta_coef)*beta ))

  return(num - dem)
}

compute_neg_log_prof_llh <- function(par, z) {
  alpha <- par[1]
  beta <- par[2]

  result <- sapply(1:nrow(z), function(i) {
    tt <- sapply(1:ncol(z), function(j) {
      compute_log_cond_prob(z, i, j, 0, beta)
    })
    sum(tt)
  })

  return(-sum(result))
}

compute_neg_log_prof_llh_beta <- function(par, z) {
  beta <- par[1]

  result <- sapply(1:nrow(z), function(i) {
    tt <- sapply(1:ncol(z), function(j) {
      compute_log_cond_prob(z, i, j, 0, beta)
    })
    sum(tt)
  })

  return(-sum(result))
}

compute_neg_log_prof_llh_alpha <- function(par, z) {
  alpha <- par[1]

  result <- sapply(1:nrow(z), function(i) {
    tt <- sapply(1:ncol(z), function(j) {
      compute_log_cond_prob(z, i, j, alpha, 0)
    })
    sum(tt)
  })

  return(-sum(result))
}


dat = LSBs[[7]]$lsb[, , 1]
opt.r <- optim(par = c(0.5), compute_neg_log_prof_llh_beta, z = dat, method = "L-BFGS-B")

system.time({
  opt.r.alpha <- optim(par = c(0.5),
                       compute_neg_log_prof_llh_alpha, z = dat,
                       method = "L-BFGS-B", lower = 0.000001)
})

alpha.hat <- opt.r.alpha$par
p <- exp(alpha.hat) / (exp(alpha.hat) + 1)


# sample from binary
# smp = matrix(rbinom(nrow(dat)^2, 1, 0.5), nrow = nrow(dat))
# opt.rs <- optim(par = c(0, 0), compute_neg_log_prof_llh, z = smp, method = "L-BFGS-B")


# n <- nrow(dat)
# beta <- 1
#
# mask <- matrix(1,n,n)
# neigh <- getNeighbors(mask, c(2,2,0,0)) # the neighborhood structure# 1st order neighborhood in 2D
# block <- getBlocks(mask, 2)
# k <- 2 #(number of classes, k=2 makes Potts’ to be an Ising model)
# system.time(result <- swNoData(beta = beta,k = k,neigh = neigh, block = block))
# z <- matrix(max.col(result$z)[1:nrow(neigh)], nrow=nrow(mask))
