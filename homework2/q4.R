library(tidyverse)

library(bayesImageS)
lsb = load("../STAT680/local/LSB.rda")
# lsb = load("./LSB.rda")

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

compute_neg_log_prof_llh <- function(par, z, alpha_zero = FALSE, beta_zero = FALSE) {
  if(length(par) == 2){
    alpha <- par[1]
    beta <- par[2]
  } else {
    if(alpha_zero) {
      alpha <- 0
      beta <- par
    } else if (beta_zero) {
      alpha <- par
      beta <- 0
    } else{
      stop("both parameters are 0")
    }
  }

  result <- sapply(1:nrow(z), function(i) {
    tt <- sapply(1:ncol(z), function(j) {
      compute_log_cond_prob(z, i, j, alpha, beta)
    })
    sum(tt)
  })

  return(-sum(result))
}

# compute_neg_log_prof_llh_beta <- function(par, z) {
#   beta <- par[1]
#
#   result <- sapply(1:nrow(z), function(i) {
#     tt <- sapply(1:ncol(z), function(j) {
#       compute_log_cond_prob(z, i, j, 0, beta)
#     })
#     sum(tt)
#   })
#
#   return(-sum(result))
# }
#
# compute_neg_log_prof_llh_alpha <- function(par, z) {
#   alpha <- par[1]
#
#   result <- sapply(1:nrow(z), function(i) {
#     tt <- sapply(1:ncol(z), function(j) {
#       compute_log_cond_prob(z, i, j, alpha, 0)
#     })
#     sum(tt)
#   })
#
#   return(-sum(result))
# }

count_neighbor <- function(z) {
  grid <- expand.grid(i = 1:dim(z)[1], j = 1:dim(z)[2])

  grid <- grid %>% mutate(
    xi = purrr::map2_dbl(i, j, function(i, j) z[i,j]),
    xi_neighbor = purrr::map2(i, j, function(i, j) {
      sapply(find_neighbor(z, i, j), function(nn) {
        z[nn[1], nn[2]]
      })
    }),
    beta_coef = purrr::map2_dbl(xi, xi_neighbor, function(xi, xi_neighbor) {
      sum(xi == xi_neighbor)
    }),
    type = purrr::map_dbl(xi_neighbor, function(xi_neighbor) {length(xi_neighbor)})
  )

  result <- grid %>% select(-i, -j, -xi_neighbor) %>% group_by(xi, beta_coef, type) %>%
    summarise(
      count = n(),
      .groups = "drop"
    )

  return(result)
}

compute_neg_log_prof_llh.neighbor <- function(par, z.neighbor, alpha_zero = FALSE, beta_zero = FALSE) {
  if(length(par) == 2){
    alpha <- par[1]
    beta <- par[2]
  } else {
    if(alpha_zero) {
      alpha <- 0
      beta <- par
    } else if (beta_zero) {
      alpha <- par
      beta <- 0
    } else{
      stop("both parameters are 0")
    }
  }

  z.neighbor <- z.neighbor %>% mutate(
    log_cond_prob = purrr::pmap_dbl(list(xi, beta_coef, type), function(xi, beta_coef, type) {
      num <- xi * alpha + beta_coef * beta
      denom <- log(exp(num) + exp( (1-xi) * alpha + (type - beta_coef) * beta ))
      num - denom
    })
  )

  result <- with(z.neighbor, {
    -sum(count * log_cond_prob)
  })

  return(result)
}

#### code for testing
n <- 50
beta <- 1
mask <- matrix(1,n,n) # basically the grid
neigh <- getNeighbors(mask, c(2,2,0,0)) # the neighborhood structure
# 1st order neighborhood in 2D
block <- getBlocks(mask, 2)
k <- 2 #(number of classes, k=2 makes Potts’ to be an Ising model)
result <- swNoData(beta = beta,k = k,neigh = neigh, block = block)
z <- matrix(max.col(result$z)[1:nrow(neigh)], nrow=nrow(mask))
z <- z - 1

z.neighbor <- count_neighbor(z)
par <- c(0.5, 1.5)

compute_neg_log_prof_llh(par, z)
compute_neg_log_prof_llh.neighbor(par, z.neighbor)

# run on the data
dat = LSBs[[7]]$lsb[, , 1]

system.time({
  opt.r.beta <- optim(par = c(0.5), compute_neg_log_prof_llh, z = dat, alpha_zero = TRUE,
                      method = "L-BFGS-B")
})
# user  system elapsed
# 252.109   1.147 254.841
opt.r.beta$par
opt.r.beta$value

system.time({
  dat.neighbor <- count_neighbor(dat)
  opt.r.beta.neighbor <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
                               z.neighbor = dat.neighbor, alpha_zero = TRUE,
                               method = "L-BFGS-B")
})
# user  system elapsed
# 15.058   0.065  15.133
opt.r.beta.neighbor$par
opt.r.beta.neighbor$value

# alpha.hat <- opt.r.alpha$par
# p <- exp(alpha.hat) / (exp(alpha.hat) + 1)



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
