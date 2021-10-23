library(tidyverse)
library(parallel)
library(bayesImageS)
library(cubelyr)
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



count_neighbor_dl <- function(z, cores = 6) {
  i.max <- dim(z)[1]
  j.max <- dim(z)[2]

  grid <- array(0, dim = c(2, 5, 3),
                dimnames = list(xi=sprintf("%d", 0:1),
                                beta_coef=sprintf("%d", 0:4),
                                type=sprintf("%d", 2:4)))

  for(i in 1:i.max) {
    for(j in 1:j.max)  {
      xi <- z[i,j]
      xi_neighbor <- sapply(find_neighbor(z, i, j), function(nn) { z[nn[1], nn[2]] })

      beta_coef <- sum(xi == xi_neighbor)
      type <- length(xi_neighbor)

      grid[xi+1, beta_coef+1, type-1] <- grid[xi+1, beta_coef+1, type-1] + 1
    }
  }

  result <- grid %>% as.tbl_cube(met_name = "count") %>% as_tibble()
  result <- result %>% filter(count != 0)

  return(result)
}

compute_neg_log_prof_llh.neighbor <- function(par, z.neighbor, alpha_zero = FALSE, beta_zero = FALSE) {
  if(length(par) == 2){
    alpha <- par[1]
    beta <- par[2]
  } else {
    if(alpha_zero & !beta_zero) {
      alpha <- 0
      beta <- par
    } else if (beta_zero & !alpha_zero) {
      alpha <- par
      beta <- 0
    } else{
      alpha <- 0
      beta <- 0
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

compute_MPLE_ratio <- function(dat.neighbor, q) {

  opt.full <- optim(par = c(0.5, 0.5), compute_neg_log_prof_llh.neighbor,
                    z.neighbor = dat.neighbor,
                    method = "L-BFGS-B", lower = c(0.0001, 0.0001))

  if(q == 0){
    null.par <- 0
    null.val <- -sum(dat.neighbor$count) * log(0.5)
  } else if (q == 1) {
    opt.null <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
                      z.neighbor = dat.neighbor, beta_zero = TRUE,
                      method = "L-BFGS-B", lower = 0.0001)
    null.par <- opt.null$par
    null.val <- opt.null$value
  } else {
    opt.null <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
                      z.neighbor = dat.neighbor, alpha_zero = TRUE,
                      method = "L-BFGS-B", lower = 0.0001)
    null.par <- opt.null$par
    null.val <- opt.null$value
  }

  return(list(full.par = opt.full$par, full.val = opt.full$value,
              null.par = null.par, null.val = null.val))
}

fun2 <- function(LSBs, id, layer, B = 50, q = 0, FUN = count_neighbor) {

  dat.original <- LSBs[[id]]$lsb[, , layer]
  dat.original.neighbor <- FUN(dat.original)

  result <- list()
  result[[1]] <- compute_MPLE_ratio(dat.original.neighbor, q)


  if(q == 0) {
    for(b in 1:B) {
      smp = matrix(rbinom(nrow(dat.original)*ncol(dat.original), 1, 0.5), nrow = nrow(dat.original))
      dat.neighbor <- FUN(smp)

      result[[b+1]] <- compute_MPLE_ratio(dat.neighbor, q)
    }

  } else if (q == 1) {
    opt.r <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
                   z.neighbor = dat.original.neighbor, beta_zero = TRUE,
                   method = "L-BFGS-B", lower = 0.0001)
    p = exp(opt.r$par)/(1 + exp(opt.r$par))
    for(b in 1:B) {
      smp = matrix(rbinom(nrow(dat.original)*ncol(dat.original), 1, p), nrow = nrow(dat.original))
      dat.neighbor <- FUN(smp)

      result[[b+1]] <- compute_MPLE_ratio(dat.neighbor, q)
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
    for(b in 1:B){
      result <- swNoData(beta = beta,k = k,neigh = neigh, block = block)
      z <- matrix(max.col(result$z)[1:nrow(neigh)], nrow=nrow(mask))
      smp = z - 1
      dat.neighbor <- FUN(smp)

      result[[b+1]] <- compute_MPLE_ratio(dat.neighbor, q)
    }
  }

  return(result)
}

fun <- function(LSBs, B = 50, id, q = 0, FUN = count_neighbor) {
  res = c()
  res2 = c()
  # foreach (b=1:B) %dopar% {
  for(b in 1:B) {
    par = c()
    val = c()
    for(i in 1:3) {
      dat = LSBs[[id]]$lsb[, , i]
      dat.neighbor <- FUN(dat)
      opt <- optim(par = c(0.5, 0.5), compute_neg_log_prof_llh.neighbor,
                   z.neighbor = dat.neighbor,
                   method = "L-BFGS-B", lower = c(0.0001, 0.0001))
      par = c(par, opt$par)
      val = c(val, opt$value)
      if(q == 0) {
        smp = matrix(rbinom(nrow(dat)*ncol(dat), 1, 0.5), nrow = nrow(dat))
        dat.neighbor <- FUN(smp)
        opt <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
                     z.neighbor = dat.neighbor, alpha_zero = TRUE,
                     method = "L-BFGS-B", lower = 0.0001)
      } else if (q == 1) {
        opt.r <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
                       z.neighbor = dat.neighbor, beta_zero = TRUE,
                       method = "L-BFGS-B", lower = 0.0001)
        p = exp(opt.r$par)/(1 + exp(opt.r$par))
        smp = matrix(rbinom(nrow(dat)*ncol(dat), 1, p), nrow = nrow(dat))
        dat.neighbor <- FUN(smp)
        opt <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
                     z.neighbor = dat.neighbor, beta_zero = TRUE,
                     method = "L-BFGS-B", lower = 0.0001)
      } else {
        n <- nrow(dat)
        m <- ncol(dat)
        opt.r <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
                       z.neighbor = dat.neighbor, alpha_zero = TRUE,
                       method = "L-BFGS-B", lower = 0.0001)
        beta <- opt.r$par
        mask <- matrix(1,n,m)
        neigh <- getNeighbors(mask, c(2,2,0,0))
        block <- getBlocks(mask, 2)
        k <- 2
        result <- swNoData(beta = beta,k = k,neigh = neigh, block = block)
        z <- matrix(max.col(result$z)[1:nrow(neigh)], nrow=nrow(mask))
        smp = z - 1
        dat.neighbor <- FUN(smp)
        opt <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
                     z.neighbor = dat.neighbor, alpha_zero = TRUE,
                     method = "L-BFGS-B", lower = 0.0001)
      }
      par = c(par, opt$par)
      val = c(val, opt$value)
    }
    res = rbind(res, par)
    res2 = rbind(res2, val)
  }
  return(list(res, res2))
}

### When beta = 0
system.time({
  res = fun(B = 1, id = 7, q = 1, LSBs = LSBs, FUN = count_neighbor_mc)
})

system.time({
  # res2 = lapply(1:3, function(layer) {
  res2 <- fun2(B = 1, id = 5, q = 0, layer = 1, LSBs = LSBs, FUN = count_neighbor_mc)
  # })
})


system.time({
  res = fun(B = 1, id = 7, q = 1, LSBs = LSBs)
})

### When alpha = 0
res = fun(B = 1, id = 7, q = 2, LSBs = LSBs)
#### code for testing
# n <- 50
# beta <- 1
# mask <- matrix(1,n,n) # basically the grid
# neigh <- getNeighbors(mask, c(2,2,0,0)) # the neighborhood structure
# # 1st order neighborhood in 2D
# block <- getBlocks(mask, 2)
# k <- 2 #(number of classes, k=2 makes Potts’ to be an Ising model)
# result <- swNoData(beta = beta,k = k,neigh = neigh, block = block)
# z <- matrix(max.col(result$z)[1:nrow(neigh)], nrow=nrow(mask))
# z <- z - 1

# z.neighbor <- count_neighbor(z)
# par <- c(0.5, 1.5)

# compute_neg_log_prof_llh(par, z)
# compute_neg_log_prof_llh.neighbor(par, z.neighbor)

# # run on the data
dat = LSBs[[7]]$lsb[, , 1]
dat.neighbor <- count_neighbor_mc(dat)
opt.r.beta <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor, z.neighbor = dat.neighbor,
                    alpha_zero = TRUE, beta_zero = TRUE,
                    method = "L-BFGS-B", lower = c(0.0001, 0.0001))

# system.time({
#   opt.r.beta <- optim(par = c(0.5), compute_neg_log_prof_llh, z = dat, alpha_zero = TRUE,
#                       method = "L-BFGS-B")
# })
# # user  system elapsed
# # 252.109   1.147 254.841
# opt.r.beta$par
# opt.r.beta$value

# system.time({
#   dat.neighbor <- count_neighbor(dat)
#   opt.r.beta.neighbor <- optim(par = c(0.5), compute_neg_log_prof_llh.neighbor,
#                                z.neighbor = dat.neighbor, alpha_zero = TRUE,
#                                method = "L-BFGS-B")
# })
# user  system elapsed
# 15.058   0.065  15.133
# opt.r.beta.neighbor$par
# opt.r.beta.neighbor$value

# alpha.hat <- opt.r.alpha$par
# p <- exp(alpha.hat) / (exp(alpha.hat) + 1)

rbenchmark::benchmark(
  neighbor1 = {
    dat.neighbor1 <- count_neighbor(dat)
  },
  neighbor_mc = {
    dat.neighbor2 <- count_neighbor_mc(dat, cores = 6)
  },
  neighbor_dl = {
    dat.neighbor3 <- count_neighbor_dl(dat)
  },
  replications = 1,
  columns = c("test", "replications", "elapsed", "relative")
)

# sample from binary
# smp = matrix(rbinom(nrow(dat)^2, 1, 0.5), nrow = nrow(dat))
# opt.rs <- optim(par = c(0, 0), compute_neg_log_prof_llh, z = smp, method = "L-BFGS-B")

nrow(dat)
n <- 50
beta <- 1

mask <- matrix(1,n,n)
neigh <- getNeighbors(mask, c(2,2,0,0)) # the neighborhood structure# 1st order neighborhood in 2D
block <- getBlocks(mask, 2)
k <- 2 #(number of classes, k=2 makes Potts’ to be an Ising model)
system.time(result <- swNoData(beta = beta,k = k,neigh = neigh, block = block))
z2 <- matrix(max.col(result$z)[1:nrow(neigh)], nrow=nrow(mask))
z <- z2 - 1
