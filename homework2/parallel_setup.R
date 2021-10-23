library(tidyverse)
library(parallel)
library(doParallel)
library(bayesImageS)
library(cubelyr)


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


count_neighbor_dl <- function(z) {
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


fun3 <- function(LSBs, id, layer, B = 50, q = 0, FUN = count_neighbor) {

  dat.original <- LSBs[[id]]$lsb[, , layer]
  dat.original.neighbor <- FUN(dat.original)

  result <- list()
  result[[1]] <- compute_MPLE_ratio(dat.original.neighbor, q)

  if(q == 0) {
    tmp <- foreach(b=1:B) %dopar% {
      source("homework2/parallel_setup.R")
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
    tmp <- foreach(b = 1:B) %dopar% {
      source("homework2/parallel_setup.R")
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
    tmp <- foreach(b = 1:B) %dopar% {
      source("homework2/parallel_setup.R")
      result <- swNoData(beta = beta,k = k,neigh = neigh, block = block)
      z <- matrix(max.col(result$z)[1:nrow(neigh)], nrow=nrow(mask))
      smp = z - 1
      dat.neighbor <- FUN(smp)

      compute_MPLE_ratio(dat.neighbor, q)
      # result[[b+1]] <- compute_MPLE_ratio(dat.neighbor, q)
    }
  }
  result <- c(result, tmp)

  return(result)
}
