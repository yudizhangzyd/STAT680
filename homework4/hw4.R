library(osDesign)
library(parallel)
library(doParallel)

k <- c(4,3,2,4)
x <- c(1,0,1,4)
n <- sum(x)
all <- sum(k)
mcp <- function(B){
  pval <- prod(choose(k,x))/choose(all,n)
  iter <- 0
  while (iter < B) {
    iter <- iter + 1
    x <- rmvhyper(k,n)
    u <- prod(choose(k,x))/choose(all,n)
    pval <- c(pval,u)
  }
  ranks <- rank(pval, ties.method = "random")
  1 - (B + 2 - ranks[1])/(B + 1)
}
mcp(999)

# question 2

x_neighbor <- function(i, j, dim.z, type = "h") {
  if(i > dim.z[1] | j > dim.z[2]) return(NULL)

  neigh.idx <- list()
  count <- 0

  if(type == "h") {
    if(j - 1 >= 1) {
      count <- count + 1
      neigh.idx[[count]] <- c(i, j - 1)
    }

    if(j + 1 <= dim.z[2]) {
      count <- count + 1
      neigh.idx[[count]] <- c(i, j + 1)
    }

    return(neigh.idx)
  } else {
    if(i - 1 >= 1) {
      count <- count + 1
      neigh.idx[[count]] <- c(i-1, j)
    }

    if(i + 1 <= dim.z[1]) {
      count <- count + 1
      neigh.idx[[count]] <- c(i + 1, j)
    }

    return(neigh.idx)
  }

}

x_prob <- function(i, j, z, par) {
  alpha <- par[1]
  beta.h <- par[2]
  beta.v <- par[3]

  dim.z <- dim(z)

  neighbor.h <- sapply(x_neighbor(i,j, dim.z, "h"),
                       function(idx) {
                         z[idx[1], idx[2]]
                       })
  neighbor.v <- sapply(x_neighbor(i,j, dim.z, "v"),
                       function(idx) {
                         z[idx[1], idx[2]]
                       })

  neighbor.one.h <- neighbor.h == 1
  neighbor.one.v <- neighbor.v == 1

  one.prob <- exp(
    alpha + sum(neighbor.one.h) * beta.h + sum(neighbor.one.v) * beta.v
  )

  zero.prob <- exp(
    beta.h * sum(!neighbor.one.h) + beta.v * sum(!neighbor.one.v)
  )

  return(one.prob / (one.prob + zero.prob))

}

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

x_prob_vec <- Vectorize(x_prob, vectorize.args = c("i", "j"))
rbinom.vec <- Vectorize(rbinom, vectorize.args = "prob")

sampler.step <- function(z, par, n.core) {
  dim.z <- dim(z)
  nr <- dim.z[1]
  nc <- dim.z[2]

  grid <- expand.grid(i = 1:nr, j = 1:nc)
  grid$type <- abs(grid$i - grid$j) %% 2
  grid$rowidx <- as.numeric(rownames(grid))

  gp.first <- grid[grid$type == 0,]
  gp.second <- grid[grid$type == 1,]

  z.next <- rep(NA, nr*nc)

  if( runif(1) > 0.5) {
    gp1 <- gp.first
    gp2 <- gp.second
  } else {
    gp1 <- gp.second
    gp2 <- gp.first
  }

  # prob.gp1 <- x_prob_vec(gp1$i, gp1$j, z, par)
  prob.gp1 <- x_prob_parallel(gp1$i, gp1$j, z, par, n.core = n.core)

  z.next[gp1$rowidx] <- rbinom.vec(1, 1, prob.gp1)
  z.next2 <- matrix(z.next, nrow = nr, ncol = nc)

  # prob.gp2 <- x_prob_vec(gp2$i, gp2$j, z.next2, par)
  prob.gp2 <- x_prob_parallel(gp2$i, gp2$j, z.next2, par, n.core = n.core)
  z.next[gp2$rowidx] <- rbinom.vec(1, 1, prob.gp2)

  z.next <- matrix(z.next, nrow = nr, ncol = nc)

  return(z.next)

}

gibbs.sampler.ising <- function(z, par, n.core, sample.size = 100) {
  result <- list()
  for(i in 1:sample.size) {
    result[[i]] <- sampler.step(z, par, n.core)
    z <- result[[i]]
  }

  return(result)
}

z <- matrix(rbinom(16*16, 1, 0.5), ncol = 16)
par <- c(0.25, 0.5, 0.5)

z <- matrix(1, nrow = 16, ncol=16)
z <- sampler.step(z, par)
z.list <- gibbs.sampler.ising(z, par, 10)

#### (b) perfect simulation
find.perfect.m <- function(m, par, dim.z, n.core, seed = 17381128) {
  nr <- dim.z[1]
  nc <- dim.z[2]
  zeros <- matrix(0, nrow = nr, ncol = nc)
  ones <- matrix(1, nrow = nr, ncol = nc)

  set.seed(seed)
  one.list <- gibbs.sampler.ising(ones, par, n.core = n.core, sample.size = m)

  set.seed(seed)
  zero.list <- gibbs.sampler.ising(zeros, par, n.core = n.core, sample.size = m)

  if(all(one.list[[m]] == zero.list[[m]])) {
    return(list(m = m, x.zero = one.list[[m]]))
  } else {
    return(NULL)
  }

}


find.perfect.sample <- function(dim.z, par, n.core, seed = 17381128){
  m <- 0
  while(TRUE) {
    m <- m + 1
    result <- find.perfect.m(m, par, dim.z, n.core, seed = seed)
    if(!is.null(result)) {
      break
    }
  }
  return(result)
}


result <- find.perfect.sample(c(16, 16), par, n.core = n.core, m.start = 26)
result$m

m <- result$m
set.seed(17381128)
t1 <- gibbs.sampler.ising(ones, par, sample.size = m)[[m]]
set.seed(17381128)
t2 <- gibbs.sampler.ising(zeros, par, sample.size = m)[[m]]

all(t1 == t2)
all(t1 == result$x.zero)

# result$x.zero is the perfect sample

set.seed(17381128)
t3 <- gibbs.sampler.ising(zeros, par, sample.size = m+1)
all(t3[[m]] == t1)

### (c)
dim.z <- c(512, 512)
par <- c(0, 1, 1)

cl <- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)

nr <- dim.z[1]
nc <- dim.z[2]
zeros <- matrix(0, nrow = nr, ncol = nc)
ones <- matrix(1, nrow = nr, ncol = nc)

set.seed(seed)
one.list <- gibbs.sampler.ising(ones, par, sample.size = 1)

set.seed(seed)
zero.list <- gibbs.sampler.ising(zeros, par, sample.size = m)

Rprof()
find.perfect.m(2, par, dim.z)
Rprof(NULL)

parallel::stopCluster(cl)

result.c <- find.perfect.sample(dim.z, par, seed = 18271128)
saveRDS(result.c, file = "homework4/c_perfect_sample_512.rds")


!is.null(find.perfect.m(100, par, dim.z))


