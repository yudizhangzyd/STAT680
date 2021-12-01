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


z <- matrix(rbinom(16*16, 1, 0.5), ncol = 16)
par <- c(0.25, 0.5, 0.5)

z <- matrix(1, nrow = 16, ncol=16)
z <- sampler.step(z, par, n.core = 6)
z.list <- gibbs.sampler.ising(z, par, n.core = 6, sample.size = 10)


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
source("funcs.R")

dim.z <- c(16,16)
par <- c(0, 1, 1)

n.core <- 6

cl <- parallel::makeCluster(n.core)
doParallel::registerDoParallel(cl)

nr <- dim.z[1]
nc <- dim.z[2]
zeros <- matrix(0, nrow = nr, ncol = nc)
ones <- matrix(1, nrow = nr, ncol = nc)


seed <- 17381128
set.seed(seed)
system.time({
  one.list <- gibbs.sampler.ising(ones, par, n.core = 6, sample.size = 5000)
})
zero.list <- gibbs.sampler.ising(zeros, par, n.core = 6, sample.size = 1000)

identical(one.list[[1000]], zero.list[[1000]])
one.list[[1000]]
zero.list[[1000]]


parallel::stopCluster(cl)

result.c <- find.perfect.sample(dim.z, par, n.core = n.core,
                                m.start = 0, m.by = 5000)
saveRDS(result.c, file = "homework4/c_perfect_sample_512.rds")


!is.null(find.perfect.m(100, par, dim.z))


z <- matrix(rbinom(16*16, 1, 0.5), ncol = 16)



rr <- multi_grid_method(result.c$x.zero, par, n.core=4, sample.size = 100)


example.16 <- perfect.gibbs.multi.grid.sampler(c(16, 16), c(0.25, 0.5, 0.5), n.core = 4, m.by = 10)

source("IsingMPLE.R")
mple.RES <- get.MPLE(example.16[[100]], homogenous = FALSE)
mple.RES





