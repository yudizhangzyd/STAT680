library(osDesign)
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

x_prob_vec <- Vectorize(x_prob, vectorize.args = c("i", "j"))
rbinom.vec <- Vectorize(rbinom, vectorize.args = "prob")

sampler.step <- function(z, par) {
  dim.z <- dim(z)
  nr <- dim.z[1]
  nc <- dim.z[2]

  grid <- expand.grid(i = 1:nr, j = 1:nc)
  grid$type <- abs(grid$i - grid$j) %% 2
  grid$rowidx <- as.numeric(rownames(grid))

  gp.first <- grid %>% filter(type == 0)
  gp.second <- grid %>% filter(type == 1)

  z.next <- rep(NA, nr*nc)

  if( runif(1) > 0.5) {
    gp1 <- gp.first
    gp2 <- gp.second
  } else {
    gp1 <- gp.second
    gp2 <- gp.first
  }

  prob.gp1 <- x_prob_vec(gp1$i, gp1$j, z, par)

  z.next[gp1$rowidx] <- rbinom.vec(1, 1, prob.gp1)
  z.next2 <- matrix(z.next, nrow = nr, ncol = nc)

  prob.gp2 <- x_prob_vec(gp2$i, gp2$j, z.next2, par)
  z.next[gp2$rowidx] <- rbinom.vec(1, 1, prob.gp2)

  z.next <- matrix(z.next, nrow = nr, ncol = nc)

  return(z.next)

}

z <- matrix(rbinom(16*16, 1, 0.5), ncol = 16)
par <- c(0.25, 0.5, 0.5)

z <- matrix(1, nrow = 16, ncol=16)
z <- sampler.step(z, par)
z






