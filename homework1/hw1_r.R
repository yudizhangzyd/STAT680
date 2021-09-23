library(tidyverse)
library(magrittr, include.only = c("%<>%", "%$%"))

# Q1

## a)
n <- 1000
u <- runif(n)
ee <- -log(1-u)

alpha <- 1.89

tt <- ee^(alpha - 1)
var(tt) / n

(mean(ee^(2*alpha-2)) - (mean(tt))^2) / n

e.gamma <- ee^(alpha - 1) %>% mean
e.gamma
gamma(alpha)

## b)

## c)
estimate_gamma <- function(alpha, n = 1000) {
  u <- runif(n)
  ee <- -log(1-u)

  tt <- ee^(alpha - 1)
  est <- mean(tt)
  sd <- sqrt(var(tt) / n)
  sd2 <- sqrt((mean(ee^(2*alpha-2)) - (mean(tt))^2) / n)

  result <- c(est, sd, sd2)
  names(result) <- c("estimate", "sd", "sd2")
  return(result)
}

alpha <- 1.5
n <- 1000
estimate_gamma(alpha, n)
gamma(alpha)

## d)
estimate_gamma_antithetic <- function(alpha, n = 1000) {
  u <- runif(n)
  u2 <- 1 - u
  u <- c(u, u2)
  ee <- -log(1-u)

  tt <- ee^(alpha - 1)
  est <- mean(tt)
  sd <- sqrt(var(tt) / (2*n))
  sd2 <- sqrt((mean(ee^(2*alpha-2)) - (mean(tt))^2) / (2*n))

  result <- c(est, sd, sd2)
  names(result) <- c("estimate", "sd", "sd2")
  return(result)
}

alpha <- 1.5
n <- 10000
estimate_gamma(alpha, n)
estimate_gamma_antithetic(alpha, n)
gamma(alpha)

## e)
tr_val <- 1 - exp(-1)
n <- 100000

### i)
u <- runif(n)
tt <- as.numeric(u < exp(-(1:n) / n))
mean(tt)
sqrt(var(tt) / n)

### ii)
u <- runif(n)
tt <- exp(-u)
mean(tt)
sqrt(var(tt) / n)

### iii)
u <- runif(n)
u2 <- 1 - u
u <- c(u, u2)
tt <- exp(-u)
mean(tt)
sqrt(var(tt) / (2*n))

### iv)
is_beta <- function(a=2, b=2, n=1000) {
  bb <- rbeta(n, a, b)
  tt <- exp(-bb) / dbeta(bb, a, b)
  est <- mean(tt)
  sd <- sqrt(var(tt) / n)
  result <- c(est, sd)
  names(result) <- c("estimate", "sd")
  return(result)
}

grid <- expand.grid(a=seq(0.2,5,by=.2), b=seq(0.2,5,by=.2))

grid %<>% as_tibble %>% mutate(
  is_beta_result = purrr::map2(a, b, function(aa, bb) {
    rr <- is_beta(aa, bb, n = 100000) %>% t %>% as_tibble
    colnames(rr) <- c("estimate", "sd")
    rr
  })
) %>% unnest(is_beta_result) %>%
  mutate(
    diff = abs(tr_val - estimate)
  )

grid$sd %>% which.min() %>% {grid[.,]}
grid$diff %>% which.min() %>% {grid[.,]}
grid %>% filter(a==1, b==1)
is_beta(1, 1, n = 100000)

saveRDS(grid, "homework1/grid.rds")


# Q2

## a)
library(boot)
require(fitdistrplus)
data("aircondit")

fit.gamma <- fitdist(aircondit$hours, distr = "gamma", method = "mle")
summary(fit.gamma)

plot(fit.gamma)

tb.find <- function(alpha, data) {
  xbar <- mean(data)
  result <- log(alpha) - log(xbar) - digamma(alpha) + mean(log(data))
  return(result)
}

alpha.est.result <- uniroot(tb.find, interval = c(0.001, 4), data = aircondit$hours)
alpha.est.result
alpha.est <- alpha.est.result$root
beta.est <- alpha.est / mean(aircondit$hours)
c(alpha.est, beta.est)
# https://www.math.arizona.edu/~jwatkins/O3_mle.pdf

## b)
