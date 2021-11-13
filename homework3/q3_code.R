library(plot3D)

library(ars)
library(armspp)

d = oro.nifti::readNIfTI(fname = "local/new_phantom.nii.gz")
data = d@.Data
data = drop(data)
r1 = data[40:53,95:114 , ]
image2D(r1)
r2 = data[81:90, 81:90, ]
image2D(r2)
r3 = data[100:167, 100:167, ]
image2D(r3)

f<-function(x,mu=0,sigma=1){-1/(2*sigma^2)*(x-mu)^2}
fprima<-function(x,mu=0,sigma=1){-1/sigma^2*(x-mu)}
mysample<-ars(200,f,fprima,mu=0,sigma=1)
mysample
hist(mysample)


logf <- function(x, r, rho, sigma1, sigma2, tau) {
  -1/(2*(1-rho^2)) * (r - tau)^2 *
    (cos(x)^2/sigma1^2 - (2*rho*cos(x)*sin(x))/(sigma1*sigma2) + sin(x)^2/sigma2^2)
}

logfprime <- function(x, r, rho, sigma1, sigma2, tau) {
  -1/(2*(1-rho^2)) * (r - tau)^2 *
    (
      (-2*sin(x)*cos(x))/sigma1^2 -
        (2*rho*cos(2*x))/(sigma1 * sigma2) +
        (2*sin(x)*cos(x))/sigma2^2
    )
}

mysample<-ars(200, logf, logfprime, x=0, m=1, xlb=-1, xub=1, lb=TRUE, ub=TRUE,
              r=r, rho=rho, sigma1=sigma1, sigma2=sigma2, tau=tau)
mysample
hist(mysample)

r <- 1300
rho <- 0.4
sigma1 <- 30
sigma2 <- 30
tau <- 1499

x <- seq(-pi, pi, by=0.01)
y <- logf(x, r=r, rho = rho, sigma1 = sigma1, sigma2 = sigma2, tau = tau)
plot(x,exp(y))

mysample <- arms(10000, logf, -pi, pi, metropolis = TRUE,
                 arguments = list(r=r, rho=rho, sigma1=sigma1, sigma2=sigma2, tau=tau))
hist(mysample)

# par c(sigma1, sigma2, rho, tau)



Qfunc.each <- function(par, par.cur, r, sample.size = 10000, FUNC = logf){
  theta.cond <- arms(sample.size, FUNC, -pi, pi, metropolis = TRUE,
                      arguments = list(
                        r=r, sigma1=par.cur[1], sigma2=par.cur[2], rho=par.cur[3], tau=par.cur[4]
                      ))
  sigma1 <- par[1]
  sigma2 <- par[2]
  rho <- par[3]
  tau <- par[4]

  intt <-
    cos(theta.cond)^2/sigma1^2 -
    (2*rho*sin(theta.cond)*cos(theta.cond))/(sigma1*sigma2) +
    sin(theta.cond)^2/sigma2^2
  intt <- mean(intt)

  result <-
    -log(sigma1) - log(sigma2) - 0.5*log(1-rho^2) + log(r) -
    1/(2*(1-rho^2)) * (r-tau)^2 * intt

  result
}

Qfunc <- function(par, par.cur, r.vec, sample.size=1000) {
  result <- sapply(r.vec, function(r) {
    Qfunc.each(par, par.cur, r, sample.size)
  })

  return(-sum(result))
}

Qfunc(par.cur, par.cur, c(r2[,,1]), sample.size = 1000)

rr <- c(r2[,,1])
par.cur <- c(sd(rr), sd(rr), 0.5, mean(rr))

opt <- optim(par.cur, Qfunc, par.cur=par.cur, r.vec=rr, sample.size=1000,
             method = "L-BFGS-B", lower = c(0.0001,0.0001,-0.99,0))

em_q3 <- function(r, iter.max=40, tol=1e-4, sample.size=1000){
  par.cur <- c(sd(r), sd(r), 0.5, mean(r))

  iter <- 0
  while(TRUE){
    iter <- iter + 1

    par.old <- par.cur
    # obtain the estimates
    opt <- optim(par.cur, Qfunc, par.cur=par.cur, r.vec=r, sample.size=sample.size,
                 method = "L-BFGS-B", lower = c(0.0001,0.0001,-0.99,0.0001),
                 upper=c(1000, 1000, 0.99, 5000))

    # update the current value
    par.cur <- opt$par

    l2.diff <- sum((par.cur-par.old)^2)

    cat(paste("iter:", iter, "l2.diff:", l2.diff, "Q value:", opt$value, "\n"))

    if(l2.diff < tol){
      cat("less than the tolerance. exit\n")
      break
    } else if (iter > iter.max) {
      cat("reach max iteration. exit\n")
      break
    }
  }

  return(par.cur)

}

em_q3(rr, tol = 1e-9, sample.size = 100)


selected.region <- list(r1, r2, r3)
save(selected.region, file="local/selected_region.rda")


ests <- list()
for(i in 1:18) {
  ests[[i]] <- list()
  for (j in 1:3) {
    ests[[i]][[j]] <- readRDS(
      paste("homework3/results_q3/result-", j, "-", i, "-est.rds", sep = "")
      )
  }
}

variance.est.each <- function(par, r, sample.size = 1000, FUNC = logf){
  sigma1 <- par[1]
  sigma2 <- par[2]
  rho <- par[3]
  tau <- par[4]

  theta.cond <- arms(sample.size, FUNC, -pi, pi, metropolis = TRUE,
                     arguments = list(
                       r=r, sigma1=sigma1, sigma2=sigma2, rho=rho, tau=tau
                     ))

  part.sigma1 <- mean(-2*cos(theta.cond)^2/sigma1^3 + (2*rho*sin(theta.cond)*cos(theta.cond))/(sigma1^2*sigma2))
  part.sigma2 <- mean(-2*sin(theta.cond)^2/sigma2^3 + (2*rho*sin(theta.cond)*cos(theta.cond))/(sigma1*sigma2^2))
  part.rho1 <- mean(cos(theta.cond)^2/sigma1^2 - (2*rho*cos(theta.cond)*sin(theta.cond))/(sigma1*sigma2) +
    sin(theta.cond)^2/sigma2^2)
  part.rho2 <- mean(-2*sin(theta.cond)*cos(theta.cond)/(sigma1*sigma2))
  part.tau <- part.rho1

  dev.sigma1 <- -1/sigma1 - 1/(2*(1-rho^2))*(r-tau)^2*part.sigma1
  dev.sigma2 <- -1/sigma2 - 1/(2*(1-rho^2))*(r-tau)^2*part.sigma2
  dev.rho <- rho/(1-rho^2) + (r-tau)^2 * (-rho/(1-rho^2)^2 * part.rho1 - 1/(2*(1-rho^2)) * part.rho2)
  dev.tau <- (r-tau)/(1-rho^2) * part.tau

  q <- c(dev.sigma1, dev.sigma2, dev.rho, dev.tau)

  I.matrix <- q %*% t(q)

  return(I.matrix)
}

load("~/STAT_680/STAT680/local/selected_region.rda")

I.est <- variance.est.each(ests[[1]][[2]], selected.region[[2]][,,1][1])


variance.est <- function(par, r.vec, sample.size = 1000, return_inv = TRUE){
  result <- lapply(r.vec, function(r) {
    variance.est.each(par, r, sample.size)
  })

  result <- Reduce("+", result)

  if(return_inv) {
    return(solve(result))
  } else {
    return(result)
  }

}

ests.var <- list()
for(i in 1:18) {
  ests.var[[i]] <- list()
  for(j in 1:3) {
    ests.var[[i]][[j]] <- variance.est(par=ests[[i]][[j]], r.vec=c(selected.region[[j]][,,i]) )
    cat("finish for i:", i, "j:", j, "\n")
  }
}

library(tidyverse)
result.container <- expand.grid(region=1:3, image=1:18) %>% as_tibble()

result.container <- result.container %>% mutate(
  estimates = purrr::map2(region, image, function(j, i) {
    ests[[i]][[j]]
  }),
  variance.matrix = purrr::map2(region, image, function(j, i) {
    ests.var[[i]][[j]]
  }),
  variance = purrr::map(variance.matrix, diag)
)

result.container %>% mutate(
  tau.value = purrr::map_dbl(estimates, function(est) est[4]),
  sigma.ratio = purrr::map_dbl(estimates, function(est) est[1] / est[2])
) %>% pivot_longer(cols = c(tau.value, sigma.ratio), names_to = "type", values_to = "value") %>%
  ggplot() +
  geom_line(aes(x=image, y=value, color=factor(region))) +
  facet_wrap(~type, scales = "free")


result.container.all <- result.container %>% unnest(c(estimates, variance))
result.container.all$parameter <- rep(c("sigma1", "sigma2", "rho", "tau"), 18 * 3)

result.container.all <- result.container.all %>% mutate(
  lower = estimates - qnorm(.975) * sqrt(variance),
  upper = estimates + qnorm(.975) * sqrt(variance)
)

result.container %>% mutate(
  var.sig1_sig2 = purrr::map_dbl(variance.matrix, function(vv) {
    t(c(1,-1,0,0)) %*% vv %*% c(1,-1,0,0)
  }),
  sig1_sig2 = purrr::map_dbl(estimates, function(est) est[1] - est[2]),
  lower = sig1_sig2 - qnorm(.975) * sqrt(var.sig1_sig2),
  upper = sig1_sig2 + qnorm(.975) * sqrt(var.sig1_sig2)
)


save(result.container, ests, ests.var, file="homework3/q3_all_results.rda")







