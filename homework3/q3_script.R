library(parallel)
library(doParallel)
library(armspp)

load("data/selected_region.rda")

logf <- function(x, r, rho, sigma1, sigma2, tau) {
  -1/(2*(1-rho^2)) * (r - tau)^2 *
    (cos(x)^2/sigma1^2 - (2*rho*cos(x)*sin(x))/(sigma1*sigma2) + sin(x)^2/sigma2^2)
}

# par c(sigma1, sigma2, rho, tau)
Qfunc.each <- function(par, par.cur, r, sample.size = 10000, FUNC = logf){
  theta.cond <- armspp::arms(sample.size, FUNC, -pi, pi, metropolis = TRUE,
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


cl <- parallel::makeCluster(9)
doParallel::registerDoParallel(cl)

args = commandArgs(trailingOnly=TRUE)

start.idx <- as.numeric(args[1])
end.idx <- as.numeric(args[2])

print("hello1")

result <- foreach(b=start.idx:end.idx) %dopar% {

  for(i in 1:3) {
    rr <- c(selected.region[[i]][,,b])

    print(paste("compute region", i, "for image", b, "\n"))
    result <- em_q3(rr, tol = 1e-8, sample.size = 100)

    saveRDS(result, paste("results_q3/result", i, b, "est.rds", sep="-" ))
  }

}

parallel::stopCluster(cl)

