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
