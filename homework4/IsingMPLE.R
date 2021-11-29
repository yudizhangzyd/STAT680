convert_table2array <- function(tabl, edge.type = "nonedge"){
    dimarr <- switch(edge.type,
                     nonedge = c(2, 3, 3),
                     inneredge.v = c(2, 3, 2),
                     inneredge.h = c(2, 2, 3),
                     edge = c(2, 2, 2))
    if (all(dim(tabl) == dimarr))
        tabl else {
                  arr <- array(0, dim = dimarr)
                  dims <- dim(tabl)
                  dimns <- dimnames(tabl)
                  for (i in 1:dims[1])
                      for (j in 1:dims[[2]])
                          for (k in 1:dims[[3]])
                              arr[as.integer(dimns[[1]][i]) + 1,
                                  as.integer(dimns[[2]][j]) + 1,
                                  as.integer(dimns[[3]][k]) + 1] <- tabl[i,j,k]
                  arr
             }
}

count_all_1stordr_nhbrs  <- function(x) {
    likes.x <- 1*(x[-nrow(x),] == x[-1,])
    nonedge.likes.x <- (likes.x[-nrow(likes.x),]+likes.x[-1,])
    edge.likes.x <- rbind(likes.x[1,], likes.x[nrow(likes.x),])
    likes.y <- 1*(x[,-ncol(x)] == x[,-1])
    nonedge.likes.y <- (likes.y[,-ncol(likes.y)]+likes.y[,-1])
    edge.likes.y <- cbind(likes.y[,1], likes.y[,ncol(likes.y)])
    tmp <- table(x[-c(1,nrow(x)),-c(1,ncol(x))], nonedge.likes.x[,-c(1,ncol(likes.x))], nonedge.likes.y[-c(1, nrow(likes.y)),])
    num.nonedge.combos <- convert_table2array(tmp)
    tmp  <- table(x[c(1,nrow(x)),-c(1,ncol(x))], edge.likes.x[,-c(1, ncol(edge.likes.x))],nonedge.likes.y[c(1,nrow(nonedge.likes.y)),])
    num.inneredge.combos.x <- convert_table2array(tmp, edge.type = "inneredge.h")
    tmp <- table(x[-c(1,nrow(x)),c(1,ncol(x))], nonedge.likes.x[,c(1,ncol(nonedge.likes.x))], edge.likes.y[-c(1, nrow(edge.likes.y)),])
    num.inneredge.combos.y <- convert_table2array(tmp, edge.type = "inneredge.v")
    tmp <- table(x[c(1,nrow(x)),c(1,ncol(x))], edge.likes.x[,c(1,ncol(edge.likes.x))], edge.likes.y[c(1,nrow(edge.likes.y)),])
    num.edge.combos <- convert_table2array(tmp, edge.type = "edge")
    list(nonedge =  num.nonedge.combos, inner.edge.v = num.inneredge.combos.y, inner.edge.h = num.inneredge.combos.x, edge = num.edge.combos)
}

cntrbtns2lple <- function(params, x, n.like.nhbrs.horiz, n.like.nhbrs.vert, max.nhbrs.horiz, max.nhbrs.vert) {
    params[1] -> alpha
    params[2] -> beta.horiz
    params[3] -> beta.vert
    -log(1 + exp(alpha * (1 - 2 * x) + beta.horiz * (max.nhbrs.horiz - 2 *  n.like.nhbrs.horiz) + beta.vert * (max.nhbrs.vert - 2 * n.like.nhbrs.vert)))
}

totalcntrbtns2lple <- function(params, cts.arr, edge.type = "nonedge") {
    max.nhbrs.horiz <- ifelse((edge.type == "inner.edge.h") | (edge.type == "edge"), 1, 2)
    max.nhbrs.vert <- ifelse((edge.type == "inner.edge.v") | (edge.type == "edge"), 1, 2)
    sum <- 0 -> sumstat
    for (i in 0:1) {
        for (j in 0:max.nhbrs.horiz) {
            for (k in 0:max.nhbrs.vert) {
                num.combos  <- cts.arr[i+1, j+1, k+1]
                sum <- sum + num.combos * cntrbtns2lple(params, x = i, n.like.nhbrs.horiz = k, n.like.nhbrs.vert = j, max.nhbrs.horiz, max.nhbrs.vert)
                sumstat  <- sumstat  + num.combos
            }
        }
    }
    sum
}

negmple <- function(params, count.all.1stordr.nhbrs, homogenous = TRUE) {
    if (homogenous)
        params <- params[c(1:2,2)]
    neg.llsum  <- 0
    for (edge_type in c("nonedge", "inner.edge.h", "inner.edge.v", "edge")) {
        neg.llsum <- neg.llsum - totalcntrbtns2lple(params, cts.arr = count.all.1stordr.nhbrs[[eval(quote(expr = edge_type))]], edge.type = edge_type)
    }
    neg.llsum
}

get.MPLE <- function(x, homogenous = TRUE) {
    count.all.1stordr.nhbrs <- count_all_1stordr_nhbrs(x) 
    if (homogenous) pars <- rep(0, 2) else
                                           pars <- rep(0, 3)
    out <- optim(par = pars, fn = negmple, count.all.1stordr.nhbrs = count.all.1stordr.nhbrs, homogenous = homogenous)
    list(MPLES = out$par, neg.MPLvalue = out$value)
}
