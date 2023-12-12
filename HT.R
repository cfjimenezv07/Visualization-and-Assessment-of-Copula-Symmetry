###########################################################################
# Install packages
###########################################################################
library(fda)
library(foreach)
library(doParallel)
library(copula)
library(VineCopula)
library(stats)
###########################################################################
# Set a working directory
###########################################################################
# wd.

###########################################################################
# Load Auxiliary files
###########################################################################
source("get_data.R")
source("tfplot.R")
source("copula_simulation.R")


set.seed(123)
#### subroutines
# get the test statistic W
get.W <- function(y, y.s){
  set.seed(123)
  n <- dim(y)[2]
  n.s <- dim(y.s)[2]
  dps <- fbplot(cbind(y, y.s), plot=F)$depth
  W <- sum(rank(dps, ties.method = "random")[1:n])
  return(W)
}

# exclude outliers in functional data 
non.out <- function(data){
  n <- dim(data)[2]
  idx.out <- fbplot(data, plot=F)$outpoint
  idx <- setdiff(1:n, idx.out)
  if (length(idx) < 0.9*n){
    return(data)
  } else {
    return(data[,idx])
  }
}



# test for symmetry

# U <- nt.nor[,tests[,1]]
# y <- get.y(U, 1200, t)
# N.b = 1000
test.sym <- function(U, y=NULL, p=100, sig=0.05, N.b=250, n.curves=NULL){
  set.seed(123)
  # setup for time parameter 
  e <- 10^(-2)
  t <- seq(e, 1 - e, len=p)
  n <- dim(U)[1]
  n.s <- n
  if (is.null(n.curves)){
    n.curves <- 1.2*n
  }
  
  # hypothesis testing procedure
  if (is.null(y)) {
    y <- get.y(U, n.curves, t)
  } else {
    p <- dim(y)[1]
    e <- 10^(-2)
    t <- seq(e, 1 - e, len=p)
  }
  y <- non.out(y)
  
  # use U to construct simulated sample 
 
  ec <- empCopula(U)
  U.s <- simulate.sym_1(n.s, ec)
  
  y.s <- get.y(U.s, n.curves, t)
  y.s <- non.out(y.s)

  
  
  W <- get.W(y, y.s)
  
  # bootstrap samples of the test statistic W
  n.cores <- parallel::detectCores() - 2
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  # clusterSetRNGStream(cl, 123) 

  Ws.b = foreach(i=1:(N.b), .combine="c",
                 .packages=c("fda", "copula", "VineCopula", "deSolve", "locfit", "fdaoutlier"),
                 .export=c("f.sym", "simulate.sym", "get.y", "get.W", "non.out")) %dopar% {
                   set.seed(123+i)
                   U.b <- simulate.sym(n, ec)
                   y.b <- get.y(U.b, n.curves, t)
                   y.b <- non.out(y.b)
                   
                   # use U to construct simulated sample 
                   ec.b <- empCopula(U.b)
                   U.bs <- simulate.sym(n.s, ec.b)
                   y.bs <- get.y(U.bs, n.curves, t)
                   y.bs <- non.out(y.bs)
                   
                   W.b <- get.W(y.b, y.bs)
                   
                   return (W.b)
                 }
  
  stopCluster(cl)
  
  # calculate p value
  p.val <- sum(Ws.b <= W)/N.b
  return (p.val)
}








# test for radial symmetry
test.rsym <- function(U, y=NULL, p=100, sig=0.05, N.b=250, n.curves=NULL){
  set.seed(123)
  # setup for time parameter 
  n <- dim(U)[1]
  n.s <- n
  if (is.null(n.curves)){
    n.curves <- ceiling(1.2*n)
  }
  
  # hypothesis testing procedure
  if (is.null(y)){
    e <- 10^(-2)
    t <- seq(e, 1 - e, len=p)
    y <- get.y.r(U, n.curves, t)
  } else {
    p <- dim(y)[1]
    e <- 10^(-2)
    t <- seq(e, 1 - e, len=p)
  }
  y <- non.out(y)
  
  # use U to construct simulated sample 
  ec <- empCopula(U)
  U.s <- simulate.rsym_1(n.s, ec)
  y.s <- get.y.r(U.s, n.curves, t)
  y.s <- non.out(y.s)
  
  W <- get.W(y, y.s)
  
  # bootstrap samples of the test statistic W
  n.cores <- parallel::detectCores() - 2
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  Ws.b = foreach(i=1:(N.b), .combine="c",
                 .packages=c("fda", "copula", "VineCopula", "deSolve", "locfit", "fdaoutlier"),
                 .export=c("f.rsym", "simulate.rsym", "get.y.r", "get.W", "non.out")) %dopar% {
                   set.seed(123+i)
                   U.b <- simulate.rsym(n, ec)
                   y.b <- get.y.r(U.b, n.curves, t)
                   y.b <- non.out(y.b)
                   
                   # use U to construct simulated sample 
                   ec.b <- empCopula(U.b)
                   U.bs <- simulate.rsym(n.s, ec.b)
                   y.bs <- get.y.r(U.bs, n.curves, t)
                   y.bs <- non.out(y.bs)
                   
                   W.b <- get.W(y.b, y.bs)
                   
                   return (W.b)
                 }
  
  stopCluster(cl)
  
  # calculate p value
  p.val <- sum(Ws.b <= W)/N.b
  return (p.val)
}

# test for joint symmetry
test.jsym.1 <- function(U, y=NULL, p=100, sig=0.05, N.b=250, n.curves=NULL){
  set.seed(123)
  # setup for time parameter 
  n <- dim(U)[1]
  n.s <- n
  if (is.null(n.curves)){
    n.curves <- 1.2*n
  }
  
  # hypothesis testing procedure
  if (is.null(y)){
    e <- 10^(-2)
    t <- seq(e, 1 - e, len=p)
    y <- get.y.j1(U, n.curves, t)
  } else {
    p <- dim(y)[1]
    e <- 10^(-2)
    t <- seq(e, 1 - e, len=p)
  }
  y <- non.out(y)
  
  # use U to construct simulated sample 
  ec <- empCopula(U)
  U.s <- simulate.jsym_1(n.s, ec)
  y.s <- get.y.j1(U.s, n.curves, t)
  y.s <- non.out(y.s)
  
  W <- get.W(y, y.s)
  
  # bootstrap samples of the test statistic W
  n.cores <- parallel::detectCores() - 2
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  Ws.b = foreach(i=1:(N.b), .combine="c",
                 .packages=c("fda", "copula", "VineCopula", "deSolve", "locfit", "fdaoutlier"),
                 .export=c("f.jsym.1", "simulate.jsym", "get.y.j1", "get.W", "non.out")) %dopar% {
                   set.seed(123+i)
                   U.b <- simulate.jsym(n, ec)
                   y.b <- get.y.j1(U.b, n.curves, t)
                   y.b <- non.out(y.b)
                   
                   # use U to construct simulated sample 
                   ec.b <- empCopula(U.b)
                   U.bs <- simulate.jsym(n.s, ec.b)
                   y.bs <- get.y.j1(U.bs, n.curves, t)
                   y.bs <- non.out(y.bs)
                   
                   W.b <- get.W(y.b, y.bs)
                   
                   return (W.b)
                 }
  
  stopCluster(cl)
  
  # calculate p value
  p.val <- sum(Ws.b <= W)/N.b
  return (p.val)
}

test.jsym.2 <- function(U, y=NULL, p=100, sig=0.05, N.b=250, n.curves=NULL){
  set.seed(123)
  # setup for time parameter 
  n <- dim(U)[1]
  n.s <- n
  if (is.null(n.curves)){
    n.curves <- 1.2*n
  }
  
  # hypothesis testing procedure
  if (is.null(y)){
    e <- 10^(-2)
    t <- seq(e, 1 - e, len=p)
    y <- get.y.j2(U, n.curves, t)
  } else {
    p <- dim(y)[1]
    e <- 10^(-2)
    t <- seq(e, 1 - e, len=p)
  }
  y <- non.out(y)
  
  # use U to construct simulated sample 
  ec <- empCopula(U)
  U.s <- simulate.jsym_1(n.s, ec)
  y.s <- get.y.j2(U.s, n.curves, t)
  y.s <- non.out(y.s)
  
  W <- get.W(y, y.s)
  
  # bootstrap samples of the test statistic W
  n.cores <- parallel::detectCores() - 2
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  Ws.b = foreach(i=1:(N.b), .combine="c",
                 .packages=c("fda", "copula", "VineCopula", "deSolve", "locfit", "fdaoutlier"),
                 .export=c("f.jsym.2", "simulate.jsym", "get.y.j2", "get.W", "non.out")) %dopar% {
                   set.seed(123+i)
                   U.b <- simulate.jsym(n, ec)
                   y.b <- get.y.j2(U.b, n.curves, t)
                   y.b <- non.out(y.b)
                   
                   # use U to construct simulated sample 
                   ec.b <- empCopula(U.b)
                   U.bs <- simulate.jsym(n.s, ec.b)
                   y.bs <- get.y.j2(U.bs, n.curves, t)
                   y.bs <- non.out(y.bs)
                   
                   W.b <- get.W(y.b, y.bs)
                   
                   return (W.b)
                 }
  
  stopCluster(cl)
  
  # calculate p value
  p.val <- sum(Ws.b <= W)/N.b
  return (p.val)
}

