
###########################################################################
# Set a working directory
###########################################################################
setwd("~/Desktop/PhD project 3/Rcodes")
# define a folder within wd where the sizes and powers will be saved.
dirl.r <- "~/Desktop/PhD project 3/Simulation_results/Joint_Symmetry/"

###########################################################################
# Install packages
##########################################################################
packages <- c("devtools", "deSolve", "locfit","foreach","doParallel", 
              "copula", "VineCopula","fdaoutlier","rainbow","fda")

## Now load or install&load all
package_check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


###########################################################################
# Load Auxiliary files
###########################################################################
source("Aux_sym_test.R")

# N stands for the number of simulations to estimate the empirical size / power 
size <- function(cpl, n, sig=0.05, p=100, N=1000, N.b=250, n.curves=NULL,name=c("Indep","Frank","Gausian","Clayton","Gumbel"),tau=NULL){
  
  n.cores <- parallel::detectCores() - 5
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  # setup for time parameter
  e <- 10^(-2)
  t <- seq(e, 1 - e, len=p)
  n.s <- n
  if (is.null(n.curves)){
    n.curves <- ceiling(1.2*n)
  }
  
  # hypothesis testing procedure
  res = foreach(a=1:N, .combine="rbind", .inorder=FALSE,
                .packages=c("copula", "VineCopula", "deSolve", "locfit", "fdaoutlier", "foreach", "fda"),
                .export=c("f.jsym.1", "f.jsym.2", "simulate.jsym", "simulate.jsym_1", "get.y.1", "get.y.2", "get.W", "non.out")) %dopar% {
                  set.seed(123+a)
                  U <- rCopula(n, cpl)
                  y.1 <- get.y.1(U, n.curves, t)
                  y.2 <- get.y.2(U, n.curves, t)
                  
                  # use U to construct simulated sample
                  ec <- empCopula(U)
                  U.s <- simulate.jsym(n.s, ec)
                  y.s.1 <- get.y.1(U.s, n.curves, t)
                  y.s.2 <- get.y.2(U.s, n.curves, t)
                  
                  W.1 <- get.W(y.1, y.s.1)
                  W.2 <- get.W(y.2, y.s.2)
                  
                  # bootstrap samples of the test statistic W
                  Ws.b = foreach(b=1:(N.b), .combine="rbind", .inorder=FALSE,
                                 .packages=c("fda","copula", "VineCopula", "deSolve", "locfit", "fdaoutlier"),
                                 .export=c("f.jsym.1", "f.jsym.2", "simulate.jsym", "get.y.1", "get.y.2", "get.W", "non.out")) %do% {
                                   set.seed(123+b)
                                   U.b <- simulate.jsym(n, ec)
                                   y.b.1 <- get.y.1(U.b, n.curves, t)
                                   y.b.2 <- get.y.2(U.b, n.curves, t)
                                   
                                   # use U to construct simulated sample
                                   ec.b <- empCopula(U.b)
                                   U.bs <- simulate.jsym(n.s, ec.b)
                                   y.bs.1 <- get.y.1(U.bs, n.curves, t)
                                   y.bs.2 <- get.y.2(U.bs, n.curves, t)
                                   
                                   W.b.1 <- get.W(y.b.1, y.bs.1)
                                   W.b.2 <- get.W(y.b.2, y.bs.2)
                                   return (c(W.b.1, W.b.2))	
                                 }
                  
                  # calculate p value
                  p.val.1 <- sum(Ws.b[,1] <= W.1)/N.b
                  p.val.2 <- sum(Ws.b[,2] <= W.2)/N.b
                  p.vals <- p.adjust(c(p.val.1, p.val.2), method="BY")
                  return (p.vals)
                }
  
  stopCluster(cl)
  #print(res)
  
  res.min <- apply(res, 1, min)
  size <- sum(res.min <= sig)/N
  # print(res)
  saveRDS(res,paste0(dirl.r,"res_",n,"_",tau,"_",name,".rds"))
  saveRDS(size,paste0(dirl.r,"size_",n,"_",tau,"_",name,".rds"))
  return (size)
}

### experiments

##### experiment setup
#n <- c(100, 250, 500, 1000)
n <- c(1000)
# paramters - test of radial/joint symmetry
tau <- c(1/4, 1/2, 3/4)


print("#### Joint Symmetry ####")
# independnet
print("independent copula")
ic <- indepCopula()
for (k in 1:length(n)){
print(sprintf("n = %f", n[k]))
system.time(print(size(ic, n[k],name="Indep",tau=0)))
}

# frank
print("frank copula")
for (i in 1:length(tau)){
 for (j in 1:length(n)){
   print(sprintf("tau = %f, n = %f", tau[i], n[j]))
   theta <- BiCopTau2Par(5, tau[i])
   fc <- frankCopula(theta)
   system.time(print(size(fc, n[j],name="Frank",tau=tau[i])))
 }
}

## gaussian
print("gaussian copula")
for (i in 1:length(tau)){
  for (j in 1:length(n)){
    print(sprintf("tau = %f, n = %f", tau[i], n[j]))
    theta <- BiCopTau2Par(1, tau[i])
    nc <- normalCopula(theta)
    system.time(print(size(nc, n[j],name="Gaussian",tau=tau[i])))
  }
}

# clayton
print("clayton copula")
for (i in 1:length(tau)){
  for (j in 1:length(n)){
    print(sprintf("tau = %f, n = %f", tau[i], n[j]))
    theta <- BiCopTau2Par(3, tau[i])
    cc <- claytonCopula(theta)
    system.time(print(size(cc, n[j],name="Clayton",tau=tau[i])))
  }
}

## gumbel
print("gumbel copula")
for (i in 1:length(tau)){
  for (j in 1:length(n)){
    print(sprintf("tau = %f, n = %f", tau[i], n[j]))
    theta <- BiCopTau2Par(4, tau[i])
    gc <- gumbelCopula(theta)
    system.time(print(size(gc, n[j],name="Gumbel",tau=tau[i])))
  }
}