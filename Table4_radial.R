###########################################################################
# Set a working directory
###########################################################################
setwd("~/Desktop/PhD project 3/Rcodes")

# define a folder within wd where the sizes and powers will be saved.
dirl.r <- "~/Desktop/PhD project 3/Simulation_results/Radial_symmetry/"

###########################################################################
# Install packages
###########################################################################
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


# get the empirical size/power of experiments
size <- function(cpl, n, sig=0.05, p=100, N=1000, N.b=250, n.curves=NULL,name=c("Indep","Frank","Gausian","Clayton","Gumbel"),tau=NULL){
  
  n.cores <- parallel::detectCores() - 4
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
  res = foreach(a=1:N, .combine="c",
                .packages=c("fda", "copula", "VineCopula", "deSolve", "locfit", "fdaoutlier", "foreach", "fda"),
                .export=c("f.rsym", "simulate.rsym", "simulate.rsym_1", "get.y", "get.W")) %dopar% {
                  set.seed(123+a)
                  U <- rCopula(n, cpl)
                  y <- get.y.r(U, n.curves, t)
                  
                  # use U to construct simulated sample
                  ec <- empCopula(U)
                  U.s <- simulate.rsym(n.s, ec)
                  y.s <- get.y.r(U.s, n.curves, t)
                  
                  W <- get.W(y, y.s)
                  
                  # bootstrap samples of the test statistic W
                  Ws.b = foreach(b=1:(N.b), .combine="c",
                                 .packages=c("fda", "copula", "VineCopula", "deSolve", "locfit", "fdaoutlier"),
                                 .export=c("f.rsym", "simulate.rsym", "get.y", "get.W")) %do% {
                                   set.seed(123+b)
                                   U.b <- simulate.rsym(n, ec)
                                   y.b <- get.y.r(U.b, n.curves, t)
                                   
                                   # use U to construct simulated sample
                                   ec.b <- empCopula(U.b)
                                   U.bs <- simulate.rsym(n.s, ec.b)
                                   y.bs <- get.y.r(U.bs, n.curves, t)
                                   
                                   W.b <- get.W(y.b, y.bs)
                                   return (W.b)	
                                 }
                  
                  # calculate p value
                  p.val <- sum(Ws.b <= W)/N.b
                  return (p.val)
                }
  
  stopCluster(cl)
  size <- sum(res <= sig)/N
  # print(res)
  saveRDS(res,paste0(dirl.r,"res_",n,"_",tau,"_",name,".rds"))
  saveRDS(size,paste0(dirl.r,"size_",n,"_",tau,"_",name,".rds"))
  return (size)
}

### experiments

##### experiment setup
n <- c(100, 250, 500, 1000)
# paramters - test of radial/joint.rsymmetry
tau <- c(1/4, 1/2, 3/4)


# name=c("Indep","Clayton","Gausian","Gumbel")
#### test simulated.rsymmetric copula
ic <- indepCopula()
for (k in 1:length(n)){
  print(sprintf("n = %f", n[k]))
  system.time(print(size(ic, n[k],name="Indep",tau=0)))
}

print("#### Radial Symmetry ####")
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