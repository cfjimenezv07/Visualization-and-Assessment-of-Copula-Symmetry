###########################################################################
# Set a working directory
###########################################################################
setwd("~/Desktop/PhD project 3/Rcodes")

# define a folder within wd where the sizes and powers will be saved.
dirl.r <- "~/Desktop/PhD project 3/Simulation_results/Symmetry/"

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
size <- function(cpl, n, sig=0.05, p=100, N=1000, N.b=250, n.curves=NULL,name=c("Indep","Clayton","Gausian","Gumbel")){
  
  av <- parallel::detectCores()	
  n.cores <- max(floor(av/2) + 6, av - 4)
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  # setup for time parameter
  e <- 10^(-2)
  t <- seq(e, 1 - e, len=p)
  n.s <- n
  
  # number of curves 
  if (is.null(n.curves)) {
    n.curves <- ceiling(1.2*n)
  }
  
  # hypothesis testing procedure
  res = foreach(a=1:N, .combine="c", .inorder=FALSE,
                .packages=c("copula", "VineCopula", "deSolve", "locfit", "fda", "foreach", "fdaoutlier"),
                .export=c("f.sym", "simulate.sym", "simulate.sym_1", "get.y", "get.W") ) %dopar% {
                  set.seed(123+a)
                  U <- rCopula(n, cpl)
                  y <- get.y(U, n.curves, t)
                  
                  # use U to construct simulated sample
                  ec <- empCopula(U)
                  U.s <- simulate.sym(n.s, ec)
                  y.s <- get.y(U.s, n.curves, t)
                  
                  W <- get.W(y, y.s)
                  
                  # bootstrap samples of the test statistic W
                  Ws.b = foreach(b=1:(N.b), .combine="c", .inorder=FALSE,
                                 .packages=c("copula", "VineCopula", "deSolve", "locfit", "fda", "fdaoutlier"),
                                 .export=c("f.sym", "simulate.sym", "get.y", "get.W") ) %do% {
                                   
                                   set.seed(123+b)
                                   U.b <- simulate.sym(n, ec)
                                   y.b <- get.y(U.b, n.curves, t)
                                   
                                   # use U to construct simulated sample
                                   ec.b <- empCopula(U.b)
                                   U.bs <- simulate.sym(n.s, ec.b)
                                   y.bs <- get.y(U.bs, n.curves, t)
                                   
                                   W.b <- get.W(y.b, y.bs)               
                                 }
                  
                  p.val <- sum(Ws.b <= W)/N.b
                  return (p.val)
                }
  
  stopCluster(cl)
  size <- sum(res <= sig)/N
  print(res)
  # saveRDS(res,paste0(dirl.r,"res_",n,"_",name,".rds"))
  # saveRDS(size,paste0(dirl.r,"size_",n,"_",name,".rds"))
  return (size)
}

### experiments

##### experiment setup
n <- c(100, 250, 500, 1000)

# paramaters - test of symmetry
delta <- c(0, 1/4, 1/2, 3/4)
# delta <- c(3/4)
tau.asym <- c(0.5, 0.7, 0.9)
tau.sym <- c(1/4, 1/2, 3/4)
# paramters - test of radial/joint symmetry
tau <- c(1/4, 1/2, 3/4)


# name=c("Indep","Clayton","Gausian","Gumbel")



#### test simulated symmetric copula
print("#### Symmetry ####")
# n<-100
print("independent copula")
ic <- indepCopula()
for (k in 1:length(n)){
print(sprintf("n = %f", n[k]))
system.time(print(size(ic, n[k],name="Indep")))
}


#
### clayton
print("clayton copula")
for (i in 1:length(delta)){
 for (j in 1:3){
   for (k in 1:length(n)){
     if (delta[i] == 0){
       print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.sym[j], n[k]))
       theta <- BiCopTau2Par(3, tau.sym[j])
     
       
       system.time(print(size(cc, n[k],name="Clayton")))
writeLines('\n')
     } else {
       print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.asym[j], n[k]))
       theta <- BiCopTau2Par(3, tau.asym[j])
       cc <- claytonCopula(theta)
       kcc <- khoudrajiCopula(copula1 = ic, copula2 = cc, shapes = c((1 - delta[i]), 1))
       system.time(print(size(kcc, n[k],name="Clayton")))
       writeLines('\n')
     }
   }
 }
}



### gaussian
print("gaussian copula")
for (i in 1:length(delta)){
 for (j in 1:3){
   for (k in 1:length(n)){
     if (delta[i] == 0){
       print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.sym[j], n[k]))
       theta <- BiCopTau2Par(1, tau.sym[j])
       nc <- normalCopula(theta)
       system.time(print(size(nc, n[k],name="Gausian")))
       writeLines('\n')
     } else {
       print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.asym[j], n[k]))
       theta <- BiCopTau2Par(1, tau.asym[j])
       nc <- normalCopula(theta)
       knc <- khoudrajiCopula(copula1 = ic, copula2 = nc, shapes = c((1 - delta[i]), 1))
       system.time(print(size(knc, n[k],name="Gausian")))
       writeLines('\n')
     }
   }
 }
}

## gumbel
print("gumbel copula")
for (i in 1:length(delta)){
  for (j in 1:3){
    for (k in 1:length(n)){
      if (delta[i] == 0){
        print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.sym[j], n[k]))
        theta <- BiCopTau2Par(4, tau.sym[j])
        gc <- gumbelCopula(theta)
        system.time(print(size(gc, n[k],name="Gumbel")))
        writeLines('\n')
      } else {
        print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.asym[j], n[k]))
        theta <- BiCopTau2Par(4, tau.asym[j])
        gc <- gumbelCopula(theta)
        kgc <- khoudrajiCopula(copula1 = ic, copula2 = gc, shapes = c((1 - delta[i]), 1))
        system.time(print(size(kgc, n[k],name="Gumbel")))
        writeLines('\n')
      }
    }
  }
}
