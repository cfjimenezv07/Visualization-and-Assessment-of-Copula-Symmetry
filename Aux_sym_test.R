# Aux#### subroutines
###########################################################################
# Set a working directory
###########################################################################
# wd
setwd("~/Desktop/PhD project 3/Rcodes/")


###################################
# Symmetry
###################################
# resampling method
simulate.sym <- function(n, cpl){
  
  smpl <- rCopula(n, cpl)
  dcs <- runif(n)
  more <- (dcs > 0.5) 
  smpl[more,] <- smpl[more, c(2,1)] 
  
  return(smpl)
}


simulate.sym_1 <- function(n, cpl){
  set.seed(123)
  smpl <- rCopula(n, cpl)
  dcs <- runif(n)
  ind <- (dcs > 0.5)
  smpl[ind,] <- smpl[ind, c(2,1)]
  
  return(smpl)
}

# testing function
f.sym <- function(U){
  function (t, u2){
    C.n(cbind(t, u2), U) - C.n(cbind(u2, t), U)
  }
}

# get data matrix for functional boxplot
get.y <- function(U, n.curves, t){
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves != n){
    u1 <- runif(n.curves)
    u2 <- runif(n.curves)
  } else {
    u1 <- U[,1]  
    u2 <- U[,2]
    n.curves <- n
  }
  y <- matrix(f.sym(U)(rep(t, n.curves), rep(u2, each=p)), nrow = p)
  ind <- sample(1:n.curves, floor(n.curves/2), replace=FALSE)
  y[,ind] <- (-1)*y[,ind]
  return(y)
}

# get the test statistic W
get.W <- function(y, y.s){
  n <- dim(y)[2]
  n.s <- dim(y.s)[2]
  
  data <- t(cbind(y, y.s))
  dps.mbd <- modified_band_depth(data) 
  W.mbd <- sum(rank(dps.mbd, ties.method = "random")[1:n])
  
  return(W.mbd)
}

###################################
# Radial
###################################

# resampling method
simulate.rsym <- function(n, cpl){
  
  smpl <- rCopula(n, cpl)
  dcs <- runif(n)
  more <- (dcs > 0.5) 
  smpl[more,] <- 1 - smpl[more, c(1,2)] 
  
  return(smpl)
}

simulate.rsym_1 <- function(n, cpl){
  set.seed(123)
  smpl <- rCopula(n, cpl)
  dcs <- runif(n)
  ind <- (dcs > 0.5)
  smpl[ind,] <- 1 - smpl[ind, c(1,2)]
  
  return(smpl)
}

# testing function
f.rsym = function(U){
  function (t, u2){
    C.n(cbind(t, u2), U) - C.n(cbind(1 - t, 1 - u2), U) + 1 - t - u2
  }
}

# get data matrix for functional boxplot
get.y.r <- function(U, n.curves, t){
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves < n){
    u2 <- sample(U[,2], n.curves, replace = FALSE)
  } else {
    u2 <- U[,2]
    n.curves <- n
  }
  y <- matrix(f.rsym(U)(rep(t, n.curves), rep(u2, each=p)), nrow = p)
  ind <- sample(1:n.curves, floor(n.curves/2), replace=FALSE)
  y[,ind] <- (-1)*y[,ind]
  return(y)
}



###################################
# Joint
###################################

#### subroutines
# resampling method
simulate.jsym <- function(n, cpl){
  
  smpl <- rCopula(n, cpl)
  dcs <- runif(n)
  ind.1 <- (dcs > 0.25 & dcs <= 0.5)
  ind.2 <- (dcs > 0.5 & dcs <= 0.75)
  ind.3 <- (dcs > 0.75) 
  
  smpl[ind.1,] <- cbind(1 - smpl[ind.1, 1], smpl[ind.1, 2])
  smpl[ind.2,] <- cbind(smpl[ind.2, 1], 1 - smpl[ind.2, 2])
  smpl[ind.3,] <- 1 - smpl[ind.3, c(1,2)] 
  
  return(smpl)
}

simulate.jsym_1<- function(n, cpl){
  set.seed(123)
  smpl <- rCopula(n, cpl)
  dcs <- runif(n)
  ind.1 <- (dcs > 0.25 & dcs <= 0.5)
  ind.2 <- (dcs > 0.5 & dcs <= 0.75)
  ind.3 <- (dcs > 0.75)
  
  smpl[ind.1,] <- cbind(1 - smpl[ind.1, 1], smpl[ind.1, 2])
  smpl[ind.2,] <- cbind(smpl[ind.2, 1], 1 - smpl[ind.2, 2])
  smpl[ind.3,] <- 1 - smpl[ind.3, c(1,2)]
  
  return(smpl)
}
# testing function
f.jsym.1 = function(U){
  function (t, u2){
    C.n(cbind(t, u2), U) + C.n(cbind(t, 1 - u2), U) - t
  }
}

f.jsym.2 = function(U){
  function (t, u2){
    C.n(cbind(t, u2), U) + C.n(cbind(1 - t, u2), U) - u2
  }
}

# get data matrix for functional boxplot
get.y.1 <- function(U, n.curves, t){
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves != n){
    u2 <- runif(n.curves)
  } else {
    u2 <- U[,2]
    n.curves <- n
  }
  y <- matrix(f.jsym.1(U)(rep(t, n.curves), rep(u2, each=p)), nrow = p)
  ind <- sample(1:n.curves, floor(n.curves/2), replace=FALSE)
  y[,ind] <- (-1)*y[,ind]
  return(y)
}

get.y.2 <- function(U, n.curves, t){
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves < n){
    u2 <- sample(U[,2], n.curves, replace = FALSE)
  } else {
    u2 <- U[,2]
    n.curves <- n
  }
  y <- matrix(f.jsym.2(U)(rep(t, n.curves), rep(u2, each=p)), nrow = p)
  ind <- sample(1:n.curves, floor(n.curves/2), replace=FALSE)
  y[,ind] <- (-1)*y[,ind]
  return(y)
}

