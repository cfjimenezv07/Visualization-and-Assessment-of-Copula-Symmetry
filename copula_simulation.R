###########################################################################
# Set a working directory
###########################################################################
# wd
setwd("~/Desktop/PhD project 3/Rcodes/")

###########################################################################
# Install packages
###########################################################################
library(copula)



###########################################################################
# simulate a reflection symmetric copula 
simulate.sym <- function(n, cpl){

  smpl <- rCopula(n, cpl)
  dcs <- runif(n)
  ind <- (dcs > 0.5)
  smpl[ind,] <- smpl[ind, c(2,1)]
  
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

# simulate a radially symmetric copula 
simulate.rsym <- function(n, cpl){
  smpl <- rCopula(n, cpl)
  dcs <- runif(n)
  ind <- (dcs > 0.5)
  smpl[ind,] <- 1 - smpl[ind, c(1,2)]
  
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


# simulate a jointly symmetric copula 
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