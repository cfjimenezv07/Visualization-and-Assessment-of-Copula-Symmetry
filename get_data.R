###########################################################################
# Install packages
###########################################################################
library(copula)
###########################################################################
# Set a working directory
###########################################################################
# wd

# testing function for (reflection) symmetry
f.sym <- function(U){
  function (t, v){
    C.n(cbind(t, v), U) - C.n(cbind(v, t), U)
  }
}

get.y <- function(U, n.curves, t){
set.seed(123)
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves != n){
    v <- runif(n.curves)
  } else {
    v <- U[,2]
    n.curves <- n
  }
  y <- matrix(f.sym(U)(rep(t, n.curves), rep(v, each=p)), nrow = p)
  ind <- sample(1:n.curves, floor(n.curves/2), replace=FALSE)
  y[,ind] <- (-1)*y[,ind]

  return(y)
}

# testing function for radial symmetry 
f.rsym <- function(U){
  function (t, v){
    C.n(cbind(t, v), U) - C.n(cbind(1-t, 1-v), U) + 1 - t - v
  }
}

get.y.r <- function(U, n.curves, t){
  set.seed(123)
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves != n){
    v <- runif(n.curves)
  } else {
    v <- U[,2]
    n.curves <- n
  }
  y <- matrix(f.rsym(U)(rep(t, n.curves), rep(v, each=p)), nrow = p)
  ind <- sample(1:n.curves, floor(n.curves/2), replace=FALSE)
  y[,ind] <- (-1)*y[,ind]
  return(y)
}

# testing function for joint symmetry 
f.jsym.1 <- function(U){
  function (t, v){
    C.n(cbind(t, v), U) + C.n(cbind(t, 1-v), U) - t
  }
}

f.jsym.2 <- function(U){
  function (t, v){
    C.n(cbind(t, v), U) + C.n(cbind(1-t, v), U) - v
  }
}

get.y.j1 <- function(U, n.curves, t){
  set.seed(123)
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves != n){
    v <- runif(n.curves)
  } else {
    v <- U[,2]
    n.curves <- n
  }
  y <- matrix(f.jsym.1(U)(rep(t, n.curves), rep(v, each=p)), nrow = p)
  ind <- sample(1:n.curves, floor(n.curves/2), replace=FALSE)
  y[,ind] <- (-1)*y[,ind]
  return(y)
}

get.y.j2 <- function(U, n.curves, t){
  set.seed(123)
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves != n){
    v <- runif(n.curves)
  } else {
    v <- U[,2]
    n.curves <- n
  }
  y <- matrix(f.jsym.2(U)(rep(t, n.curves), rep(v, each=p)), nrow = p)
  ind <- sample(1:n.curves, floor(n.curves/2), replace=FALSE)
  y[,ind] <- (-1)*y[,ind]
  return(y)
}

