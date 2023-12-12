###########################################################################
# Install packages
###########################################################################
library(fda)
library(copula)
library(MASS)
library(mnormt)
library(VC2copula)


###########################################################################
# Set a working directory
###########################################################################


##  generate n random curves with various bivariate copulas 

d <- 2
n <- 500
p <- 100

set.seed(5)
temp <- runif(n) # u2 fixed
U2 <- rep(temp, each = p)
e <- 10^(-2)
t <- seq(0, 1, len=p)
rt <- rep(t, n)

m <- 1000

## general subroutines 
t.plot = function(y.emp, y.theo, l=-0.1, u=0.1){
  par(mfrow=c(1,2))
  fbplot(y.emp, x=t, method="MBD", xlim=c(0, 1),
         main="Empirical")
  abline(h = 0, col="darkblue", lty="dotted")
  fbplot(y.theo, x=t, method="MBD", xlim=c(0,1), 
         main="Theoretical")
  abline(h = 0, col="darkblue", lty="dotted")
  par(mfrow=c(1,1))
}

t.plotMstable = function(y.emp, y.theo, l=-0.1, u=0.1){
  par(mfrow=c(1,2))
  fbplot(y.emp, x=t, main="Empirical")
  abline(h = 0, col="darkblue", lty="dotted")
  fbplot(y.theo, x=t, main="Theoretical")
  abline(h = 0, col="darkblue", lty="dotted")
  par(mfrow=c(1,1))
}

###### test for symmetry ######

## auxiliary functions 
f.sym = function(cpl, u2=U2){
  function (t){
    pCopula(cbind(t, u2), cpl) - pCopula(cbind(u2, t), cpl)
  }
} 

t.sym = function(cpl, u2=U2){
  U<- rCopula(m, cpl)
  y.emp <- matrix( C.n(cbind(rt, u2), U) - C.n(cbind(u2, rt), U), p, n)
  y.theo <- matrix(f.sym(cpl, u2)(rt), p, n)
  t.plot(y.emp, y.theo)
}

c.contour = function(cpl){
  contourplot2(cpl, FUN = pCopula, xlab="u", ylab="v")
}

# Ali-Makhail-Haq 
set.seed(1)
amhc <- amhCopula(-0.8)
c.contour(amhc)
t.sym(amhc)

# Clayton
set.seed(2)
cc <- claytonCopula(4)
c.contour(cc)
t.sym(cc)

# independent 
set.seed(3)
ic <- indepCopula()
c.contour(ic)
t.sym(ic)

# Gumbel-Hougaard
set.seed(7)
ghc <- gumbelCopula(3)
c.contour(ghc)
t.sym(ghc)

# frank 
set.seed(5)
fc <- frankCopula(-9)
c.contour(fc)
t.sym(fc)

# Joe 
set.seed(51)
jc <- joeCopula(2)
c.contour(jc)
t.sym(jc)

# gaussian
set.seed(6)
# gaussianc <- normalCopula()
# c.contour(gaussianc)
# t.sym(gaussianc)
mu <- rep(0, d)
sigma <- diag(2)
X <- mvrnorm(m, mu, sigma)
U <- pnorm(X, 0, 1)
y.emp <- matrix( C.n(cbind(rt, U2), U) - C.n(cbind(U2, rt), U), p, n)
y.theo <- matrix(pmnorm(cbind(qnorm(rt), qnorm(U2)), mu, sigma)-pmnorm(cbind(qnorm(U2), qnorm(rt)), mu, sigma), p, n)
t.plot(y.emp, y.theo)

# student t 
set.seed(7)
tc <- tCopula(0.2)
c.contour(tc)
t.sym(tc)

# frechet Hoeffding lower bound 
set.seed(8)
U <- runif(m)
U <- cbind(U, 1-U)
y.emp <- matrix( C.n(cbind(rt, U2), U) - C.n(cbind(U2, rt), U), p, n)
y.theo <- matrix(pmax(rt+U2-1, 0) - pmax(U2+rt-1, 0), p, n)
t.plot(y.emp, y.theo)

# frechet Hoeffding upper bound
set.seed(9)
U <- runif(m)
U <- cbind(U, U)
y.emp <- matrix( C.n(cbind(rt, U2), U) - C.n(cbind(U2, rt), U), p, n)
y.theo <- matrix(pmin(rt, U2) - pmin(U2, rt), p, n)
t.plot(y.emp, y.theo)

# farlie-Gumbel-Morgenstern 
set.seed(10)
fgmc <- fgmCopula(-1)
c.contour(fgmc)
t.sym(fgmc)

# Plackett 
set.seed(11)
plkc <- plackettCopula(0.5)
c.contour(plkc)
t.sym(plkc)

# Galambos 
set.seed(12)
glbc <- galambosCopula(2)
c.contour(glbc)
t.sym(glbc)

# Cuadras-Auge
set.seed(18)
cac <- moCopula(c(0.3, 0.3))
c.contour(cac)
t.sym(cac)


# Tawn 
## type 1
set.seed(17)
tct1 <- tawnT1Copula(c(10, 0.9))
c.contour(tct1)
t.sym(tct1)
par(mfrow=c(1,1))
curve(f.sym(tct1, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.sym(tct1, u2=temp[i])(x), ylab="", add=TRUE)
}

## type 2
tct2 <- tawnT2Copula(c(10, 0.9))
c.contour(tct2)
t.sym(tct2)
curve(f.sym(tct2, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.sym(tct2, u2=temp[i])(x), ylab="", add=TRUE)
}

# Marshall-Olkin 
set.seed(20)
moc <- moCopula(c(0.1, 0.8))
c.contour(moc)
t.sym(moc)
curve(f.sym(moc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.sym(moc, u2=temp[i])(x), ylab="", add=TRUE)
}

moc <- moCopula(c(0.4, 0.6))
c.contour(moc)
t.sym(moc)
curve(f.sym(moc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.sym(moc, u2=temp[i])(x), ylab="", add=TRUE)
}

moc <- moCopula(c(0.3, 0.7))
c.contour(moc)
t.sym(moc)
curve(f.sym(moc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.sym(moc, u2=temp[i])(x), ylab="", add=TRUE)
}

# Khoudragi's device 
## independent + Clayton/0.25
icc.25 <- khoudrajiCopula(copula1 = ic, copula2 = cc, shapes = c(1/4, 1))
c.contour(icc.25)
t.sym(icc.25)
curve(f.sym(icc.25, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.sym(icc.25, u2=temp[i])(x), ylab="", add=TRUE)
}
## independent + Clayton/0.5
icc.5 <- khoudrajiCopula(copula1 = ic, copula2 = cc, shapes = c(1/2, 1))
c.contour(icc.5)
t.sym(icc.5)
curve(f.sym(icc.5, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.sym(icc.5, u2=temp[i])(x), ylab="", add=TRUE)
}
## independent + Clayton/0.75
icc.75 <- khoudrajiCopula(copula1 = ic, copula2 = cc, shapes = c(3/4, 1))
c.contour(icc.75)
t.sym(icc.75)
curve(f.sym(icc.75, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.sym(icc.75, u2=temp[i])(x), ylab="", add=TRUE)
}
## independent + Gumbel-H/0.25
ighc.25 <- khoudrajiCopula(copula1 = ic, copula2 = ghc, shapes = c(1/4, 1))
c.contour(ighc.25)
t.sym(ighc.25)
curve(f.sym(ighc.25, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.sym(ighc.25, u2=temp[i])(x), ylab="", add=TRUE)
}
## independent + Gumbel-H/0.5
ighc.5 <- khoudrajiCopula(copula1 = ic, copula2 = ghc, shapes = c(1/2, 1))
c.contour(ighc.5)
t.sym(ighc.5)
curve(f.sym(ighc.5, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.sym(ighc.5, u2=temp[i])(x), ylab="", add=TRUE)
}
## independent + Gumbel-H/0.75
ighc.75 <- khoudrajiCopula(copula1 = ic, copula2 = ghc, shapes = c(3/4, 1))
c.contour(ighc.75)
t.sym(ighc.75)
curve(f.sym(ighc.75, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.sym(ighc.75, u2=temp[i])(x), ylab="", add=TRUE)
}


###### test for radial symmetry ######
## auxiliary functions 
f.rsym = function(cpl, u2=U2){
  function (t){
    (pCopula(cbind(t, u2), cpl) - pCopula(cbind(1 - t, 1- u2), cpl) + 1 - t - u2)^2
  }
} 
t.rsym = function(cpl, u2=U2){
  U <- rCopula(m, cpl)
  y.emp <- matrix( (C.n(cbind(rt, u2), U) - C.n(cbind(1-rt, 1-u2), U) + 1 - rt - u2)^2, p, n)
  y.theo <- matrix(f.rsym(cpl, u2)(rt), p, n)
  t.plot(y.emp, y.theo)
}

# Ali-Makhail-Haq 
set.seed(101)
amhc <- amhCopula(-0.8)
t.rsym(amhc)
curve(f.rsym(amhc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n,30)){
  curve(f.rsym(amhc, u2=temp[i])(x), ylab="", add=TRUE)
}

# Clayton
set.seed(102)
cc <- claytonCopula(4)
t.rsym(cc)
curve(f.rsym(cc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(cc, u2=temp[i])(x), ylab="", add=TRUE)
}

# independent 
set.seed(103)
ic <- indepCopula()
t.rsym(ic)

# Gumbel-Hougaard
set.seed(104)
ghc <- gumbelCopula(3)
t.rsym(ghc)
curve(f.rsym(ghc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(ghc, u2=temp[i])(x), ylab="", add=TRUE)
}

# frank 
set.seed(105)
fc <- frankCopula(-9)
t.rsym(fc)

# Joe 
set.seed(106)
jc <- joeCopula(2)
t.rsym(jc)
curve(f.rsym(jc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(jc, u2=temp[i])(x), ylab="", add=TRUE)
}

# gaussian
set.seed(107)
# gaussianc <- normalCopula()
# c.contour(gaussianc)
# t.sym(gaussianc)
mu <- rep(0, d)
sigma <- diag(2)
X <- mvrnorm(m, mu, sigma)
U <- pnorm(X, 0, 1)
y.emp <- matrix( C.n(cbind(rt, U2), U) - C.n(cbind(1- rt, 1 - U2), U) + 1 - rt - U2, p, n)
y.theo <- matrix(pmnorm(cbind(qnorm(rt), qnorm(U2)), mu, sigma)
                 -pmnorm(cbind(qnorm(1-rt), qnorm(1-U2)), mu, sigma) + 1 - rt - U2, p, n)
t.plot(y.emp, y.theo)

# student t 
set.seed(108)
tc <- tCopula(0.2)
t.rsym(tc)

# frechet Hoeffding lower bound 
set.seed(109)
U <- runif(m)
U <- cbind(U, 1-U)
y.emp <- matrix( C.n(cbind(rt, U2), U) - C.n(cbind(1-rt, 1-U2), U) 
                 + 1 - rt - U2, p, n)
y.theo <- matrix(pmax(rt+U2-1, 0) - pmax(1-rt-U2, 0) + 1 - rt - U2, p, n)
t.plot(y.emp, y.theo)

# frechet Hoeffding upper bound
set.seed(110)
U <- runif(m)
U <- cbind(U, U)
y.emp <- matrix( C.n(cbind(rt, U2), U) - C.n(cbind(1-rt, 1-U2), U)
                 + 1 - rt - U2, p, n)
y.theo <- matrix(pmin(rt, U2) - pmin(1-rt, 1-U2) + 1 - rt - U2, p, n)
t.plot(y.emp, y.theo)

# farlie-Gumbel-Morgenstern 
set.seed(111)
fgmc <- fgmCopula(-1)
t.rsym(fgmc)

# Plackett 
set.seed(112)
plkc <- plackettCopula(0.5)
t.rsym(plkc)

# Galambos 
set.seed(113)
glbc <- galambosCopula(2)
t.rsym(glbc)
curve(f.rsym(glbc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(glbc, u2=temp[i])(x), ylab="", add=TRUE)
}

# Cuadras-Auge
set.seed(114)
cac <- moCopula(c(0.3, 0.3))
t.rsym(cac)
curve(f.rsym(cac, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(cac, u2=temp[i])(x), ylab="", add=TRUE)
}

# Tawn 
## type 1
set.seed(115)
tct1 <- tawnT1Copula(c(10, 0.9))
t.rsym(tct1)
par(mfrow=c(1,1))
curve(f.rsym(tct1, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(tct1, u2=temp[i])(x), ylab="", add=TRUE)
}

## type 2
set.seed(116)
tct2 <- tawnT2Copula(c(10, 0.9))
t.rsym(tct2)
curve(f.rsym(tct2, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(tct2, u2=temp[i])(x), ylab="", add=TRUE)
}

# Marshall-Olkin 
set.seed(117)
moc <- moCopula(c(0.5, 0.8))
t.rsym(moc)
curve(f.rsym(moc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(moc, u2=temp[i])(x), ylab="", add=TRUE)
}

set.seed(118)
moc <- moCopula(c(0.4, 0.6))
t.rsym(moc)
curve(f.rsym(moc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(moc, u2=temp[i])(x), ylab="", add=TRUE)
}

set.seed(119)
moc <- moCopula(c(0.3, 0.7))
t.rsym(moc)
curve(f.rsym(moc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(moc, u2=temp[i])(x), ylab="", add=TRUE)
}

###### test for joint symmetry ######
## auxiliary functions 
f.jsym = function(cpl, u2=U2){
  function (t){
    5*(pCopula(cbind(t, u2), cpl) + pCopula(cbind(t, 1-u2), cpl) - t)^2 
    + 5*(pCopula(cbind(t, u2), cpl) + pCopula(cbind( 1-t, u2), cpl) - u2)^2
  }
} 
t.jsym = function(cpl, u2=U2, lb=-0.05, ub=0.15){
  U <- rCopula(m, cpl)
  y.emp <- matrix(5*(C.n(cbind(rt, u2), U) + C.n(cbind(rt, 1-u2), U) - rt)^2 
                  + 5*(C.n(cbind(rt, u2), U) + C.n(cbind(1-rt, u2), U) - u2)^2, p, n)
  y.theo <- matrix(f.jsym(cpl, u2)(rt), p, n)
  t.plot(y.emp, y.theo, l=lb, u=ub)
}

# Ali-Makhail-Haq 
set.seed(1001)
amhc <- amhCopula(-0.8)
t.jsym(amhc)
curve(f.jsym(amhc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.05, 0.15), ylab="")
for (i in 2:min(n,30)){
  curve(f.jsym(amhc, u2=temp[i])(x), ylab="", add=TRUE)
}

# Clayton
set.seed(1002)
cc <- claytonCopula(4)
t.jsym(cc, ub=1.5)
curve(f.jsym(cc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.05, 1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.jsym(cc, u2=temp[i])(x), ylab="", add=TRUE)
}

# independent 
set.seed(1003)
ic <- indepCopula()
t.jsym(ic)

# Gumbel-Hougaard
set.seed(104)
ghc <- gumbelCopula(3)
t.jsym(ghc, ub=1.5)
curve(f.jsym(ghc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.05, 1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.jsym(ghc, u2=temp[i])(x), ylab="", add=TRUE)
}

# frank 
set.seed(105)
fc <- frankCopula(-9)
t.rsym(fc)

# Joe 
set.seed(106)
jc <- joeCopula(2)
t.rsym(jc)
curve(f.rsym(jc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(jc, u2=temp[i])(x), ylab="", add=TRUE)
}

# gaussian
set.seed(107)
# gaussianc <- normalCopula()
# c.contour(gaussianc)
# t.sym(gaussianc)
mu <- rep(0, d)
sigma <- diag(2)
X <- mvrnorm(m, mu, sigma)
U <- pnorm(X, 0, 1)
y.emp <- matrix( C.n(cbind(rt, U2), U) - C.n(cbind(1- rt, 1 - U2), U) + 1 - rt - U2, p, n)
y.theo <- matrix(pmnorm(cbind(qnorm(rt), qnorm(U2)), mu, sigma)
                 -pmnorm(cbind(qnorm(1-rt), qnorm(1-U2)), mu, sigma) + 1 - rt - U2, p, n)
t.plot(y.emp, y.theo)

# student t 
set.seed(108)
tc <- tCopula(0.2)
t.rsym(tc)

# frechet Hoeffding lower bound 
set.seed(109)
U <- runif(m)
U <- cbind(U, 1-U)
y.emp <- matrix( C.n(cbind(rt, U2), U) - C.n(cbind(1-rt, 1-U2), U) 
                 + 1 - rt - U2, p, n)
y.theo <- matrix(pmax(rt+U2-1, 0) - pmax(1-rt-U2, 0) + 1 - rt - U2, p, n)
t.plot(y.emp, y.theo)

# frechet Hoeffding upper bound
set.seed(110)
U <- runif(m)
U <- cbind(U, U)
y.emp <- matrix( C.n(cbind(rt, U2), U) - C.n(cbind(1-rt, 1-U2), U)
                 + 1 - rt - U2, p, n)
y.theo <- matrix(pmin(rt, U2) - pmin(1-rt, 1-U2) + 1 - rt - U2, p, n)
t.plot(y.emp, y.theo)

# farlie-Gumbel-Morgenstern 
set.seed(111)
fgmc <- fgmCopula(-1)
t.rsym(fgmc)

# Plackett 
set.seed(112)
plkc <- plackettCopula(0.5)
t.rsym(plkc)

# Galambos 
set.seed(113)
glbc <- galambosCopula(2)
t.rsym(glbc)
curve(f.rsym(glbc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(glbc, u2=temp[i])(x), ylab="", add=TRUE)
}

# Cuadras-Auge
set.seed(114)
cac <- moCopula(c(0.3, 0.3))
t.rsym(cac)
curve(f.rsym(cac, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(cac, u2=temp[i])(x), ylab="", add=TRUE)
}

# Tawn 
## type 1
set.seed(115)
tct1 <- tawnT1Copula(c(10, 0.9))
t.rsym(tct1)
par(mfrow=c(1,1))
curve(f.rsym(tct1, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(tct1, u2=temp[i])(x), ylab="", add=TRUE)
}

## type 2
set.seed(116)
tct2 <- tawnT2Copula(c(10, 0.9))
t.rsym(tct2)
curve(f.rsym(tct2, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(tct2, u2=temp[i])(x), ylab="", add=TRUE)
}

# Marshall-Olkin 
set.seed(117)
moc <- moCopula(c(0.1, 0.8))
t.rsym(moc)
curve(f.rsym(moc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(moc, u2=temp[i])(x), ylab="", add=TRUE)
}

set.seed(118)
moc <- moCopula(c(0.4, 0.6))
t.rsym(moc)
curve(f.rsym(moc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(moc, u2=temp[i])(x), ylab="", add=TRUE)
}

set.seed(119)
moc <- moCopula(c(0.3, 0.7))
t.rsym(moc)
curve(f.rsym(moc, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(moc, u2=temp[i])(x), ylab="", add=TRUE)
}

###### test for max-stability ##### 
## auxiliary functions 
f.mstable = function(cpl, U, P=p){
  function(x){
    rep(pCopula(cbind(U[,1], U[,2]), cpl), each=P)-((pCopula(cbind(rep(U[,1], each=P)^(1/x), rep(U[,2], each=P)^(1/x)), cpl))^x)
  }
} 
t.mstable = function(cpl, U, lb=-0.1, ub=0.1, P=p){
  y.emp <- matrix(rep(C.n(cbind(U[,1], U[,2]), U), each=P) 
                  - ((C.n(cbind(rep(U[,1], each=P)^(1/rt), rep(U[,2], each=P)^(1/rt)), U))^rt), P, n)
  y.theo <- matrix(f.mstable(cpl, U)(rt), P, n)
  t.plotMstable(y.emp, y.theo, l=lb, u=ub)
}

# experiment setup 
t <- 2:(p+1)
rt <- rep(t, n)


# Ali-Makhail-Haq 
set.seed(2001)
amhc <- amhCopula(-0.8)
U <- rCopula(n, amhc)
t.mstable(amhc, U)
# curve(f.mstable(amhc, matrix(U[1,], 1, 2), P = 1)(x), ylab="")
# for (i in 2:min(n,30)){
#   curve(f.mstable(amhc, matrix(U[i,], 1, 2), P =1)(x), ylab="", add=TRUE)
# }

# Galambos 
set.seed(213)
glbc <- galambosCopula(2)
U <- rCopula(n, glbc)
t.mstable(glbc, U)

# Cuadras-Auge
set.seed(214)
cac <- moCopula(c(0.3, 0.3))
U <- rCopula(n, cac)
t.mstable(cac, U)

# Tawn 
# type 1
set.seed(115)
tct1 <- tawnT1Copula(c(10, 0.9))
t.rsym(tct1)
par(mfrow=c(1,1))
curve(f.rsym(tct1, u2=temp[1])(x), xlim=c(e, 1-e), ylim=c(-0.1, 0.1), ylab="")
for (i in 2:min(n, 30)){
  curve(f.rsym(tct1, u2=temp[i])(x), ylab="", add=TRUE)
}



