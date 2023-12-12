###########################################################################
# Install packages
###########################################################################
library(copula)
library(VineCopula)
library(VC2copula)
###########################################################################
# Set a working directory
###########################################################################
# wd.

# define a folder where plots must be saved
dirl.p <- "~/My Drive/Spring 2023/STAT 397/PhD project 3/Plots_paper/"

###########################################################################
# Load Auxiliary files
###########################################################################
source("get_data.R")
source("tfplot.R")


###########################################################################
tau <- 0.5 
p <- 100
t <- seq(0.01, 0.99, len=100)

pdf(file = paste0(dirl.p,"Figure1_1.pdf"),
    width = 10, 
    height = 8)

layout(matrix(c(1:3, 3:19), ncol=4, byrow=TRUE), heights = c(2, 9, 10, 10, 12.8))

# plot title 
par(mar = c(0.2, 5, 0.2, 0.5))
plot.new()
text(0.5, 0.5, expression(bold("Symmetry (S)")), cex = 1.6, font=2)

plot.new()
text(0.5, 0.5, expression(bold("Radial Symmetry (R)")), cex = 1.6, font=1)

plot.new()
text(0.5, 0.5, expression(bold("Joint Symmetry (J)")), cex = 1.6, font=1)

# Clayton copula 
theta <- BiCopTau2Par(3, tau)
cc <- claytonCopula(theta)
U <- rCopula(1000, cc)
par(mar = c(0.5, 4.3, 0.2, 0.5))
y <- get.y(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="", xaxt="n",
       y.label="Clayton (S)", yaxt="n", p.value=1, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

par(mar = c(0.5, 4.15, 0.2, 0.5))
y <- get.y.r(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))
y <- get.y.j1(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))
y <- get.y.j2(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt ="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))

# Frank copula 
theta <- BiCopTau2Par(5, tau)
fc <- frankCopula(theta)
U <- rCopula(1000, fc)
par(mar = c(0.5, 4.3, 1.5, 0.5))
y <- get.y(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="", xaxt="n",
       y.label="Frank (S,R)", yaxt="n", p.value=1, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))
par(mar = c(0.5, 4.15, 1.5, 0.5))
y <- get.y.r(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=1, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))
y <- get.y.j1(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))
y <- get.y.j2(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))

# Gaussian copula 
theta <- BiCopTau2Par(1, tau)
nc <- normalCopula(theta)
U <- rCopula(1000, nc)
par(mar = c(0.5, 4.3, 1.5, 0.5))
y <- get.y(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="", xaxt="n",
       y.label="Gaussian (S,R)", yaxt="n", p.value=1, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))
par(mar = c(0.5, 4.15, 1.5, 0.5))
y <- get.y.r(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=1, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))
y <- get.y.j1(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))
y <- get.y.j2(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))

# independent copula 
ic <- indepCopula()
U <- rCopula(1000, ic)
par(mar = c(4.2, 4.3, 1.5, 0.5))
y <- get.y(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="v",
       y.label="Independent (S,R,J)", yaxt="n", p.value=1, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))
par(mar = c(4.2, 4.15, 1.5, 0.5))
y <- get.y.r(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="v",
       y.label="", yaxt="n", p.value=1, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))
y <- get.y.j1(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="v",
       y.label="", yaxt="n", p.value=1, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))
y <- get.y.j2(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="v", 
       y.label="", yaxt="n", p.value=1, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

dev.off()

pdf(file = "~/My Drive/Spring 2023/STAT 397/PhD project 3/Plots_paper/Figure1_2.pdf",
    width = 10, 
    height = 8)

layout(matrix(c(1:3, 3:19), ncol=4, byrow=TRUE), heights = c(2, 9, 10, 10, 12.8))

# plot title 
par(mar = c(0.2, 5, 0.2, 0.5))
plot.new()
text(0.5, 0.5, expression(bold("Symmetry (S)")), cex = 1.6, font=2)

plot.new()
text(0.5, 0.5, expression(bold("Radial Symmetry (R)")), cex = 1.6, font=1)

plot.new()
text(0.5, 0.5, expression(bold("Joint Symmetry (J)")), cex = 1.6, font=1)

# Gumbel copula
theta <- BiCopTau2Par(4, tau)
gc <- gumbelCopula(theta)
U <- rCopula(1000, gc)
par(mar = c(0.5, 4.3, 0.2, 0.5))
y <- get.y(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="", xaxt="n",
       y.label="Gumbel (S)", yaxt="n", p.value=1, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))
par(mar = c(0.5, 4.15, 0.2, 0.5))
y <- get.y.r(U, 500, t)
tfplot(y, x=t, ylim=c(-0.06, 0.06), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))
y <- get.y.j1(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))
y <- get.y.j2(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))

# Marshall-Olkin copula
moc <- moCopula(c(0.55, 0.85))
U <- rCopula(1000, moc)
y <- get.y(U, 500, t)
par(mar = c(0.5, 4.3, 1.5, 0.5))
tfplot(y, x=t, ylim=c(-0.08, 0.08), x.label="", xaxt="n",
       y.label="Marshal-Olkin (AS)", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.08, -0.03, 0.03, 0.08))
par(mar = c(0.5, 4.15, 1.5, 0.5))
y <- get.y.r(U, 500, t)
tfplot(y, x=t, ylim=c(-0.08, 0.08), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.08, -0.03, 0.03, 0.08))
y <- get.y.j1(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))
y <- get.y.j2(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))

# Tawn Type 1 
tct1 <- tawnT1Copula(c(4.28, 0.6))
U <- rCopula(1000, tct1)
y <- get.y(U, 500, t)
par(mar = c(0.5, 4.3, 1.5, 0.5))
tfplot(y, x=t, ylim=c(-0.08, 0.08), x.label="", xaxt="n",
       y.label="Tawn T1 (AS)", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.08, -0.03, 0.03, 0.08))
par(mar = c(0.5, 4.15, 1.5, 0.5))
y <- get.y.r(U, 500, t)
tfplot(y, x=t, ylim=c(-0.08, 0.08), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.08, -0.03, 0.03, 0.08))
y <- get.y.j1(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))
y <- get.y.j2(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="", xaxt="n",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))

# Tawn Type 2 
tct2 <- tawnT2Copula(c(4.28, 0.6))
U <- rCopula(1000, tct2)
y <- get.y(U, 500, t)
par(mar = c(4.2, 4.3, 1.5, 0.5))
tfplot(y, x=t, ylim=c(-0.08, 0.08), x.label="v", 
       y.label="Tawn T2 (AS)", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.08, -0.03, 0.03, 0.08))
par(mar = c(4.2, 4.15, 1.5, 0.5))
y <- get.y.r(U, 500, t)
tfplot(y, x=t, ylim=c(-0.08, 0.08), x.label="v",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.08, -0.03, 0.03, 0.08))
y <- get.y.j1(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="v",
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))
y <- get.y.j2(U, 500, t)
tfplot(y, x=t, ylim=c(-0.3, 0.3), x.label="v", 
       y.label="", yaxt="n", p.value=0.001, p.value.show=F, cex.lab=1.8)
axis(side = 2, at=c(-0.3, -0.1, 0.1, 0.3))

dev.off()

