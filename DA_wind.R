#Visualization Wind data_1
###########################################################################
# Set a working directory
###########################################################################
# wd
setwd("~/Desktop/PhD project 3/Rcodes/")

# define a folder where plots must be saved
dirl.p <- "~/Desktop/PhD project 3/Plots_paper/DA_wind_1/"
###########################################################################
# Install packages
###########################################################################
library(copula)
library(VineCopula)
library(VC2copula)

###########################################################################
# Load Auxiliary files
###########################################################################
source("get_data.R")
source("tfplot.R")
source("HT.R")


###########################################################################

tau <- 0.5 
p <- 100
t <- seq(0.01, 0.99, len=100)
n.curves <- 1200
N.b <- 1000
set.seed(123)

# explore and clean the current data 
wd <- readRDS(file="Wind_Dataset1.rds")
head(wd)
dim(wd$ws)
# the dataset has no header 
ws <- t(wd$ws)
ds <- dim(ws)

# remove the effect of marginals   
n <- ds[1]
ws.nor <- pobs(as.matrix(ws))
head(ws.nor)

# run data analysis 
tests <- combn(1:3, 2)

# # scatter plots 
par(mfrow=c(2, 3))
for (i in 1:(dim(tests)[2])){
  pair <- tests[,i]
  v1 <- sprintf("Position %.0f", pair[1])
  v2 <- sprintf("Position %.0f", pair[2])
  title <-  paste(v1, " v.s ", v2)
  plot(ws.nor[, pair[1]], ws.nor[, pair[2]], pch = 16, main=title, xlab = v1, ylab = v2)
}

##### plots for the wind_data_set_1

pdf(file = paste0(dirl.p,"DA_wind_1_b1000.pdf"),
    width = 12, 
    height = 13)

layout(matrix(c(1:3, 3:23), ncol=4, byrow=TRUE), heights = c(2, 9, 10, 10, 10, 12))

# plot title 
par(mar = c(0.2, 5, 0.2, 0.5))
plot.new()
text(0.5, 0.5, expression(bold("Symmetry (S)")), cex = 1.6, font=2)

plot.new()
text(0.5, 0.5, expression(bold("Radial Symmetry (R)")), cex = 1.6, font=1)

plot.new()
text(0.5, 0.5, expression(bold("Joint Symmetry (J)")), cex = 1.6, font=1)

#####################
# Position 1 vs Position 2
par(mar = c(0.5, 5, 0.2, 0.5))
U <- ws.nor[,tests[,1]]


f.y <- paste0(dirl.p,"data_matrices/P1_P2_sym.txt")
f.p <- paste0(dirl.p,"p_values/P1_P2_sym.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U,  n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U=U, y=y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", yaxt="n", ylim=c(-0.06, 0.06),
       y.label="Position 1 / Position 2", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

par(mar = c(0.5, 4.3, 0.2, 1.2))

f.y <- paste0(dirl.p,"data_matrices/P1_P2_rsym.txt")
f.p <- paste0(dirl.p,"p_values/P1_P2_rsym.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", yaxt="n", ylim=c(-0.06, 0.06), p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))


f.y <- paste0(dirl.p,"data_matrices/P1_P2_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/P1_P2_jsym_1.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

f.y <- paste0(dirl.p,"data_matrices/P1_P2_jsym_2.txt")
f.p <- paste0(dirl.p,"p_values/P1_P2_jsym_2.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)


####################
# Position 1 vs Position 3 
par(mar = c(0.5, 5, 1.5, 0.5))
U <- ws.nor[,tests[,2]]

f.y <- paste0(dirl.p,"data_matrices/P1_P3_sym.txt")
f.p <- paste0(dirl.p,"p_values/P1_P3_sym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", yaxt="n", ylim=c(-0.06, 0.06),
       y.label="Position 1 / Position 3", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

par(mar = c(0.5, 4.3, 1.5, 1.2))

f.y <- paste0(dirl.p,"data_matrices/P1_P3_rsym.txt")
f.p <- paste0(dirl.p,"p_values/P1_P3_rsym.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", yaxt="n", ylim=c(-0.06, 0.06), p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))


f.y <- paste0(dirl.p,"data_matrices/P1_P3_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/P1_P3_jsym_1.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)


f.y <- paste0(dirl.p,"data_matrices/P1_P3_jsym_2.txt")
f.p <- paste0(dirl.p,"p_values/P1_P3_jsym_2.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

######################
# Position 2 vs Position 3
par(mar = c(0.5, 5, 1.5, 0.5))
U <- ws.nor[,tests[,3]]

f.y <- paste0(dirl.p,"data_matrices/P2_P3_sym.txt")
f.p <- paste0(dirl.p,"p_values/P2_P3_sym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", yaxt="n", ylim=c(-0.06, 0.06),
       y.label="Position 2 / Position 3", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

par(mar = c(0.5, 4.3, 1.5, 1.2))
f.y <- paste0(dirl.p,"data_matrices/P2_P3_rsym.txt")
f.p <- paste0(dirl.p,"p_values/P2_P3_rsym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", yaxt="n", ylim=c(-0.06, 0.06), p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

f.y <- paste0(dirl.p,"data_matrices/P2_P3_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/P2_P3_jsym_1.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
  
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)
f.y <- paste0(dirl.p,"data_matrices/P2_P3_jsym_2.txt")

f.p <- paste0(dirl.p,"p_values/P2_P3_jsym_2.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  
  
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = 1000, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

dev.off()

