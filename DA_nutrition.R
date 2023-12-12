
###########################################################################
# Set a working directory
###########################################################################
# wd
setwd("~/Desktop/PhD project 3/Rcodes/")

# define a folder where plots must be saved
dirl.p <- "~/Desktop/PhD project 3/Plots_paper/DA_Nutrition/"
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
n.curves <- 2000
N.b <- 2000

set.seed(123)

# explore and clean the current data 
nt <- read.table(file="nutrients.txt", header = TRUE)
head(nt)
drops <- c("ID")
nt <- nt[, !(names(nt) %in% drops)]
head(nt)
ds <- dim(nt)

# remove the effect of marginals   
n <- ds[1]
nt.nor <- as.data.frame(apply(nt, 2, function(c) rank(c, ties.method = "max")/(n+1)))
head(nt.nor)

# run data analysis 
tests <- combn(names(nt.nor), 2)
n.pair <- dim(tests)[2]          # 10 pairs in total

##### scatterplot 
pdf(file = paste0(dirl.p,"DA_Nutrition_scatter.pdf"),
    width = 13.5, 
    height = 6)

par(mfrow=c(2, 5))
for (i in 1:n.pair){
  U <- nt.nor[,tests[,i]]
  plot(U[,1], U[,2], xlab=tests[1,i], ylab=tests[2,i])
}

dev.off()

##### plots for the first 5 pairs 
pdf(file = paste0(dirl.p,"DA_Nutrition_1_b1000.pdf"),
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
# calcium vs iron 
par(mar = c(0.5, 5, 0.2, 0.5))
U <- nt.nor[,tests[,1]]




f.y <- paste0(dirl.p,"data_matrices/CalciumIron_sym.txt")
f.p <- paste0(dirl.p,"p_values/CalciumIron_sym.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}

tfplot(y, x=t, x.label="", xaxt="n", yaxt="n", ylim=c(-0.06, 0.06),
       y.label="calcium / iron", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

par(mar = c(0.5, 4.3, 0.2, 1.2))

f.y <- paste0(dirl.p,"data_matrices/CalciumIron_rsym.txt")
f.p <- paste0(dirl.p,"p_values/CalciumIron_rsym.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", yaxt="n", ylim=c(-0.06, 0.06), p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))


f.y <- paste0(dirl.p,"data_matrices/CalciumIron_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/CalciumIron_jsym_1.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)


f.y <- paste0(dirl.p,"data_matrices/CalciumIron_jsym_2.txt")
f.p <- paste0(dirl.p,"p_values/CalciumIron_jsym_2.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)


####################
# calcium vs protein 
par(mar = c(0.5, 5, 1.5, 0.5))
U <- nt.nor[,tests[,2]]


f.y <- paste0(dirl.p,"data_matrices/CalciumProtein_sym.txt")
f.p <- paste0(dirl.p,"p_values/CalciumProtein_sym.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", yaxt="n", ylim=c(-0.06, 0.06),
       y.label="calcium / protein", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

par(mar = c(0.5, 4.3, 1.5, 1.2))

f.y <- paste0(dirl.p,"data_matrices/CalciumProtein_rsym.txt")
f.p <- paste0(dirl.p,"p_values/CalciumProtein_rsym.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", yaxt="n", ylim=c(-0.06, 0.06), p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))


f.y <- paste0(dirl.p,"data_matrices/CalciumProtein_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/CalciumProtein_jsym_1.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)


f.y <- paste0(dirl.p,"data_matrices/CalciumProtein_jsym_2.txt")
f.p <- paste0(dirl.p,"p_values/CalciumProtein_jsym_2.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

######################
# calcium vs vitamin A
par(mar = c(0.5, 5, 1.5, 0.5))
U <- nt.nor[,tests[,3]]

f.y <- paste0(dirl.p,"data_matrices/CalciumVitaminA_sym.txt")
f.p <- paste0(dirl.p,"p_values/CalciumVitaminA_sym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", yaxt="n", ylim=c(-0.06, 0.06),
       y.label="calcium / vitamin A", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

par(mar = c(0.5, 4.3, 1.5, 1.2))


f.y <- paste0(dirl.p,"data_matrices/CalciumVitaminA_rsym.txt")
f.p <- paste0(dirl.p,"p_values/CalciumVitaminA_rsym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", yaxt="n", ylim=c(-0.06, 0.06), p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))


f.y <- paste0(dirl.p,"data_matrices/CalciumVitaminA_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/CalciumVitaminA_jsym_1.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)


f.y <- paste0(dirl.p,"data_matrices/CalciumVitaminA_jsym_2.txt")
f.p <- paste0(dirl.p,"p_values/CalciumVitaminA_jsym_2.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

######################
# calcium vs vitamin C
par(mar = c(0.5, 5, 1.5, 0.5))
U <- nt.nor[,tests[,4]]


f.y <- paste0(dirl.p,"data_matrices/CalciumVitaminC_sym.txt")
f.p <- paste0(dirl.p,"p_values/CalciumVitaminC_sym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", yaxt="n", ylim=c(-0.06, 0.06),
       y.label="calcium / vitamin C", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

par(mar = c(0.5, 4.3, 1.5, 1.2))

f.y <- paste0(dirl.p,"data_matrices/CalciumVitaminC_rsym.txt")
f.p <- paste0(dirl.p,"p_values/CalciumVitaminC_rsym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="v", xaxt="n", y.label="", yaxt="n", ylim=c(-0.06, 0.06), p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))


f.y <- paste0(dirl.p,"data_matrices/CalciumVitaminC_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/CalciumVitaminC_jsym_1.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="v", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)


f.y <- paste0(dirl.p,"data_matrices/CalciumVitaminC_jsym_2.txt")
f.p <- paste0(dirl.p,"p_values/CalciumVitaminC_jsym_2.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="v", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

#################
# iron vs protein
par(mar = c(4.4, 5, 1.5, 0.5))
U <- nt.nor[,tests[,5]]


f.y <- paste0(dirl.p,"data_matrices/IronProtein_sym.txt")
f.p <- paste0(dirl.p,"p_values/IronProtein_sym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="v", yaxt="n", ylim=c(-0.06, 0.06),
       y.label="iron / protein", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

par(mar = c(4.4, 4.3, 1.5, 1.2))

f.y <- paste0(dirl.p,"data_matrices/IronProtein_rsym.txt")
f.p <- paste0(dirl.p,"p_values/IronProtien_rsym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="v", y.label="", yaxt="n", ylim=c(-0.06, 0.06), p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))


f.y <- paste0(dirl.p,"data_matrices/IronProtein_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/IronProtein_jsym_1.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="v", y.label="", p.value=p.val, cex.lab=1.8)


f.y <- paste0(dirl.p,"data_matrices/IronProtein_jsym_2.txt")
f.p <- paste0(dirl.p,"p_values/IronProtein_jsym_2.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="v", y.label="", p.value=p.val, cex.lab=1.8)

dev.off()




#################################
######### second part of the plot 
pdf(file = paste0(dirl.p,"DA_Nutrition_2_b1000.pdf"),
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

###################
# iron vs vitamin A
par(mar = c(0.5, 5, 0.2, 0.5))
U <- nt.nor[,tests[,6]]

f.y <- paste0(dirl.p,"data_matrices/IronVitaminA_sym.txt")
f.p <- paste0(dirl.p,"p_values/IronVitaminA_sym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n",
       y.label="iron / vitamin A", p.value=p.val, cex.lab=1.8)

par(mar = c(0.5, 4.3, 0.2, 1.2))

f.y <- paste0(dirl.p,"data_matrices/IronVitaminA_rsym.txt")
f.p <- paste0(dirl.p,"p_values/IronVitaminA_rsym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", ylim=c(-0.06, 0.06), 
       y.label="", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))


f.y <- paste0(dirl.p,"data_matrices/IronVitaminA_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/IronVitaminA_jsym_1.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)


f.y <- paste0(dirl.p,"data_matrices/IronVitaminA_jsym_2.txt")
f.p <- paste0(dirl.p,"p_values/IronVitaminA_jsym_2.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

###################
# iron vs vitamin C
par(mar = c(0.5, 5, 1.5, 0.5))
U <- nt.nor[,tests[,7]]

f.y <- paste0(dirl.p,"data_matrices/IronVitaminC_sym.txt")
f.p <- paste0(dirl.p,"p_values/IronVitaminC_sym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n",
       y.label="iron / vitamin C", p.value=p.val, cex.lab=1.8)

par(mar = c(0.5, 4.3, 1.5, 1.2))

f.y <- paste0(dirl.p,"data_matrices/IronVitaminC_rsym.txt")
f.p <- paste0(dirl.p,"p_values/IronVitaminC_rsym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

f.y <- paste0(dirl.p,"data_matrices/IronVitaminC_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/IronVitaminC_jsym_1.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)


f.y <- paste0(dirl.p,"data_matrices/IronVitaminC_jsym_2.txt")
f.p <- paste0(dirl.p,"p_values/IronVitaminC_jsym_2.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

######################
# protein vs vitamin A
par(mar = c(0.5, 5, 1.5, 0.5))
U <- nt.nor[,tests[,8]]

f.y <- paste0(dirl.p,"data_matrices/ProteinVitaminA_sym.txt")
f.p <- paste0(dirl.p,"p_values/ProteinVitaminA_sym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n",
       y.label="protein / vitamin A", p.value=p.val, cex.lab=1.8)

par(mar = c(0.5, 4.3, 1.5, 1.2))

f.y <- paste0(dirl.p,"data_matrices/ProteinVitaminA_rsym.txt")
f.p <- paste0(dirl.p,"p_values/ProteinVitaminA_rsym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)


f.y <- paste0(dirl.p,"data_matrices/ProteinVitaminA_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/ProteinVitaminA_jsym_1.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)


f.y <- paste0(dirl.p,"data_matrices/ProteinVitaminA_jsym_2.txt")
f.p <- paste0(dirl.p,"p_values/ProteinVitaminA_jsym_2.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

######################
# protein vs vitamin C
par(mar = c(0.5, 5, 1.5, 0.5))
U <- nt.nor[,tests[,9]]

f.y <- paste0(dirl.p,"data_matrices/ProteinVitaminC_sym.txt")
f.p <- paste0(dirl.p,"p_values/ProteinVitaminC_sym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", yaxt="n", ylim=c(-0.06, 0.06),
       y.label="protein / vitamin C", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

par(mar = c(0.5, 4.3, 1.5, 1.2))

f.y <- paste0(dirl.p,"data_matrices/ProteinVitaminC_rsym.txt")
f.p <- paste0(dirl.p,"p_values/ProteinVitaminC_rsym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

f.y <- paste0(dirl.p,"data_matrices/ProteinVitaminC_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/ProteinVitaminC_jsym_1.txt")
if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

f.y <- paste0(dirl.p,"data_matrices/ProteinVitaminC_jsym_2.txt")
f.p <- paste0(dirl.p,"p_values/ProteinVitaminC_jsym_2.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="", xaxt="n", y.label="", p.value=p.val, cex.lab=1.8)

########################
# vitamin A vs vitamin C
par(mar = c(4.4, 5, 1.5, 0.5))
U <- nt.nor[,tests[,10]]


f.y <- paste0(dirl.p,"data_matrices/VitaminsAC_sym.txt")
f.p <- paste0(dirl.p,"p_values/VitaminsAC_sym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.sym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="v", yaxt="n", ylim=c(-0.06, 0.06),
       y.label="vitamin A / vitamin C", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

par(mar = c(4.4, 4.3, 1.5, 1.2))

f.y <- paste0(dirl.p,"data_matrices/VitaminsAC_rsym.txt")
f.p <- paste0(dirl.p,"p_values/VitaminsAC_rsym.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.r(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.rsym(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="v", ylim=c(-0.06, 0.06), y.label="", yaxt="n", p.value=p.val, cex.lab=1.8)
axis(side = 2, at=c(-0.06, -0.02, 0.02, 0.06))

f.y <- paste0(dirl.p,"data_matrices/VitaminsAC_jsym_1.txt")
f.p <- paste0(dirl.p,"p_values/VitaminsAC_jsym_1.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j1(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.1(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}
tfplot(y, x=t, x.label="v", y.label="", p.value=p.val, cex.lab=1.8)

f.y <- paste0(dirl.p,"data_matrices/VitaminsAC_jsym_2.txt")



f.p <- paste0(dirl.p,"p_values/VitaminsAC_jsym_2.txt")

if (file.exists(f.y)){
  y <- as.matrix(read.table(f.y, header=F, sep=" "))
  p.val <- as.numeric(read.table(f.p, header=F, sep=" "))
} else {
  y <- get.y.j2(U, n.curves=n.curves, t)
  write.table(y, file=f.y, row.names=F, col.names=F)
  p.val <- test.jsym.2(U, y, N.b = N.b, n.curves=n.curves)
  write.table(p.val, file=f.p, row.names=F, col.names=F)
}


















tfplot(y, x=t, x.label="v", y.label="", p.value=p.val, cex.lab=1.8)

dev.off()






