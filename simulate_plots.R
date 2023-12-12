
###########################################################################
# Install packages
###########################################################################
library(copula)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(egg)
###########################################################################
# Set a working directory
###########################################################################
# wd.

# define a folder where plots must be saved
dirl.p <- "~/My Drive/Spring 2023/STAT 397/PhD project 3/Plots_paper/"
###########################################################################
# Load Auxiliary files
###########################################################################
source("copula_simulation.R")

################
# Clayton copula 
theta <- BiCopTau2Par(3, tau)
cc <- claytonCopula(theta)
U <- rCopula(1000, cc)
ec <- empCopula(U)
U <- as.data.frame(U)

p.cl <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.text.y=element_text(angle=90,size=8),
  axis.title=element_text(size=12), 
  plot.title = element_text(hjust = 0.5))+ggtitle("Original")+
  xlab("")+ylab("Clayton (S)")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.cl <- set_panel_size(p.cl, width=unit(4, "cm"), height=unit(4, "cm"))

U.sym <- simulate.sym(1000, ec)
U <- as.data.frame(U.sym)
p.cl.sym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12), 
  plot.title = element_text(hjust = 0.5))+ggtitle("Symmetry (S)")+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.cl.sym <- set_panel_size(p.cl.sym, width=unit(4, "cm"), height=unit(4, "cm"))

U.rsym <- simulate.rsym(1000, ec)
U <- as.data.frame(U.rsym)
p.cl.rsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12), 
  plot.title = element_text(hjust = 0.5))+ggtitle("Radial Symmetry (R)")+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.cl.rsym <- set_panel_size(p.cl.rsym, width=unit(4, "cm"), height=unit(4, "cm"))

U.jsym <- simulate.jsym(1000, ec)
U <- as.data.frame(U.jsym)
p.cl.jsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12), 
  plot.title=element_text(hjust = 0.5))+ggtitle("Joint Symmetry (J)")+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.cl.jsym <- set_panel_size(p.cl.jsym, width=unit(4, "cm"), height=unit(4, "cm"))

################
# Frank 
theta <- BiCopTau2Par(5, tau)
fc <- frankCopula(theta)
U <- rCopula(1000, fc)
ec <- empCopula(U)
U <- as.data.frame(U)

p.fr <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.y = element_text(angle=90,size=8),
  axis.title=element_text(size=12), 
  plot.title = element_blank())+
  xlab("")+ylab("Frank (S, R)")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.fr <- set_panel_size(p.fr, width=unit(4, "cm"), height=unit(4, "cm"))

U.sym <- simulate.sym(1000, ec)
U <- as.data.frame(U.sym)
p.fr.sym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.fr.sym <- set_panel_size(p.fr.sym, width=unit(4, "cm"), height=unit(4, "cm"))

U.rsym <- simulate.rsym(1000, ec)
U <- as.data.frame(U.rsym)
p.fr.rsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12), 
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.fr.rsym <- set_panel_size(p.fr.rsym, width=unit(4, "cm"), height=unit(4, "cm"))

U.jsym <- simulate.jsym(1000, ec)
U <- as.data.frame(U.jsym)
p.fr.jsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.fr.jsym <- set_panel_size(p.fr.jsym, width=unit(4, "cm"), height=unit(4, "cm"))


###################
### Gaussian copula 
theta <- BiCopTau2Par(1, tau)
nc <- normalCopula(theta)
U <- rCopula(1000, nc)
ec <- empCopula(U)
U <- as.data.frame(U)

p.ga <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.y = element_text(angle=90,size=8),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("Gaussian (S, R)")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.ga <- set_panel_size(p.ga, width=unit(4, "cm"), height=unit(4, "cm"))

U.sym <- simulate.sym(1000, ec)
U <- as.data.frame(U.sym)
p.ga.sym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.ga.sym <- set_panel_size(p.ga.sym, width=unit(4, "cm"), height=unit(4, "cm"))

U.rsym <- simulate.rsym(1000, ec)
U <- as.data.frame(U.rsym)
p.ga.rsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.ga.rsym <- set_panel_size(p.ga.rsym, width=unit(4, "cm"), height=unit(4, "cm"))

U.jsym <- simulate.jsym(1000, ec)
U <- as.data.frame(U.jsym)
p.ga.jsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.ga.jsym <- set_panel_size(p.ga.jsym, width=unit(4, "cm"), height=unit(4, "cm"))


###################
### Independent 
ic <- indepCopula()
U <- rCopula(1000, ic)
ec <- empCopula(U)
U <- as.data.frame(U)

p.ind <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_text(size=8),
  axis.text.y=element_text(angle=90, size=8),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab(expression(U[1]))+ylab("Independent (S, R, J)")+
  scale_x_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))+
  scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.ind <- set_panel_size(p.ind, width=unit(4, "cm"), height=unit(4, "cm"))

U.sym <- simulate.sym(1000, ec)
U <- as.data.frame(U.sym)
p.ind.sym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_text(size=8),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab(expression(U[1]))+ylab("")+
  scale_x_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))+
  scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.ind.sym <- set_panel_size(p.ind.sym, width=unit(4, "cm"), height=unit(4, "cm"))

U.rsym <- simulate.rsym(1000, ec)
U <- as.data.frame(U.rsym)
p.ind.rsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_text(size=8),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab(expression(U[1]))+ylab("")+
  scale_x_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))+
  scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.ind.rsym <- set_panel_size(p.ind.rsym, width=unit(4, "cm"), height=unit(4, "cm"))

U.jsym <- simulate.jsym(1000, ec)
U <- as.data.frame(U.jsym)
p.ind.jsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_text(size=8),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab(expression(U[1]))+ylab("")+
  scale_x_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))+
  scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.ind.jsym <- set_panel_size(p.ind.jsym, width=unit(4, "cm"), height=unit(4, "cm"))

##################
### build plots 
p.cl.jsym <- annotate_figure(p.cl.jsym, right=text_grob(expression(paste(U[2], "  ")), rot=-90))
p.fr.jsym <- annotate_figure(p.fr.jsym, right=text_grob(expression(paste(U[2], "  ")), rot=-90))
p.ga.jsym <- annotate_figure(p.ga.jsym, right=text_grob(expression(paste(U[2], "  ")), rot=-90))
p.ind.jsym <- annotate_figure(p.ind.jsym, right=text_grob(expression(paste(U[2], "  ")), rot=-90))

p.ori <- grid.arrange(p.cl, p.fr, p.ga, p.ind, ncol=1)
p.sym <- grid.arrange(p.cl.sym, p.fr.sym, p.ga.sym, p.ind.sym, ncol=1)
p.rsym <- grid.arrange(p.cl.rsym, p.fr.rsym, p.ga.rsym, p.ind.rsym, ncol=1)
p.jsym <- grid.arrange(p.cl.jsym, p.fr.jsym, p.ga.jsym, p.ind.jsym, ncol=1)


pdf(file = paste0(dirl.p,"simulated_copulas_1.pdf"),
    width = 9, 
    height = 10)
grid.arrange(p.ori, p.sym, p.rsym, p.jsym, nrow=1)
dev.off()


################# 
# plot 2 
# Gumbel copula 
theta <- BiCopTau2Par(4, tau)
gc <- gumbelCopula(theta)
U <- rCopula(1000, gc)
ec <- empCopula(U)
U <- as.data.frame(U)

p.gb <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.text.y=element_text(angle=90,size=8),
  axis.title=element_text(size=12), 
  plot.title = element_text(hjust = 0.5))+ggtitle("Original")+
  xlab("")+ylab("Gumbel (S)")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.gb <- set_panel_size(p.gb, width=unit(4, "cm"), height=unit(4, "cm"))

U.sym <- simulate.sym(1000, ec)
U <- as.data.frame(U.sym)
p.gb.sym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12), 
  plot.title = element_text(hjust = 0.5))+ggtitle("Symmetry (S)")+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.cl.sym <- set_panel_size(p.gb.sym, width=unit(4, "cm"), height=unit(4, "cm"))

U.rsym <- simulate.rsym(1000, ec)
U <- as.data.frame(U.rsym)
p.gb.rsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border = element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12), 
  plot.title = element_text(hjust = 0.5))+ggtitle("Radial Symmetry (R)")+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.gb.rsym <- set_panel_size(p.gb.rsym, width=unit(4, "cm"), height=unit(4, "cm"))

U.jsym <- simulate.jsym(1000, ec)
U <- as.data.frame(U.jsym)
p.gb.jsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12), 
  plot.title=element_text(hjust = 0.5))+ggtitle("Joint Symmetry (J)")+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.gb.jsym <- set_panel_size(p.gb.jsym, width=unit(4, "cm"), height=unit(4, "cm"))

################
# Marshall-Olkin
moc <- moCopula(c(0.55, 0.85))
U <- rCopula(1000, moc)
ec <- empCopula(U)
U <- as.data.frame(U)

p.mo <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.y = element_text(angle=90,size=8),
  axis.title=element_text(size=12), 
  plot.title = element_blank())+
  xlab("")+ylab("Marshall-Olkin (AS)")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.mo <- set_panel_size(p.mo, width=unit(4, "cm"), height=unit(4, "cm"))

U.sym <- simulate.sym(1000, ec)
U <- as.data.frame(U.sym)
p.mo.sym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.mo.sym <- set_panel_size(p.mo.sym, width=unit(4, "cm"), height=unit(4, "cm"))

U.rsym <- simulate.rsym(1000, ec)
U <- as.data.frame(U.rsym)
p.mo.rsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12), 
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.mo.rsym <- set_panel_size(p.mo.rsym, width=unit(4, "cm"), height=unit(4, "cm"))

U.jsym <- simulate.jsym(1000, ec)
U <- as.data.frame(U.jsym)
p.mo.jsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.mo.jsym <- set_panel_size(p.mo.jsym, width=unit(4, "cm"), height=unit(4, "cm"))


###################
### Tawn Type 1 
tct1 <- tawnT1Copula(c(4.28, 0.6))
U <- rCopula(1000, tct1)
ec <- empCopula(U)
U <- as.data.frame(U)

p.tt1 <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.y = element_text(angle=90,size=8),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("Tawn T1 (AS)")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.tt1 <- set_panel_size(p.tt1, width=unit(4, "cm"), height=unit(4, "cm"))

U.sym <- simulate.sym(1000, ec)
U <- as.data.frame(U.sym)
p.tt1.sym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.tt1.sym <- set_panel_size(p.tt1.sym, width=unit(4, "cm"), height=unit(4, "cm"))

U.rsym <- simulate.rsym(1000, ec)
U <- as.data.frame(U.rsym)
p.tt1.rsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.tt1.rsym <- set_panel_size(p.tt1.rsym, width=unit(4, "cm"), height=unit(4, "cm"))

U.jsym <- simulate.jsym(1000, ec)
U <- as.data.frame(U.jsym)
p.tt1.jsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab("")+ylab("")+scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.tt1.jsym <- set_panel_size(p.tt1.jsym, width=unit(4, "cm"), height=unit(4, "cm"))


###################
# Tawn Type 2 
tct2 <- tawnT2Copula(c(4.28, 0.6))
U <- rCopula(1000, tct2)
ec <- empCopula(U)
U <- as.data.frame(U)

p.tt2 <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_text(size=8),
  axis.text.y=element_text(angle=90, size=8),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab(expression(U[1]))+ylab("Tawn T2 (AS)")+
  scale_x_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))+
  scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.tt2 <- set_panel_size(p.tt2, width=unit(4, "cm"), height=unit(4, "cm"))

U.sym <- simulate.sym(1000, ec)
U <- as.data.frame(U.sym)
p.tt2.sym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_text(size=8),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab(expression(U[1]))+ylab("")+
  scale_x_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))+
  scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.tt2.sym <- set_panel_size(p.tt2.sym, width=unit(4, "cm"), height=unit(4, "cm"))

U.rsym <- simulate.rsym(1000, ec)
U <- as.data.frame(U.rsym)
p.tt2.rsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_text(size=8),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab(expression(U[1]))+ylab("")+
  scale_x_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))+
  scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.tt2.rsym <- set_panel_size(p.tt2.rsym, width=unit(4, "cm"), height=unit(4, "cm"))

U.jsym <- simulate.jsym(1000, ec)
U <- as.data.frame(U.jsym)
p.tt2.jsym <- ggplot(data=U, aes(x = V1, y = V2))+geom_point(shape=1, size=0.8)+theme_bw()+theme(
  panel.border =element_rect(colour = "black"),
  panel.grid = element_blank(),
  #plot.margin = margin(1.5, -1, 1, -1, unit="cm"),
  axis.text.x=element_text(size=8),
  axis.ticks.y = element_blank(),
  axis.text.y=element_blank(),
  axis.title=element_text(size=12),
  plot.title = element_blank())+
  xlab(expression(U[1]))+ylab("")+
  scale_x_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))+
  scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
p.tt2.jsym <- set_panel_size(p.tt2.jsym, width=unit(4, "cm"), height=unit(4, "cm"))

##################
### build plots 
p.gb.jsym <- annotate_figure(p.gb.jsym, right=text_grob(expression(paste(U[2], "  ")), rot=-90))
p.mo.jsym <- annotate_figure(p.mo.jsym, right=text_grob(expression(paste(U[2], "  ")), rot=-90))
p.tt1.jsym <- annotate_figure(p.tt1.jsym, right=text_grob(expression(paste(U[2], "  ")), rot=-90))
p.tt2.jsym <- annotate_figure(p.tt2.jsym, right=text_grob(expression(paste(U[2], "  ")), rot=-90))

p.ori <- grid.arrange(p.gb, p.mo, p.tt1, p.tt2, ncol=1)
p.sym <- grid.arrange(p.gb.sym, p.mo.sym, p.tt1.sym, p.tt2.sym, ncol=1)
p.rsym <- grid.arrange(p.gb.rsym, p.mo.rsym, p.tt1.rsym, p.tt2.rsym, ncol=1)
p.jsym <- grid.arrange(p.gb.jsym, p.mo.jsym, p.tt1.jsym, p.tt2.jsym, ncol=1)


pdf(file = paste0(dirl.p,"simulated_copulas_2.pdf"),
    width = 9, 
    height = 10)
grid.arrange(p.ori, p.sym, p.rsym, p.jsym, nrow=1)
dev.off()

