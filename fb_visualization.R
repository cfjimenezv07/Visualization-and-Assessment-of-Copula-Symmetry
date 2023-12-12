###########################################################################
# Install packages
###########################################################################
library(fda)

###########################################################################
# Set a working directory
###########################################################################




###########################################################################
fb.visualization <- function(data, t, x.lab, y.lab, title, y.lim = c(-0.1, 0.1), alpha=0.05, p.val=NULL, 
                             p.val.pos="bottomleft", p.val.size = 1, p.val.resolution = 0.001, p.val.show=TRUE){
  
  # choose the color displayed in functional boxplot
  isNA <- is.na(p.val)
  if (isNA) {
    col <- "magenta"
  } else {
    col <- ifelse(p.val <= alpha, "red", "green3")
  }
  
  # plot functional boxplot 
  p <- dim(data)[1]
  
  # use fda::fbplot 
  fbplot(data, x = t, color=pal, xlim=c(0, 1), ylim=y.lim, main=title, cex.lab=1.5, xlab=x.lab, ylab=y.lab)
  lines(c(0, 1), c(0, 0), col = "black", lwd = 4, lty = "dotted")
  if ((!isNA) && p.val.show ){
    if (p.val > 0.001) {
      legend(p.val.pos, 
             paste0("p=", sprintf(paste0("%.", floor(log(p.val.resolution) / log(0.1)), "f"), p.val)),
             cex = p.val.size, 
             x.intersp = 0
      )
    } else {
      legend(p.val.pos, 
             paste("p<", p.val.resolution, sep = ""),
             cex = p.val.size, x.intersp = 0
      )
    }
  }
}

###########################################################################