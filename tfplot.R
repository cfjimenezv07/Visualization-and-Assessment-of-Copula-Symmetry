###########################################################################
# Install packages
###########################################################################
library(fdaoutlier)

###########################################################################
# Set a working directory
###########################################################################

tfplot <- function(fit, x = NULL, ylim = NULL,
                   x.resolution = 6, y.resolution = 6,
                   x.label = "u", y.label = "function",
                   main = "", p.value = NA,
                   p.value.resolution = 0.001,
                   p.value.show = TRUE,
                   p.value.show.pos = "topright",
                   p.value.show.size = 1, alpha = 0.05, ...) {
  #' Visualization of test functions
  #'
  #' This function visualizes the obtained test functions.
  #'
  #' @param fit Matrix. The obtained test functions. Each column is a realization of the test functions.
  #' @param x Vector. The coordinates in the x-axis to be displayed.
  #' The length should be equal to the number of rows of the input \emph{fit}.
  #' @param ylim Vector of two numeric numbers. The range of the y-axis to be displayed.
  #' @param x.resolution Interger. Number of intermediate points for one temporal lag along the x-axis to draw density
  #' @param y.resolution Interger. Number of intermediate points along the y-axis to draw density
  #' @param x.label String. The label on the x-axis.
  #' @param y.label String. The label on the y-axis.
  #' @param main String. The plot title to be displayed.
  #' @param p.value NA or numeric. The obtained p-value to visualize the test functions.
  #'          - If NA, the central region of the test functions is shown in magenta.
  #'
  #'          - If a numeric value other than NA, the p-value is displayed and
  #'          the central region of the test functions is shown in red when \eqn{p.value < alpha} (the test is rejected)
  #'          or green when \eqn{p.value >= alpha} (the test is not rejected).
  #' @param p.value.resolution Numeric. The smallest precise p-value in the bootstrap precedure, i.e, the reciprocal
  #'                           of the number of bootstraps. When \eqn{p.value < p.value.resolution} such as \eqn{p.value = 0},
  #'                           the plot only displays \eqn{p < p.value.resolution}.
  #' @param p.value.show.pos String. The position of the p-value displayed in the plot.
  #' @param p.value.show.size Numeric. The size of the text showing the p-value.
  #' @param alpha Numeric. The significance level for the test.
  #' @param ... Other arguments to be passed to the R function \emph{plot}.
  
  #---------- choose color
  if (is.na(p.value)) {
    col <- "magenta"
  } else {
    col <- ifelse(p.value < alpha, "red", "green3")
  }
  
  #---------- get dimension
  tp <- dim(fit)[1]
  n <- dim(fit)[2]
  if (is.null(x)) x <- 1:tp
  
  
  #---------- compute functional data depth
  #depth <- fMBD_not_scaled(fit)
  depth <- modified_band_depth(t(fit))
  dp_s <- sort(depth, decreasing = TRUE)
  index <- order(depth, decreasing = TRUE)
  med <- depth == max(depth)
  
  #---------- compute central area
  m <- ceiling(n * 0.5)
  center <- fit[, index[1:m]]
  out <- fit[, index[(m + 1):n]]
  inf <- apply(center, 1, min)
  sup <- apply(center, 1, max)
  
  factor <- 1.5
  dist <- factor * (sup - inf)
  upper <- sup + dist
  lower <- inf - dist
  outly <- (fit < lower) + (fit > upper)
  outcol <- colSums(outly)
  remove <- (outcol > 0)
  colum <- 1:n
  outpoint <- colum[remove == 1]
  out <- fit[, remove]
  woout <- fit
  good <- woout[, (remove == 0), drop = FALSE]
  maxcurve <- apply(good, 1, max)
  mincurve <- apply(good, 1, min)
  
  #---------- decide plot limit
  if (is.null(ylim)) {
    y_max <- max(maxcurve)
    y_min <- min(mincurve)
    offset <- (y_max - y_min) / 10
    ylim <- max(abs(y_min - offset), abs(y_max + offset))
    ylim <- c(-ylim, ylim)
  }
  
  #---------- draw background
  plot(x, rep(0, length(x)),
       lty = 1, col = "white", type = "l", ylim=1.2*ylim, yaxt="n",
       xlab = x.label, ylab = y.label, ...
  )
  axis(side = 2, at=c(round(ylim[1], 2), round(ylim[1]/3, 2), round(ylim[2]/3, 2), round(ylim[2], 2)))
  title(main, line = 1, cex.main = list(...)$cex.main)
  
  if (!is.na(p.value) && p.value.show) {
    if (p.value >= p.value.resolution) {
      legend(p.value.show.pos, paste0(
        "p=",
        sprintf(paste0("%.", floor(log(p.value.resolution) / log(0.1)), "f"), p.value)
      ),
      cex = p.value.show.size, x.intersp = 0
      )
    } else {
      legend(p.value.show.pos, paste("p<", p.value.resolution, sep = ""),
             cex = p.value.show.size, x.intersp = 0
      )
    }
  }
  
  #---------- draw density
  ## Initialize density matrix
  density.x <- seq(x[1], x[tp], length.out = x.resolution * (tp - 1) + 1)
  density.nx <- length(density.x)
  density.y <- matrix(NA, density.nx, y.resolution)
  density.mat <- matrix(0, density.nx - 1, y.resolution - 1)
  
  for (u in 1:(tp - 1))
  {
    density.y[1:x.resolution + (u - 1) * x.resolution, 1] <-
      seq(inf[u], inf[u + 1], length.out = x.resolution + 1)[1:x.resolution]
    density.y[1:x.resolution + (u - 1) * x.resolution, y.resolution] <-
      seq(sup[u], sup[u + 1], length.out = x.resolution + 1)[1:x.resolution]
  }
  density.y[density.nx, 1] <- inf[tp]
  density.y[density.nx, y.resolution] <- sup[tp]
  
  for (i in 1:density.nx) {
    density.y[i, ] <- seq(density.y[i, 1], density.y[i, y.resolution], length.out = y.resolution)
  }
  
  ## Compute density matrix
  for (c in 1:dim(center)[2])
  {
    curve <- center[, c]
    for (u in 1:(tp - 1))
    {
      curve.point <- seq(curve[u], curve[u + 1], length.out = x.resolution + 1)[1:x.resolution]
      for (i in 1:x.resolution)
      {
        index <- which(density.y[i + (u - 1) * x.resolution, ] <= curve.point[i])
        if (length(index) > 0) {
          index <- min(index[length(index)], y.resolution - 1)
          density.mat[i + (u - 1) * x.resolution, index] <- density.mat[i + (u - 1) * x.resolution, index] + 1
        }
      }
    }
  }
  
  draw_density(density.mat, density.x, density.y, col = col)
  
  #---------- draw central region border
  lines(x, sup, col = "blue", lwd = 2, lty = 1)
  lines(x, inf, col = "blue", lwd = 2, lty = 1)
  
  lines(c(x[1], x[1]), c(inf[1], sup[1]), col = "blue", lwd = 2, lty = 1)
  lines(c(x[tp], x[tp]), c(inf[tp], sup[tp]), col = "blue", lwd = 2, lty = 1)
  
  #---------- draw max and min curve
  lines(x, maxcurve, col = "blue", lwd = 2, lty = 1)
  lines(x, mincurve, col = "blue", lwd = 2, lty = 1)
  
  #---------- draw vertical bar
  barval <- (x[1] + x[tp]) / 2
  bar <- which(sort(c(x, barval)) == barval)[1]
  
  lines(c(x[bar], x[bar]), c(maxcurve[bar], sup[bar]), col = 4, lwd = 2, lty = 1)
  lines(c(x[bar], x[bar]), c(mincurve[bar], inf[bar]), col = 4, lwd = 2, lty = 1)
  
  #---------- draw zero curve
  lines(c(x[1], x[tp]), c(0, 0), col = "black", lwd = 4, lty = 3)
  # abline(h = 0, col = "black", lwd = 4, lty = 3)
  
  # return(density.mat)
}

draw_density <- function(density.mat, x, y, col = "magenta", ncolor = 32) {
  offset <- 2 # Add offset to avoid colors very close to white
  pal <- colorRampPalette(c("white", col))(ncolor + offset)
  
  # Transform densityMat to [1, ncolor]
  max.value <- max(density.mat)
  min.value <- min(density.mat)
  density.mat <- (density.mat - min.value) / (max.value - min.value) # [0,1]
  density.mat <- floor(density.mat * (ncolor - 1) + 1) # [1, ncolor]
  
  # Draw polygon
  nx <- length(x)
  ny <- dim(y)[2]
  for (i in 1:(nx - 1)) {
    for (j in 1:(ny - 1))
    {
      polygon(c(x[i], x[i + 1], x[i + 1], x[i]),
              c(y[i, j], y[i + 1, j], y[i + 1, j + 1], y[i, j + 1]),
              col = pal[density.mat[i, j] + offset],
              border = pal[density.mat[i, j] + offset]
      )
    }
  }
}
