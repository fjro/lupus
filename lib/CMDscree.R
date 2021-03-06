#' Function to draw a screeplot of cumulative
#' eigenvalues in classical metric scaling
#' 
#' @param mydata A matrix
#' @param raw T if raw eigenvalues should be used
#' @param abs
#' @param maxcomp The max number of components
CMDscreeplot <- function(mydata,
                         raw = T,
                         abs = T,
                         maxcomp = 10) {
  n = min(dim(mydata)) - 1
  x <- cmdscale(dist(mydata), k = n, eig = TRUE)

  k <- min(n, maxcomp)
  w <- c(0:k)
  
  if (raw == T)
    y <- x$eig
  
  
  if (raw == F & abs == T)
    y <- abs(x$eig)
  if (raw == F & abs == F)
    y <- x$eig * x$eig
  y <- c(0, y)
  z <- 100 * cumsum(y) / sum(y)
  
  if (raw == T) {
    plot(
      w,
      z[1:(k + 1)],
      type = "l",
      xlab = "number of dimensions",
      cex.main = 1.5,
      lwd = 3,
      col = "red",
      ylim = c(0, 100),
      ylab = "cumulative percentage of total eigenvlaues",
      main = "Scree plot of eigenvalues",
      xaxt = "n",
      yaxt = "n"
    )
  }
  
  if (raw == F &
      abs == T) {
    plot(
      w,
      z[1:(k + 1)],
      type = "l",
      xlab = "number of dimensions",
      cex.main = 1.5,
      lwd = 3,
      col = "red",
      ylim = c(0, 100),
      ylab = "cumulative percentage of total absolutes eigenvalues",
      main = "Scree plot of absolute values of eigenvalues",
      xaxt = "n",
      yaxt = "n"
    )
  }
  
  if (raw == F &
      abs == F) {
    plot(
      w,
      z[1:(k + 1)],
      type = "l",
      xlab = "number of dimensions",
      cex.main = 1.5,
      lwd = 3,
      col = "red",
      ylim = c(0, 100),
      ylab = "cumulative percentage of total variance",
      main = "Scree plot of squared eigenvalues",
      xaxt = "n",
      yaxt = "n"
    )
  }
  
  axis(1, at = w, lwd = 2)
  axis(2, at = c(0, 20, 40, 60, 80, 100), lwd = 2)
  abline(
    a = 100,
    b = 0,
    lwd = 2,
    lty = "dashed",
    col = "orange"
  )
  text(
    w,
    z[1:(k + 1)],
    labels = w,
    cex = 0.8,
    adj = c(1.2, -.1),
    col = "blue"
  )
}


# Examples of calls to it are
#default choice using raw eigenvalues and maximum of ten
CMDscreeplot(genes)



