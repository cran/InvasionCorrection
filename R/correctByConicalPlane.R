#' Correct invasion data by conical plane
#'
#' @description Correct z-component of a 3D collagen invasion essay.
#' The correction is achieved under the assumption that non-migrating cells of the essay approximately
#' form a quadratic flow profile due to frictional effects, compare law of Hagen-Poiseuille for flow in a tube.
#'
#' @param filename Name of data file in csv format. It should contain columns "Pos_X", "Pos_Y" and "Pos_Z".
#' @param nrfits Numeric, Number of randomly chosen starting points for the optimization. Choose lower values for speeding up computational time.
#' Choose higher values for more reliable optimization results.
#' @param threshold Numeric, A threshold for counting cells as being invaded or not. When cells move towards negative z-direction, threshold should be negative.
#' @param plot Boole, if TRUE exemplary 3D plots before and after the correction are plotted
#' @param write_csv, if TRUE resulting corrected values are saved as a csv file
#'
#' @return Data.frame containing input positions, corrected z-positions as well as number and percentage of invaded cells.
#' @import lattice
#' @import stats
#' @import utils
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}

correctByConicalPlane <- function(filename, nrfits = 1000, threshold = -30, plot = FALSE, write_csv = TRUE){
  xyz_data <- read.csv(filename, sep=";", dec=",")
  xdata <- xyz_data$Pos_X
  ydata <- xyz_data$Pos_Y
  zdata <- xyz_data$Pos_Z

  ### Defining the distance function that has to be minimized
  ### A conic plane of the form z=z0+a*(x-x0)**2+b*(y-y0)**2 with 5 parameters is used.
  ### "k" gives the parameter vector
  dquad <- function(k){
    a <- k[1]
    b <- k[2]
    x0 <- k[3]
    y0 <- k[4]
    z0 <- k[5]
    z <- z0 + a*(xdata-x0)**2 + b*(ydata-y0)**2
    out <- sum((zdata-z)**2)
    return(out)
  }
  localopt <- do.call(rbind,lapply(1:nrfits, function(i){
    ### Initialization of the parameters
    a <- rnorm(1, -3, 1)
    b <- rnorm(1, -3, 1)
    x0 <- rnorm(1, mean(xdata), 1/5*sd(xdata))
    y0 <- rnorm(1, mean(ydata), 1/5*sd(ydata))
    z0 <- rnorm(1, mean(zdata), 1/5*sd(zdata))
    ### Optimization with the function optim
    opt <- optim(par=c(exp(a),exp(b),x0,y0,z0), dquad)
    c(opt$value, opt$par)
  }))

  ### Change data type to data.frame and give names
  localopt <- data.frame(localopt)
  names(localopt) <- c("value", "a", "b", "x0", "y0", "z0")
  ### Order the results with respect to their values
  localopt <- localopt[order(localopt$value),]

  ### Extract best fit
  globalopt <- localopt[1,]
  a <- globalopt$a
  b <- globalopt$b
  x0 <- globalopt$x0
  y0 <- globalopt$y0
  z0 <- globalopt$z0

  ### Calculate corrected z values and add them to myInput
  z_corr <- zdata - (z0 + a*(xdata-x0)**2 + b*(ydata-y0)**2)
  output <- data.frame(xdata, ydata, zdata, z_corr=z_corr)
  Total <- nrow(output)
  ##### Count number of invaded cells
  Invaded <- nrow(output[output$z_corr < threshold,])
  Percentage <- Invaded/Total * 100
  cat("Number of cells invaded: ", as.character(Invaded), "\n")
  cat("Percentage: ", as.character(round(Percentage, 2)))

  ### Exemplary 3D plots before (P1) and after (P2) the correction
  if(plot){
    P1 <- lattice::cloud(Pos_Z ~ Pos_X * Pos_Y, data = output, zlim = c(min(zdata),max(zdata)), screen = list(z = 410, x = -70), panel.aspect = 0.75, xlab = "x", ylab = "y", zlab = "z")
    print(P1)

    P2 <- lattice::cloud(z_corr ~ Pos_X * Pos_Y, data = output, zlim = c(min(z_corr),max(z_corr)), screen = list(z = 410, x = -70), panel.aspect = 0.75, xlab = "x", ylab = "y", zlab = "z")
    print(P2)

  }


  ### Write a new csv file with the corrected values
  if(write_csv)  write.csv(output, file=paste0(strsplit(filename, ".csv")[[1]][1], "_corrected.csv"))


  return(as.data.frame(output))
}
