#' Plots the joint diagonalization. I.e. if it was successful the matrices should all
#' be diagonal.
#' 
#' @param estConnectivity Estimate for connectivity matrix returned by \code{backShift}.
#' @param X Data matrix
#' @param env Indicator of the experiment or the intervention type an 
#' observation belongs to (a numeric vector of length n).
#' @param whichEnv Indicator for the environment to be plotted.
#'  
plotDiagonalization <- function(estConnectivity, X, env, whichEnv){
  deltas <- computeDelta(X, env)$Delta
  p <- ncol(estConnectivity)
  K <- length(deltas)
  
  deltaIntVar <- vector("list", K)
  
  tmp <- (diag(p) - t(estConnectivity))%*%deltas[[whichEnv]]%*%t(diag(p) - t(estConnectivity))
  deltaIntVar <- abs(tmp)/max(abs(tmp))
  
  fields::image.plot(z = deltaIntVar, col=fields::tim.colors(), legend.shrink = 0.8,
                     main = paste("Env.", whichEnv))

}
