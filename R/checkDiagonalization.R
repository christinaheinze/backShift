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
  
  image(
    z  = deltaIntVar,
    axes = FALSE, 
    col = fields::tim.colors(),
    main = paste("Env.", whichEnv))
  axis(
    side     = 1,
    labels   = 1:p,
    at       = seq( 0, 1, length = p ),
    cex.axis = 0.85)
  axis(
    side     = 2,
    labels   = 1:p,
    at       = seq( 0, 1, length = p ),
    cex.axis = 0.85)
}

plotDiagonalizationAll <- function(estConnectivity, X, env){
  deltas <- computeDelta(X, env)$Delta
  p <- ncol(estConnectivity)
  K <- length(deltas)
  
  deltaIntVar <- vector("list", K)
  
  for(i in 1:K){
    tmp <- (diag(p) - t(estConnectivity))%*%deltas[[i]]%*%t(diag(p) - t(estConnectivity))
    deltaIntVar[[i]] <- abs(tmp)/max(abs(tmp))
  }
  
  gridSize <- ceiling(sqrt(K))
  
  par(mfrow = c(gridSize, gridSize))
  for(i in 1:K){
    image(
      z  = deltaIntVar[[i]],
      axes = FALSE, 
      col = fields::tim.colors(),
      main = paste("Environment", i))
    axis(
      side     = 1,
      labels   = 1:p,
      at       = seq( 0, 1, length = p ),
      cex.axis = 0.85)
    axis(
      side     = 2,
      labels   = 1:p,
      at       = seq( 0, 1, length = p ),
      cex.axis = 0.85)
    fields::image.plot(1:p, 1:p, deltaIntVar[[i]], legend.only=TRUE, col=fields::tim.colors())
  }
  
}
