#' Computes a simple model-based bootstrap confidence interval for success of joint diagonalization
#' procedure. The model-based bootstrap approach assumes normally distributed error terms; the
#' parameters of the noise distribution are estimated with maximum likelihood.  
#'
#' @param Ahat Estimated connectivity matrix returned by \code{backShift}.
#' @param X A (nxp)-dimensional matrix (or data frame) with n observations of p variables.
#' @param ExpInd Indicator of the experiment or the intervention type an observation belongs to. A numeric vector of length n. Has to contain at least three different unique values.
#' @param nrep Number of bootstrap samples.
#' @param alpha Significance level for confidence interval.
#' @param covariance A boolean variable. If \code{TRUE}, use only shift in covariance matrix; otherwise use shift in Gram matrix. Set only to \code{FALSE} if at most one variable has a non-zero shift in mean in the same setting (default is \code{TRUE}).
#' @param baseInd Index for baseline environment against which the intervention variances are measured. Defaults to 1.
#' @param tolerance  Precision parameter for \code{ffdiag}: the algorithm stops when the criterium difference between two iterations is less than \code{tolerance}. Default is 10^(-4).
#' @param verbose  If \code{FALSE}, messages are supressed.
#' 
#' @return A list with the following elements: 
#' \itemize{
#'  \item \code{bootsSumOffDiags} Vector of length \code{nrep} with sum of off-diagonal elements after joint diagnolization procedure in each of the bootstrap samples.
#'  \item \code{sumOffDiagsBackShift} Sum of off-diagonal elements after joint diagnolization procedure in original estimation.
#'  \item \code{jointDiagSuccess} \code{TRUE} if \code{sumOffDiagsBackShift} lies 
#'  within bootstrap confidence interval.
#'  \item \code{lower} Lower bound of bootstrap confidence interval. 
#'  \item \code{upper} Upper bound of bootstrap confidence interval. 
#'  \item \code{lowerBasic} \code{alpha/2} quantile of empirical bootstrap distribution.
#'  \item \code{upperBasic} \code{1 - alpha/2} quantile of empirical bootstrap distribution.
#'  }
#'  
bootstrapBackShift <- function(Ahat, X, ExpInd, nrep, alpha = 0.05,
                               covariance=TRUE, baseInd = 1,
                               tolerance = 1e-3, verbose = FALSE){
  
  bootsSumOffDiags <- bootstrapSumOffDiagnols(Ahat, X, ExpInd, nrep,
                                              covariance=covariance,
                                              baseInd = baseInd,
                                              tolerance = tolerance,
                                              verbose = verbose)
  sumOffDiagsBackShift <- jointDiagSumOfOffDiagnalElements(Ahat, X, ExpInd)
  
  quantile1 <- quantile(bootsSumOffDiags, alpha/2)
  quantile2 <- quantile(bootsSumOffDiags, 1-alpha/2)
  
  meanBoots <- mean(bootsSumOffDiags)
  sdBoots <- sd(bootsSumOffDiags)
  lower <- meanBoots - qnorm(1-alpha/2)*sdBoots
  upper <- meanBoots + qnorm(1-alpha/2)*sdBoots
 
  decision <- (sumOffDiagsBackShift <= upper) & (sumOffDiagsBackShift >= lower)

  names(decision) <- NULL
  
  list(bootsSumOffDiags = bootsSumOffDiags, 
       sumOffDiagsBackShift = sumOffDiagsBackShift, 
       jointDiagSuccess = decision, 
       lower = lower, 
       upper = upper, 
       quantileLower = quantile1, 
       quantileUpper = quantile2)
}

parametricNoisePars <- function(Ahat, X, env, verbose = FALSE){
  uniqueEnvs <- unique(env)
  nEnvs <- length(uniqueEnvs)  
  estimates <- vector("list", nEnvs)
  
  for(i in seq_along(uniqueEnvs)){
    if(verbose) cat(paste("\nEstimating noise parameters for env.", uniqueEnvs[i]))
    Xe <- X[env == uniqueEnvs[i],]
    ne <- nrow(Xe)
    Re <- Xe - Xe%*%Ahat
    estimates[[i]] <- mlest(Re)
  }
  estimates
}

bootstrapSumOffDiagnols <- function(Ahat, X, env, nrep, seed = 1,
                                    covariance=TRUE, baseInd = 1,
                                    tolerance = 1e-3, verbose = FALSE){
  set.seed(seed)
  uniqueEnvs <- unique(env)
  estimates <- parametricNoisePars(Ahat, X, env, verbose = verbose)
  
  n <- nrow(X)
  p <- ncol(X)
  sumOffDiagnols <- numeric(nrep)
  
  for(j in 1:nrep){
    if(verbose){
      if(j %% 10 == 0) 
        cat(paste("\nGenerating design matrix and running backShift for rep.", j))
    }
    
    Xboot <- matrix(0, nrow = n, ncol = p)
    
    for(i in seq_along(estimates)){
      ne <- sum(env == uniqueEnvs[i])
      est <- estimates[[i]]
      errorsAndInt <- mvrnorm(n = ne, mu = est$muhat, Sigma = est$sigmahat)
      Xboot[env == uniqueEnvs[i],] <- errorsAndInt%*%solve(diag(p) - Ahat)
    }
    
    Ahatj <- backShift(Xboot, env, covariance = covariance, ev=0,
                       baseSettingEnv = baseInd, tolerance = tolerance,
                       verbose = FALSE)$Ahat
    sumOffDiagnols[j] <- jointDiagSumOfOffDiagnalElements(Ahatj, Xboot, env)
  }
  
  sumOffDiagnols
}

jointDiagSumOfOffDiagnalElements <- function(estConnectivity, X, env){
  deltas <- computeDelta(X, env)$Delta
  p <- ncol(estConnectivity)
  K <- length(deltas)
  sumOffDiag <- 0

  for(i in 1:K){
    diagMat <- (diag(p) - t(estConnectivity))%*%deltas[[i]]%*%t(diag(p) - t(estConnectivity))
    sumOffDiag <- sumOffDiag + sum((diagMat - diag(diagMat))^2)
  }

  sumOffDiag
}


