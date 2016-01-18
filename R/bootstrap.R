
# bootstrappedDesignModel <- function(Ahat, X, env, seed, verbose = FALSE){
#   set.seed(seed)
#   uniqueEnvs <- unique(env)
#   n <- nrow(X)
#   p <- ncol(X)
#   Xboot <- matrix(0, nrow = n, ncol = p)
# 
#   for(i in seq_along(uniqueEnvs)){
#     if(verbose) cat(paste("\nResampling residuals for env.", uniqueEnvs[i]))
#     Xe <- X[env == uniqueEnvs[i],]
#     ne <- nrow(Xe)
#     Re <- Xe - Xe%*%Ahat
#     bootSampleResid <- matrix(sample(Re, size = ne*p, replace = TRUE), nrow = ne)
#     Xboot[env == uniqueEnvs[i],] <- bootSampleResid%*%solve(diag(p) - Ahat)
#   }
# 
#   Xboot
# }

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

# bootstrappedDesignParametric <- function(Ahat, X, env, seed, nrep, verbose = FALSE){
#   set.seed(seed)
#   uniqueEnvs <- unique(env)
#   estimates <- parametricNoisePars(Ahat, X, env, verbose = verbose)
# 
#   n <- nrow(X)
#   p <- ncol(X)
#   XbootReps <- vector("list", nrep)
#   
#   for(j in 1:nrep){
#     if(verbose) if(j %% 10 == 0) cat(paste("\nGenerating design matrix for rep.", j))
#     Xboot <- matrix(0, nrow = n, ncol = p)
#     for(i in seq_along(estimates)){
#       ne <- sum(env == uniqueEnvs[i])
#       est <- estimates[[i]]
#       errorsAndInt <- mvrnorm(n = ne, mu = est$muhat, Sigma = est$sigmahat)
#       Xboot[env == uniqueEnvs[i],] <- errorsAndInt%*%solve(diag(p) - Ahat)
#     }
#     XbootReps[[j]] <- Xboot
#   }
#   XbootReps
# }
# 
# 
# bootstrapSumOffDiagnols <- function(Ahat, X, env, nrep, seed = 1,
#                                covariance=TRUE, baseInd = 1,
#                                tolerance = 1e-3, verbose = FALSE){
#   bootstrapReps <- bootstrappedDesignParametric(Ahat, X, env, seed, nrep, verbose)
#   k <- 1
#   sumOffDiagnols <- sapply(bootstrapReps, function(Xi){
#                               if(verbose){
#                                 if(k %% 10 == 0) cat(paste("\nbackShift estimation rep.", k))
#                                 k <<- k + 1
#                               }
#                               Ahati <- backShift(Xi, env, covariance = covariance, ev=0,
#                                         baseSettingEnv = baseInd, tolerance = tolerance,
#                                         verbose = FALSE)$Ahat
#                               jointDiagSumOfOffDiagnalElements(Ahati, Xi, env)
#                             }
#                            )
#   sumOffDiagnols
# }


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
    if(verbose) if(j %% 10 == 0) cat(paste("\nGenerating design matrix and running backShift for rep.", j))
    
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

bootstrapBackShift <- function(Ahat, X, env, nrep, alpha = 0.05,
                                    covariance=TRUE, baseInd = 1,
                                    tolerance = 1e-3, verbose = FALSE){

  bootsSumOffDiags <- bootstrapSumOffDiagnols(Ahat, X, env, nrep,
                                                     covariance=covariance,
                                                     baseInd = baseInd,
                                                     tolerance = tolerance,
                                                     verbose = verbose)
  sumOffDiagsBackShift <- jointDiagSumOfOffDiagnalElements(Ahat, X, env)

  quantile1 <- quantile(bootsSumOffDiags, alpha/2)
  quantile2 <- quantile(bootsSumOffDiags, 1-alpha/2)
  lowerBasic <- 2*sumOffDiagsBackShift - quantile2
  upperBasic <- 2*sumOffDiagsBackShift - quantile1
  tmp <- names(upperBasic)
  names(upperBasic) <- names(lowerBasic)
  names(lowerBasic) <- tmp
  
  decision <- sumOffDiagsBackShift <= quantile2 & sumOffDiagsBackShift >= quantile1
  decisionBasic <- sumOffDiagsBackShift <= upperBasic & sumOffDiagsBackShift >= lowerBasic
  
  names(decision) <- NULL
  
  list(bootsSumOffDiags = bootsSumOffDiags, 
       sumOffDiagsBackShift = sumOffDiagsBackShift, 
       jointDiagSuccess = decision, 
       jointDiagSuccessBasic = decisionBasic, 
       lower = quantile1, 
       upper = quantile2, 
       lowerBasic = lowerBasic, 
       upperBasic = upperBasic)
}
