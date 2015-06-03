backShift <- function(X, ExpInd, covariance=TRUE, ev=0, threshold =0.75, nsim=100, 
                      sampleSettings=1/sqrt(2), sampleObservations=1/sqrt(2), nodewise=TRUE, 
                      tolerance = 10^(-4), baseSettingEnv = 1, verbose = FALSE){

    if(!is.matrix(X) & !is.data.frame(X)) stop("'X' must be a matrix of data frame")

    X <- as.matrix(X)
    p <- ncol(X)
    
    if(!is.vector(ExpInd)) stop("ExpInd must be a vector")
    if( length(ExpInd) != nrow(X)) stop(" 'ExpInd' must have as many entries as X has rows")
    if(ev<0) stop("ev must be non-negative")
    if( threshold <=0.5 | threshold >1) stop("threshold must be between 0.5 and 1")
    if(nsim < 2) stop("'nsim' must be at least 2 (but usually at least 20-100")
    if(length((settings <- unique(ExpInd)))<2) stop("need at least three different settings")
     if(sampleSettings<=0) stop("sampleSettings needs to be positive")
    if(sampleObservations<=0) stop("sampleObservations needs to be positive")
    if(sampleSettings>1) stop("sampleSettings needs to be at most 1")
    if(sampleObservations>1) stop("sampleObservations needs to be at most 1")
    q <- if(!nodewise) sqrt(ev*(2*threshold-1)*(p^2-p)) else sqrt(ev*(2*threshold-1))
    subs <- sampleSettings* length(settings)
    drawE <- function(x){
        z <- floor(x)
        if( runif(1) <  x-z) z <- z+1
        z <- max(3,z)
        return(z)
    }        
    
    ## point estimator
    
    Deltalist <- computeDelta(X,  ExpInd, covariance=covariance)$Delta
    Delta <- array(unlist(Deltalist), dim = c(nrow(Deltalist[[1]]), ncol(Deltalist[[1]]), length(Deltalist)))
    tryCatch({
      estimatedB <- ffdiag(Delta,eps=tolerance, itermax = 500)$B  
    }, 
    warning=function(w){
      warning("backShift -- point estimate -- diagonalization did not succeed -- result not trustworthy")
      cat("backShift -- point estimate -- diagonalization did not succeed -- result not trustworthy\n")        
    },
    error=function(e){
      warning("backShift -- point estimate -- diagonalization did not succeed. Possible model mispecification,
          returning the empty graph.")
      cat("backShift -- point estimate -- diagonalization did not succeed. Possible model mispecification,
          returning the empty graph.")
      return(list(Ahat=0*diag(p), 
                  AhatAdjacency = 0*diag(p), 
                  varianceEnv = matrix(0, nrow = length(unique(ExpInd)), ncol = p)))
    },
    finally = {
        estimatedB <- try(ffdiag(Delta,eps=tolerance, itermax = 500)$B, silent = TRUE)
      }
    )
    res.point.est <- permuteAndScale(estimatedB)
    Ahat <- res.point.est$Ahat
    varEnv <- computeVarEnv(res.point.est$rescaledDhat, Deltalist, baseSettingEnv)

    if(ev>0){
        if(verbose){
          cat("Starting stability selection... \n")
        }
        
      
        AhatList <- vector("list", nsim)
        numberOfRunsNotConverged <- 0
        
        for(i in 1:nsim){
          if(verbose){
            cat("Stability selection: Iteration", i, "... \n")
          }
            useSettings <- sample( settings, drawE(subs))
            
            useSamples <- NULL
            for(s in 1:length(useSettings)){
              ind <- which(  ExpInd %in% useSettings[s])
              useSamples <- c(useSamples, sort(sample(ind, round(length(ind)*sampleObservations))))
            }
        
            Xcurrent <- X[useSamples,]
            ExpIndcurrent <- ExpInd[useSamples]
            
            ## difference between the covariance matrices of X under the different interventions
            Deltalist <- computeDelta(Xcurrent,  ExpIndcurrent, covariance=covariance)$Delta
            Delta <- array(unlist(Deltalist), dim = c(nrow(Deltalist[[1]]), ncol(Deltalist[[1]]), length(Deltalist)))
                      
            rm(estimatedB)
            tryCatch({
              estimatedB <- ffdiag(Delta,eps=tolerance, itermax = 500)$B
            }, 
              warning=function(w){
                numberOfRunsNotConverged <<- numberOfRunsNotConverged + 1
                # warning("backShift -- stability selection -- diagonalization did not succeed")
              },
              error=function(e){
                cat("backShift -- stability selection -- diagonalization did not succeed.
                     Possible model mispecification or sample size too small.
                     Returning the point estimates only. The adjacency matrix is empty.")
                warning("backShift -- stability selection -- diagonalization did not succeed.
                     Possible model mispecification or sample size too small.
                     Returning the point estimates only. The adjacency matrix is empty.")
                return(list(Ahat=Ahat, AhatAdjacency = 0*diag(p), varianceEnv = varEnv))
              },
              finally = {
                estimatedB <- suppressWarnings(ffdiag(Delta,eps=tolerance, itermax = 500)$B)
              }
            )       
                    
            Ahat.stab.sel <- permuteAndScale(estimatedB)$Ahat            
            AhatList[[i]] <- edgeSelection(Ahat.stab.sel, q,  nodewise=nodewise)
        }

        if(verbose){
          cat("Stability selection...done! \n")
        }
          message("backShift: Percentage of runs in stability selection that converged: ", 100*(nsim-numberOfRunsNotConverged)/nsim, "%")
        
        AhatAdjacency <- edgeRetention(AhatList, threshold, p)
    }else{
        AhatAdjacency <- NULL
    }
    list(Ahat=Ahat, AhatAdjacency = AhatAdjacency, varianceEnv = varEnv)
}
