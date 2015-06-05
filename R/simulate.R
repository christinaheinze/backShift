#' Simulates data.
#' 
#' @param n Number of observations.
#' @param p Number of variables.
#' @param A Connectivity matrix A. The entry \eqn{A_{ij}} contains the edge 
#'        from node i to node j.
#' @param G Number of environments.
#' @param intervMultiplier Regulates the strength of the interventions.
#' @param noiseMult Regulates the noise variance.
#' @param nonGauss Set to \code{TRUE} to generate non-Gaussian noise.
#' @param fracVarInt
#' @param hiddenVars
#' @param knownInterventions
#' @param 
#' @return A list with the following elements: 
#' \itemize{
#'   \item First item
#'   \item Second item
#' }
#' @examples
#' add(1, 1)
#' add(10, 1)
simulateInterventions <- function(n, p, A, G, intervMultiplier, noiseMult, 
                                  nonGauss, fracVarInt, hiddenVars, 
                                  knownInterventions, simulateObs,  seed =1){
  ####### set seed
  set.seed(seed)
  
  ####### initialize
  X <- Perturb <- matrix(0,nrow=n,ncol=p)
  envVar <- matrix(0, nrow = G, ncol = p)  
  if(simulateObs) idxObs <- 1 else idxObs <- NULL
  
  ###### simulate
  if(!knownInterventions){
    ## intervention targets are unknown
    environment <- rep(1:G, each=ceiling(n/G))[1:n]
    
    # additive soft interventions
    
    ## simulate noise perturbations in each environment
    for (i in unique(environment)){
      ind <- which(environment==i)
      if(!(simulateObs & i == 1)){
        # if observational data should be simulated: no interventions in environment 1
        multiplier <- rexp(p)*intervMultiplier
        envVar[i,] <- multiplier
        Perturb[ind,] <- sweep(matrix(rnorm(length(ind)*p),ncol=p), 2, multiplier,FUN="*")
      }
    }
  
    interventions <- NULL
  }else{
    
    ## choose explicity intervention targets
    unique.interventions <- list()
    for (g in 1:G){
      unique.interventions[[g]] <- sample(1:p,round(p*fracVarInt))
      envVar[g, unique.interventions[[g]]] <- 
        rexp(length(unique.interventions[[g]]))*intervMultiplier
    }
    
    interventions <- list()
    for (i in 1:n) interventions[[i]] <- sample(unique.interventions,1)[[1]]
    environment <- match(interventions, unique(interventions))
    
    # additive soft interventions
    ## change level of noise for each intervention
    for (i in 1:n){
      if(!(simulateObs & environment[i] == 1)){
        # if observational data should be simulated: no interventions in environment 1
        Perturb[i, interventions[[i]]] <- 
          rnorm(length(interventions[[i]]))*envVar[environment[i],interventions[[i]]]
      }
    }
    
  }

  if(hiddenVars){
    ## Input of hidden variables into each variable
    gamma <- rnorm(p)*noiseMult
    W <- if(nonGauss){
      expo <- matrix(rexp(n),nrow=n)
      rs <- matrix(1, nrow = n)
      rs[sample(n, size = 0.5*n)] <- -1
      rs * expo
    }else{
      rnorm(n)
    }
    Input <- matrix(outer(W, gamma, FUN="*"), nrow = n)
  }else{
    ## Independent noise at each variable 
    Input <- if(nonGauss){
        expo <- matrix(rexp(n*p),nrow=n)
        rs <- matrix(1, nrow = n, ncol = p)
        rs[sample(n*p, size = 0.5*n*p)] <- -1
        rs * expo
      }else{
        matrix(rnorm(n*p),nrow=n)
      }
  }  
  
  inv <- solve(diag(p) - A)
  X <- (Input + Perturb)%*%inv
  
  list(X = X, trueA = A, environment = environment, 
       interventions = interventions, 
       G = G, n = n, p = p, indexObservationalData = idxObs, varInts = envVar, 
       intervMultiplier = intervMultiplier, noiseMult = noiseMult, 
       fracVarInt = fracVarInt, hiddenVars = hiddenVars, 
       knownInterventions = knownInterventions, simulateObs = simulateObs)
}


generateA <- function(p, expNumNeigh, minCoef, maxCoef, cyclic, verbose = FALSE){
  graph.obj <- randDAG(p, expNumNeigh, wFUN=list(runif, min=minCoef, max=maxCoef))
  A <- as(graph.obj, "matrix")
  
  while(sum(A) == 0){ # do not want to get empty graph
    graph.obj <- randDAG(p, expNumNeigh, wFUN=list(runif, min=minCoef, max=maxCoef))
    A <- as(graph.obj, "matrix")
  }
  
  # reverse sign of half the entries
  nz <- which(abs(A) > 0, arr.ind = F)
  reverse.ind <- sample(1:length(nz), ceiling(length(nz)/2))
  for(i in 1:length(reverse.ind)) A[nz[reverse.ind[i]]] <- - A[nz[reverse.ind[i]]]
  
  if(cyclic){
    candidate <- addCycle(A, p, minCoef, maxCoef)
    
    # check whether I-A is invertible and whether A has CP < 1
    Dhat <- t(diag(p) - candidate$A)
    while(is.singular.matrix(Dhat) && !hasCPsmallerOne(Dhat, FALSE)$success){
      candidate <- addCycle(A, p, minCoef, maxCoef)
      if(verbose) cat("I - A is singular or has CP >= 1, regenerating A...")
    }
    
    A <- candidate$A
    sizeCycle <- candidate$sizeCycle
    if(verbose) cat('Graph contains cycle of size', sizeCycle, '\n')
  }else{
    sizeCycle <- NULL
  }
  list(A = A, sizeCycle = sizeCycle)
}


hasCycles <- function(A, verbose = FALSE){ 
  G <- graph.adjacency(A,mode="directed",weighted="a")
  cyc <- !is.dag(G)
  if(verbose) if(cyc) cat('Graph is cyclic\n') else cat('Graph is acyclic\n')
  cyc
}


addCycle <- function(A, p, minCoef, maxCoef){
  # add cycles
  sample.func <- function(i) 
    sample(c(runif(i, min = minCoef, max = maxCoef), 
             runif(i, min = -maxCoef, max = -minCoef)), i)
  
  while(!hasCycles(A)){
    G <- graph.adjacency(A,mode="directed",weighted="a")
    nz <- which(abs(A) > 0, arr.ind = TRUE)
    sample.node <- sample(as.numeric(nz[,1]), 1)
    dfs.result <- graph.dfs(G, root = sample.node, unreachable = FALSE, dist = TRUE) 
    furthest.dis.from.sample.node <- max(dfs.result$dist, na.rm = TRUE) + 1
    furthest.from.sample.node <- which.max(dfs.result$dist)
    A[furthest.from.sample.node, sample.node] <- sample.func(1)
    diag(A) <- 0
  }
  list(A = A, sizeCycle = furthest.dis.from.sample.node)
}

