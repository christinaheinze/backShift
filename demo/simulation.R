# Example --------

if(!requireNamespace("pcalg", quietly = TRUE))
  stop("The package 'pcalg' is needed for the examples to 
 work. Please install it.", call. = FALSE)

seed <- 1
set.seed(seed)

## Options for method -------

# if stability selection should not be used, set EV = 0
# number of false selections for stability selection
EV <- 2
# selection threshold for stability selection
thres <- 0.75 
# use covariance matrix instead of Gram matrix
useCov <- TRUE 
# threshold for point estimate 
thres.pe <- 0.25

## Parameters for simulation ------

# number of observations
n <- 10000
# number of variables
p <- 10
# number of environments
G <- 10 
# also simulate observational data?
simulateObs <- TRUE 
# should hidden vars be included?
hidden <- FALSE 
# should the location of the interventions be known?
knownInterventions <- FALSE 
# if the location of the interventions is known, on how many vars. should
# be intervention in each environment (as a fraction of p)
fracVarInt <- 0.5
# multiplier for interventions (m_I in manuscript)
intMult <- 1.5
# multiplier for interventions (m_e in manuscript)
noiseMult <- 1
# simulate non-Gaussian noise? 
nonGauss <- FALSE

# set whether to provide the adjacency matrix A 
# (e.g. by loading the example from the
# paper or by creating an adjacency matrix yourself)
# alternatively, A is generated as well
# NOTE: A = t(B), so the entry A_ij contains the edge from 
# node i to node j
providedA <- TRUE

if(providedA){
  data("exampleAdjacencyMatrix")
  A <- exampleAdjacencyMatrix
  p <- 10
}else{
  # parameters to generate A
  
  # should A be cyclic?
  cyclic <- TRUE
  # expected number of neighbors per node
  expNumNeigh <-0.1*p 
  # range for coefficients
  minCoef <- 0.3
  maxCoef <- 0.8
  
  ## Generate A -------
  cat("Generating A...\n") 
  A.gen.result <- generateA(p, expNumNeigh, minCoef, maxCoef, cyclic)
  A <- A.gen.result$A
  cat("A has a cycle of size", A.gen.result$sizeCycle, "\n") 
}

## Simulate data -------
cat("Simulating the data...\n") 
simulation.res <- simulateInterventions(n, p, A, G, intMult, noiseMult, 
                                        nonGauss, hidden, knownInterventions, 
                                        fracVarInt, simulateObs, seed)


# extract X, environment vector and index of observational data
X <- simulation.res$X
env <- simulation.res$environment
baseInd <- simulation.res$configs$indexObservationalData

## Run backShift -------
cat("Running backShift...\n") 
backshift.res <- backShift(X, env, covariance=useCov, ev=EV, threshold=thres, 
                        baseSettingEnv = baseInd, tolerance = 1e-6, 
                        verbose = FALSE)
cat("Running backShift...done!\n") 

## Results of backShift -------

# extract estimates
Ahat <- backshift.res$Ahat
Ahat.structure <- backshift.res$AhatAdjacency

# compute performance metrics and plot result
cat("Plotting true graph... \n") 
plotGraphEdgeAttr(estimate = A, plotStabSelec = FALSE, labels = colnames(A), 
                  thres.point = 0, thres.stab = thres, main = "True graph")
cat("Plotting point estimate, thresholded at", thres.pe,"... \n") 
plotGraphEdgeAttr(estimate = Ahat, plotStabSelec = FALSE, labels = colnames(A), 
                  thres.point = thres.pe, thres.stab = thres, 
                  main = "Point estimate")
cat("Plotting stability selection result... \n") 
plotGraphEdgeAttr(estimate = Ahat.structure, plotStabSelec = TRUE, 
                  labels = colnames(A), thres.point = thres.pe, 
                  edgeWeights = Ahat, thres.stab = thres, 
                  main = "Stability selection result")

cat("\nMetrics for point estimate, thresholded at", thres.pe,"... \n") 
print(metricsThreshold(A, Ahat, thres = thres.pe))

cat("\nMetrics for stability selection result... \n") 
print(metricsThreshold(A, Ahat.structure, thres = 0))
