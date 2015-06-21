\name{backShift}
\alias{backShift}
\title{
Estimate connectivity matrix of a directed graph with linear effects and 
hidden variables. 
}
\description{
This function estimates the connectivity matrix of a directed
(possibly cyclic) graph with hidden variables. The underlying system is required 
to be linear and we assume that observations under different shift interventions 
are available. More precisely, the function takes as an input an (nxp) data matrix, 
where \eqn{n} is the sample size and \eqn{p} the number of variables. In each 
environment \eqn{j} (\eqn{j} in \{\eqn{1, \ldots, J}\}) we have observed 
\eqn{n_j} samples generated from \deqn{X_j= X_j \%*\% A + c_j + e_j}
(in case of cycles this should be understood as an equilibrium
distribution). The \eqn{c_j} is a p-dimensional random vector that is assumed 
to have a diagonal covariance matrix. The noise vector \eqn{e_j} is assumed to 
have the same distribution in all environments \eqn{j} but is allowed to have 
an arbitrary covariance matrix. The different intervention settings are provided 
to the method with the help of the vector \code{ExpInd} of length 
\eqn{n = (n_1 + ... + n_j + ... + n_J)}. The goal is to estimate the 
connectivity matrix \eqn{A}.

}
\usage{
backShift(X, ExpInd, covariance=TRUE, ev=0, threshold =0.75, nsim=100, 
          sampleSettings=1/sqrt(2), sampleObservations=1/sqrt(2), 
          nodewise=TRUE, tolerance = 10^(-4), baseSettingEnv = 1, 
          verbose = FALSE)
}
\arguments{
  \item{X}{
  A (nxp)-dimensional matrix (or data frame) with n observations of p variables.
}
  \item{ExpInd}{
  Indicator of the experiment or the intervention type an observation belongs to. 
  A numeric vector of length n. Has to contain at least three different unique values. 
}
 \item{covariance}{
A boolean variable. If \code{TRUE}, use only shift in covariance matrix; 
otherwise use shift in Gram matrix.
Set only to \code{FALSE} if at most one variable has a non-zero shift in mean 
in the same setting (default is \code{TRUE}).
}
  \item{ev}{
 The expected number of false selections for stability selection. 
 No stability selection computed if \code{ev=0}. Defaults to \code{ev=0}.
}
  \item{threshold}{
  The selection threshold for stability selection (has to be between 0.5 and 1). 
  Edges which are selected with empirical proportion higher than \code{threshold} 
  will be retained. 
}
  \item{nsim}{
  Number of resamples taken (if using stability selection). 
}
 \item{sampleSettings}{
  The proportion of unique settings to resample for each resample; has to be in [0,1].
}

 \item{sampleObservations}{
  The fraction of all samples to retain when subsampling (no replacement); has to be in [0,1].
}
\item{nodewise}{
   If \code{FALSE}, stability selection retains for each subsample the largest 
   overall entries in the connectivity matrix. If \code{TRUE}, values are 
   ordered row- and node-wise first and then the largest entries in each row 
   and column are retained. Error control is valid (under exchangeability 
   assumption) in both cases. The latter setting \code{TRUE} is perhaps more 
   robust and is the default. 
}
\item{tolerance}{
   Precision parameter for \code{ffdiag}: the algorithm stops when the criterium 
   difference between two iterations is less than \code{tolerance}. Default is 10^(-4).
}
\item{baseSettingEnv}{
Index for baseline environment against which the intervention variances are measured.
  Defaults to 1.
}
\item{verbose}{
  If \code{FALSE}, most messages are supressed.
}
}

\value{
A list with elements
  \item{Ahat}{The connectivity matrix where entry (i,j) is the effect pointing 
  from variable i to variable j.}
  \item{AhatAdjacency}{If \code{ev}>0, the connectivity matrix retained by 
  stability selection. Entries give the rounded percentage of times the edge has 
  been retained (and 0 if below the critical threshold).}
  \item{varianceEnv}{
       A (Gxp)-dimensional matrix where G is the number of unique environments.
       The ij-th entry contains the difference
       between the est. intervention variance of variable j in the base setting and the
       estimated variance of variable j in environment i.
       }
}
\references{
 Dominik Rothenhaeusler, Christina Heinze, Jonas Peters, Nicolai Meinshausen (2015):
 backShift: Learning causal cyclic graphs from unknown shift interventions.
 arXiv preprint: \url{http://arxiv.org/abs/1506.02494}
}
\author{
Christina Heinze <heinze@stat.math.ethz.ch>
}
\seealso{  \code{\link[InvariantCausalPrediction]{ICP}} and 
  \code{\link[InvariantCausalPrediction]{hiddenICP}} for reconstructing 
  the parents of a variable under interventions on all other variables.
}
\examples{
## Simulate data with connectivity matrix A

seed <- 1
# sample size n
n <- 10000
# 3 predictor variables
p  <- 3
A <- diag(p)*0
A[1,2] <- 0.8
A[2,3] <- -0.8
A[3,1] <- 0.8

# divide data into 10 different environments
G <- 10

# simulate
simulation.res <- simulateInterventions(
                    n, p, A, G, intervMultiplier = 2, 
                    noiseMult = 1, nonGauss = FALSE, 
                    fracVarInt = 0.5, hidden = TRUE, 
                    knownInterventions = FALSE, 
                    simulateObs = TRUE, seed)

environment <- simulation.res$environment
X <- simulation.res$X

## Compute feedback estimator with stability selection

network <- backShift(X, environment, ev = 1)

## Print point estimates and stable edges

# true connectivity matrix
print(A) 
# point estimate
print(network$Ahat)       
# shows empirical selection probability for stable edges
print(network$AhatAdjacency)
}

\keyword{Causality }
\keyword{Regression }
\keyword{Hidden Variables }
\keyword{Feedback }  
