\name{hiddenICE}
\alias{hiddenICE}
\title{
Estimate connectivity matrix of a directed graph with linear effects and hidden variables. 
}
\description{
Estimates the connectivity matrix of a directed graph where effects are assumed to be linear, hidden variables are allowed and interventions have to be additive in nature. 
It is assumed that equilibrium data are observed so that X = X \%*\% A + I, where X is the matrix of observations, A the connectivity matrix and I is the noise for each variable.
(The noise I can be correlated across variables, but has to change between the different settings.)
The goal is to estimate the connectivity matrix A.
}
\usage{
hiddenICE(X, ExpInd, covariance=TRUE, alpha = 0.1, 
                           threshold =0.75, nsim=100, 
                           sampleSettings=1/sqrt(2), sampleObservations=1/sqrt(2), 
                           nodewise=TRUE, warnings= FALSE )
}
\arguments{
  \item{X}{
  A (nxp)-dimensional matrix (or data frame) with n observations of p variables
}
  \item{ExpInd}{
  Indicator of the experiment or the intervention type an observation belongs to. 
  A numeric vector of length n. Has to contain at least three different unique values. 
}
 \item{covariance}{
A boolean variable. If \code{TRUE}, use only shift in covariance matrix; otherwise use shift in Gram matrix.
Set only to \code{FALSE} if at most one variable has a non-zero shift in mean in the same setting (default is \code{TRUE}).
}
  \item{alpha}{
 The expected number of false selections for stability selection. No stability selection computed if \code{alpha}=0.
}
  \item{threshold}{
  The selection threshold for stability selection (has to be between 0.5 and 1). Edges which are selected with empirical proportion higher than \code{threshold} will be retained. 
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
   If \code{FALSE}, stability selection retains for each subsample the largest overall entries in the connectivity matrix.
 If \code{FALSE}, values are ordered row- and node-wise first and then the largest entries in each row and column are retained. 
Error control is valid (under exchangeability assumption) in both cases. The latter setting \code{FALSE} is perhaps more robust and is the default. 
}
\item{warnings}{
   If \code{FALSE}, warning messages are supressed.
}
\item{tolerance}{
   Precision parameter for \code{ffdiag}: the algorithm stops when the criterium difference between two iterations is less than \code{tolerance}. Default is 10^(-4).
}
}

\value{
A list with elements
  \item{Ahat}{The connectivity matrix where entry (i,j) is the effect pointing from variable i to variable j.}
  \item{AhatAdjacency}{If \code{alpha}>0, the connectivity matrix retained by stability selection. 
       Entries give the rounded percentage of times the edge has been retained (and 0 if below the critical threshold).}
}
\references{
    None yet
}
\author{
Anonymous
}

\examples{
 ##########################################
 ####### 1st example:
 ####### Simulate data with connectivity matrix A
      set.seed(1)
    ## sample size n
      n <- 1000
    ## 3 predictor variables
      p  <- 3
      A <- diag(p)*0
      A[1,2] <- 0.8
      A[2,3] <- -0.8
      A[3,1] <- 0.8
    
    ## divide data into 10 different environments
      G <- 10
      ExpInd <- rep(1:G, each=ceiling(n/G))[1:n]
      X <- Perturb <-  matrix(0,nrow=n,ncol=p)
    ## Input of hidden variables into each variable
      alpha <- rnorm(p)
      Input <- outer(W <- rnorm(n),alpha,FUN="*")
    ## change level of noise in each setting 
      for (i in unique(ExpInd)){
        ind <- which(ExpInd==i)
        multiplier <- rexp(p)
        Perturb[ind,] <- sweep( matrix(rnorm(length(ind)*p),ncol=p), 2, multiplier,FUN="*")
     }
  
    ## iterate feedback model for get stable solution
     niter <- 100
     for (iter in 1:niter){
         X <- X \%*\% A + Input + Perturb 
     }
   

####### Compute feedback estimator with E(V)<=1
      network <- hiddenICE(X,ExpInd,alpha=1)

 ###### Print point estimates and stable edges
      print(A) ## true connectivity matrix
      print(network$Ahat)          ## point estimate
      print(network$AhatAdjacency) ## shows empirical selection probability for stable edges

}

\keyword{Causality }
\keyword{Regression }
\keyword{Hidden Variables }
\keyword{Feedback }  
