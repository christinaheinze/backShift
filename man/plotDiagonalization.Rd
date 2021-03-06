% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkDiagonalization.R
\name{plotDiagonalization}
\alias{plotDiagonalization}
\title{Plots the joint diagonalization. I.e. if it was successful the matrices should all
be diagonal.}
\usage{
plotDiagonalization(estConnectivity, X, env, whichEnv, main = NULL)
}
\arguments{
\item{estConnectivity}{Estimate for connectivity matrix returned by \code{backShift}.}

\item{X}{Data matrix}

\item{env}{Indicator of the experiment or the intervention type an 
observation belongs to (a numeric vector of length n).}

\item{whichEnv}{Indicator for the environment to be plotted.}

\item{main}{Optional title for plot; defaults to paste("Env.", whichEnv)}
}
\description{
Plots the joint diagonalization. I.e. if it was successful the matrices should all
be diagonal.
}
