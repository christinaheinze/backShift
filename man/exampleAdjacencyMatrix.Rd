\name{exampleAdjacencyMatrix}
\alias{exampleAdjacencyMatrix}
\docType{data}
\title{
Example adjacency matrix
}
\description{
An example for an adjacency matrix A to be used as input to 
\code{\link{simulateInterventions}}. The entry \eqn{A_{ij}} contains the edge
from node i to node j.
}
\usage{data("exampleAdjacencyMatrix")}
\format{
  A matrix with 10 rows and 10 columns. 
}
\references{
  Used in simulations in:

 Dominik Rothenhaeusler, Christina Heinze, Jonas Peters, Nicolai Meinshausen (2015):
 backShift: Learning causal cyclic graphs from unknown shift interventions
 %arxiv preprint \url{http://arxiv.org/abs/}
}
\examples{
data("exampleAdjacencyMatrix")
plotGraphEdgeAttr(estimate = exampleAdjacencyMatrix, plotStabSelec = FALSE, 
                  labels = colnames(exampleAdjacencyMatrix),
                  thres.point = 0, thres.stab = NULL, main = "True graph")
}
\keyword{datasets}
