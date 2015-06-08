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
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%% Used in simulations in 
}
\examples{
data("exampleAdjacencyMatrix")
plotGraphEdgeAttr(estimate = exampleAdjacencyMatrix, plotStabSelec = FALSE, 
                  labels = colnames(exampleAdjacencyMatrix),
                  thres.point = 0, thres.stab = NULL, main = "True graph")
}
\keyword{datasets}
