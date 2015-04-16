\name{simulateSteadyState}
\alias{simulateSteadyState}
\title{Simulate Steady State Data}
\description{
  This function generates simulated steady state data from a given network.
}
\usage{
simulateSteadyState(perturbationObj, tableObj, graphObj, degreeObj, wildtype=FALSE)
}
\arguments{
  \item{perturbationObj}{a matrix of perturbation experiments. Rows are genes and columns are experiments.}
  \item{tableObj}{a matrix containing the transition function tables}
  \item{graphObj}{a matrix containing the parents of each node}
  \item{degreeObj}{a vector containing the in-degree of each node}
  \item{wildtype}{if TRUE, the preturbations are assumed to be transient; if FALSE, the perturbations are assumed to be persistent.}
}
\value{
  The function creates a steadyStateObj.
}
\author{Matthew N. McCall and Anthony Almudevar}
\seealso{
  Almudevar A, McCall MN, McMurray H, Land H (2011). Fitting Boolean
  Networks from Steady State Perturbation Data, Statistical Applications
  in Genetics and Molecular Biology, 10(1): Article 47.
  }
\examples{
pObj <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)
degreeObj <- c(0,1,1)
graphObj <- matrix(nrow=1,ncol=3)
graphObj[1,1] <- 0
graphObj[1,2] <- 1
graphObj[1,3] <- 2
tableObj <- matrix(nrow=3,ncol=3)
tableObj[,1] <- rep(0,3)
tableObj[,2] <- c(-1,0,1)
tableObj[,3] <- c(-1,0,1)
ssObj <- simulateSteadyState(pObj, tableObj, graphObj, degreeObj)
}
\keyword{manip}

