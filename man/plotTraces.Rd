\name{plotTraces}
\alias{plotTraces}
\title{Network Fit Traces}
\description{
  This function plots the trace of four model parameters. 
}
\usage{
plotTraces(tfit) 
}
\arguments{
  \item{tfit}{a ternaryFit object}
}
\value{
  The function creates a 2x2 grid of the four trace plots. 
}
\author{Matthew N. McCall and Anthony Almudevar}
\seealso{
  Almudevar A, McCall MN, McMurray H, Land H (2011). Fitting Boolean
  Networks from Steady State Perturbation Data, Statistical Applications
  in Genetics and Molecular Biology, 10(1): Article 47.
  }
\examples{
ssObj <- matrix(c(1,1,1,0,1,1,0,0,1),nrow=3)
pObj <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)
tnfitObj <- tnetfit(ssObj, pObj)
plotTraces(tnfitObj)
}
\keyword{manip}

