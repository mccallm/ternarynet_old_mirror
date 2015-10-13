\name{graphPosterior}
\alias{graphPosterior}
\title{Network Topology}
\description{
  This function summarizes the topology of the ternary network using
  marginal edge probabilities.
}
\usage{
graphPosterior(tpost) 
}
\arguments{
  \item{tpost}{a ternaryPost object}
}
\value{
  The function returns a matrix of marginal posterior probabilities of
  each possible network edge -- rows are children and columns are
  parents. The first column represents no parents. 
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
rownames(ssObj) <- rownames(pObj) <- colnames(ssObj) <- colnames(pObj) <- c("Gene1","Gene2","Gene3")
tnfitObj <- tnetfit(ssObj, pObj)
tnpostObj <- tnetpost(tnfitObj, mdelta=10, msample=10)
graphPosterior(tnpostObj)
}
\keyword{manip}

