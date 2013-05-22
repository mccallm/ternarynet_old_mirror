\name{tnetpost}
\alias{tnetpost}
\title{Ternary Network Posterior Sampling}
\description{
  This function samples from the posterior density of a ternary network based on perturbation experiments.
}
\usage{
tnetpost(tfit, mdelta=as.integer(10000), msample=as.integer(2000), temperatureScale=1.0, xSeed=NA)
}
\arguments{
  \item{tfit}{a ternaryFit object}
  \item{mdelta}{number of transitions between samples}
  \item{msample}{number of samples}
  \item{temperatureScale}{the final temperature is multipled by this
    value for sampling}
  \item{xSeed}{an integer random seed. If NA, a random seed is generated.}
}
\value{
  The function returns a ternaryPost object.
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
}
\keyword{manip}

