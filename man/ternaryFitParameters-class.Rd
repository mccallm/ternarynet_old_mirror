\name{ternaryFitParameters-class}
\docType{class}
\alias{class:ternaryFitParameters}
\alias{ternaryFitParameters-class}
\alias{ternaryFitParameters}
\alias{ternaryFitParameters-methods}
\alias{initialize,ternaryFitParameters-method}
\alias{show,ternaryFitParameters-method}
\alias{perturbationType}
\alias{perturbationType,ternaryFitParameters-method}
\alias{scoreType}
\alias{scoreType,ternaryFitParameters-method}
\alias{backupStage}
\alias{backupStage,ternaryFitParameters-method}
\alias{maxStage}
\alias{maxStage,ternaryFitParameters-method}
\alias{maxTransition}
\alias{maxTransition,ternaryFitParameters-method}
\alias{epsilon}
\alias{epsilon,ternaryFitParameters-method}
\alias{beta0}
\alias{beta0,ternaryFitParameters-method}
\alias{chi0}
\alias{chi0,ternaryFitParameters-method}
\alias{delta}
\alias{delta,ternaryFitParameters-method}
\alias{ne}
\alias{ne,ternaryFitParameters-method}
\alias{m0}
\alias{m0,ternaryFitParameters-method}
\alias{maxDegree}
\alias{maxDegree,ternaryFitParameters-method}
\alias{pAddParent}
\alias{pAddParent,ternaryFitParameters-method}
\alias{pExchangeParent}
\alias{pExchangeParent,ternaryFitParameters-method}
\alias{neighborDegree}
\alias{neighborDegree,ternaryFitParameters-method}
\alias{pNeighborhood}
\alias{pNeighborhood,ternaryFitParameters-method}
\alias{rho}
\alias{rho,ternaryFitParameters-method}
\alias{edgePenalty}
\alias{edgePenalty,ternaryFitParameters-method}
\alias{perturbationType<-}
\alias{perturbationType<-,ternaryFitParameters-method}
\alias{scoreType<-}
\alias{scoreType<-,ternaryFitParameters-method}
\alias{backupStage<-}
\alias{backupStage<-,ternaryFitParameters-method}
\alias{maxStage<-}
\alias{maxStage<-,ternaryFitParameters-method}
\alias{maxTransition<-}
\alias{maxTransition<-,ternaryFitParameters-method}
\alias{epsilon<-}
\alias{epsilon<-,ternaryFitParameters-method}
\alias{beta0<-}
\alias{beta0<-,ternaryFitParameters-method}
\alias{chi0<-}
\alias{chi0<-,ternaryFitParameters-method}
\alias{delta<-}
\alias{delta<-,ternaryFitParameters-method}
\alias{ne<-}
\alias{ne<-,ternaryFitParameters-method}
\alias{m0<-}
\alias{m0<-,ternaryFitParameters-method}
\alias{maxDegree<-}
\alias{maxDegree<-,ternaryFitParameters-method}
\alias{pAddParent<-}
\alias{pAddParent<-,ternaryFitParameters-method}
\alias{pExchangeParent<-}
\alias{pExchangeParent<-,ternaryFitParameters-method}
\alias{neighborDegree<-}
\alias{neighborDegree<-,ternaryFitParameters-method}
\alias{pNeighborhood<-}
\alias{pNeighborhood<-,ternaryFitParameters-method}
\alias{rho<-}
\alias{rho<-,ternaryFitParameters-method}
\alias{edgePenalty<-}
\alias{edgePenalty<-,ternaryFitParameters-method}

\title{Ternary Network Fitting Parameters}
\description{This is a class representation of the input parameters for
  the ternary network fitting algorithm implemented in the function
  \code{tnetfit}.}

\section{Creating Objects}{

  \code{
    ternaryFitParameters()
  }  

  This creates a \code{ternaryFitParameters} object with the default
  fitting parameters. 
  }
  
\section{Slots}{
  \describe{
    \item{\code{perturbationType}:}{this parameter currently can only be
  set to 1} 
    \item{\code{scoreType}:}{the method to score networks. Can be set to
  either 1 or 2, corresponding the the score types in Almudevar et
  al. (2011).} 
    \item{\code{backupStage}:}{current fit is output periodically
  according to this parameter}
    \item{\code{maxStage}:}{the maximum number of stages
  permitted. Ideally, the actual number of stages required until
  convergence should be much less than this value.} 
    \item{\code{maxTransition}:}{This parameter provides an adaptive
  truncation of the stage sample size. The stage terminates before the
  specified fixed sample size if the number of transitions resulting in
  a strict increase of the score reaches this value. If the sampler is
  in steady state, then this count should be approximately half the
  number of transitions in which the score changes value.}
    \item{\code{epsilon}:}{Convergence tolerance.}
    \item{\code{beta0}:}{Algorithm terminates when this number of
  consecutive convergence events have occurred.} 
    \item{\code{chi0}:}{The target initial acceptance rate. This should
  be close to 1, although setting it too close will increase computation
  time.}
    \item{\code{delta}:}{The increment change in steady state
  distribution between stages (as variational distance). Larger values
  tend to decrease computation time, but too large a value will result
  in spurious convergence.}
    \item{\code{ne}:}{The fixed sample size (number of MCMC transitions) per
  stage.}
    \item{\code{m0}:}{The sample size (number of transitions) used to
  determine the initial temperature.}
    \item{\code{maxDegree}:}{Maximum number of parents per node
  permitted in model topology.}
    \item{\code{pAddParent}:}{This is the probability of adding a parent
  to a randomly selected node in the proposal function.}
    \item{\code{pExchangeParent}:}{This parameter gives the probability
  of a parent exchange in the proposal function. The AddParent operation
  takes precedence, so this probability should be interpreted as being
  conditional on the rejection of the AddParent operation.}
    \item{\code{neighborDegree}:}{Number of applications of the proposal
  function.}
    \item{\code{pNeighborhood}:}{Vector of probabilities denoted, which
  generates the random number of proposal function iterations. The
  length is one less than neighborDegree. If neighborDegree equals 1
  then no iteration is performed, and this vector is ignored.}
    \item{\code{rho}:}{Weight parameter for the exponential
  smoothing of the variance estimate. For no smoothing set to 1.}
    \item{\code{edgePenalty}:}{This parameter provides a complexity
  penalty. This number times the number of edges is added to the
  score. To apply no penalty set this parameter to 0.}
  }
}

\section{Methods}{
  All named elements can be accessed and set in the standard way
  (e.g. \code{scoreType(object)} and \code{scoreType(object)<-}).
}
\author{Matthew N. McCall and Anthony Almudevar}
\seealso{
  \code{tnetfit}, \code{ternaryFit-class}, \code{ternaryPost-class}.
  Almudevar A, McCall MN, McMurray H, Land H (2011). Fitting
  Boolean Networks from Steady State Perturbation Data, Statistical
  Applications in Genetics and Molecular Biology, 10(1): Article 47.
  }

\examples{
# create an instance
ternaryFitParameters()
}

\keyword{classes}