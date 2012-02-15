\name{ternaryPost-class}
\docType{class}
\alias{class:ternaryPost}
\alias{ternaryPost-class}
\alias{ternaryPost}
\alias{ternaryPost-methods}
\alias{initialize,ternaryPost-method}
\alias{show,ternaryPost-method}
\alias{dim,ternaryPost-method}
\alias{perturbationObj,ternaryPost-method}
\alias{steadyStateObj,ternaryPost-method}
\alias{scores}
\alias{scores,ternaryPost-method}
\alias{degreeObjs}
\alias{degreeObjs,ternaryPost-method}
\alias{graphObjs}
\alias{graphObjs,ternaryPost-method}
\alias{tableObjs}
\alias{tableObjs,ternaryPost-method}
\alias{inputParams,ternaryPost-method}
\alias{perturbationObj<-,ternaryPost-method}
\alias{steadyStateObj<-,ternaryPost-method}
\alias{scores<-}
\alias{scores<-,ternaryPost-method}
\alias{degreeObjs<-}
\alias{degreeObjs<-,ternaryPost-method}
\alias{graphObjs<-}
\alias{graphObjs<-,ternaryPost-method}
\alias{tableObjs<-}
\alias{tableObjs<-,ternaryPost-method}
\alias{inputParams<-,ternaryPost-method}

\title{Ternary Network Posterior}
\description{This is a class representation of the output of the ternary
  network posterior sampling algorithm implemented in the function
  \code{tnetpost}.} 

\section{Creating Objects}{
  While one can create their own objects using the function
  \code{ternaryPost()}, this is highly discouraged. Typically this class
  is created by running the \code{tnetpost} function.
}
  
\section{Slots}{
  \describe{
    \item{\code{perturbationObj}:}{a matrix of perturbation
  experiments. Rows are genes and columns are experiments.}
    \item{\code{steadyStateObj}:}{a matrix of steady gene expression
  observations from a perturbation experiment. Rows are genes and
  columns are experiments.}
    \item{\code{scores}:}{the score of each sample}
    \item{\code{degreeObjs}:}{the in-degree vector for each sample}
    \item{\code{graphObjs}:}{the graph matrix for each sample}
    \item{\code{tableObjs}:}{the table matrix for each sample}
    \item{\code{inputParams}:}{the ternaryFitParameters object used}
  }
}

\section{Methods}{
  All named elements can be accessed and set in the standard way
  (e.g. \code{scores(object)} and \code{scores(object)<-}).
}

\author{Matthew N. McCall and Anthony Almudevar}

\seealso{
  \code{tnetfit}, \code{ternaryFitParameters-class}, \code{ternaryFit-class}.
  Almudevar A, McCall MN, McMurray H, Land H (2011). Fitting
  Boolean Networks from Steady State Perturbation Data, Statistical
  Applications in Genetics and Molecular Biology, 10(1): Article 47.  
} 

\examples{
ssObj <- matrix(c(1,1,1,0,1,1,0,0,1),nrow=3)
pObj <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)
rownames(ssObj) <- rownames(pObj) <- colnames(ssObj) <- colnames(pObj) <- c("Gene1","Gene2","Gene3")
tnfitObj <- tnetfit(ssObj, pObj)
tnpostObj <- tnetpost(tnfitObj, mdelta=10, msample=10)
class(tnpostObj)
}

\keyword{classes}