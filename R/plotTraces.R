plotTraces <- function(tfit){

  par(mfrow=c(2,2))
  n <- length(traces(tfit)$muTrace)
  plot(x=1:n, y=traces(tfit)$muTrace, type="l", xlab="stages", ylab="mu", main="muTrace")
  plot(x=1:n, y=traces(tfit)$sigmaTrace, type="l", xlab="stages", ylab="sigma", main="sigmaTrace")
  plot(x=1:n, y=traces(tfit)$sigmaRhoTrace, type="l", xlab="stages", ylab="sigmaRho", main="sigmaRhoTrace")
  plot(x=1:n, y=traces(tfit)$temperatureTrace, type="l", xlab="stages", ylab="temperature", main="temperatureTrace")

}
