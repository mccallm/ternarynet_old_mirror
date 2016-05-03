getObservedAttractor <- function(perturbations, tableObj, graphObj, degreeObj, wildtype){
  fg <- perturbations$perturbed.genes
  fs <- perturbations$forced.states + 2
  tn.model <- list(table = tableObj+2, graph = graphObj, degree = degreeObj)
  vec0 <- rep(2, length(tn.model$degree))
  vec0[fg] <- fs
  junk <- attractor.distance.summary(vec0, tn.model, wildtype, fg, fs)
  nc <- dim(junk$attractor)[2]
  return(summarize.attractor.2(junk$attractor[,2:nc]))
}

simulateSteadyState <- function(perturbationObj, tableObj, graphObj, degreeObj, wildtype=FALSE){
  steadyStateObj <- matrix(nrow=nrow(perturbationObj), ncol=ncol(perturbationObj))
  for(k in 1:ncol(perturbationObj)){
    ind <- which(perturbationObj[,k]!=0)
    perts <- list("perturbed.genes"=ind, "forced.states"=perturbationObj[ind,k])
    steadyStateObj[,k] <- getObservedAttractor(perts, tableObj, graphObj, degreeObj, wildtype)
  }
  return(steadyStateObj)               
}
