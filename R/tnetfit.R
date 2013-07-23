tnetfit <- function(steadyStateObj, perturbationObj, params=ternaryFitParameters(), xSeed=NA){

  if(!identical(rownames(steadyStateObj),rownames(perturbationObj))) stop("steadyStateObj and perturbationObj must have identical rownames")
  if(!identical(colnames(steadyStateObj),colnames(perturbationObj))) stop("steadyStateObj and perturbationObj must have identical colnames")  
    
  if(is.null(rownames(steadyStateObj))) gNames <- paste("Gene",1:nrow(steadyStateObj),sep="") else gNames <- rownames(steadyStateObj)
  if(is.null(colnames(steadyStateObj))) expNames <- paste("Exp",1:ncol(steadyStateObj),sep="") else expNames <- colnames(steadyStateObj)
    
  nGene <- as.integer(nrow(steadyStateObj))
  nExperiment <- as.integer(ncol(steadyStateObj))

  ssObj <- as.vector(t(steadyStateObj))+2
  pObj <- as.integer(t(perturbationObj)+2)
  
  if(is.na(xSeed)) xSeed <- as.integer(floor(runif(1)*10000))
  set.seed(xSeed)
    
  tmp <- .Call("tnetfit", as.integer(perturbationType(params)), as.integer(scoreType(params)), as.integer(backupStage(params)), as.integer(maxStage(params)), as.integer(maxTransition(params)), epsilon(params), as.integer(beta0(params)), chi0(params), delta(params), as.integer(ne(params)), as.integer(m0(params)), as.integer(maxDegree(params)), pAddParent(params), pExchangeParent(params), as.integer(neighborDegree(params)), pNeighborhood(params), rho(params), nGene, nExperiment, edgePenalty(params), ssObj, pObj)

  traces <- data.frame(muTrace=tmp$muTrace, sigmaTrace=tmp$sigmaTrace, sigmaRhoTrace=tmp$sigmaRhoTrace, temperatureTrace=tmp$temperatureTrace)

  ternaryFit(perturbationObj=perturbationObj, steadyStateObj=steadyStateObj, geneNames=gNames, experimentNames=expNames, degreeObjMin=tmp$degreeObjMin, graphObjMin=matrix(tmp$graphObjMin, nrow=nGene), tableObjMin=matrix(tmp$tableObjMin, ncol=nGene), newScore=tmp$newScore, minScore=tmp$minScore, finalTemperature=tmp$scTemperature, traces=traces, stageCount=tmp$stageCount, xSeed=as.integer(xSeed), inputParams=params)
  
}
