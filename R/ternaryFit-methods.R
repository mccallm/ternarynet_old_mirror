ternaryFit <- function(perturbationObj, steadyStateObj, degreeObjMin, graphObjMin, tableObjMin, newScore, minScore, finalTemperature, traces, stageCount, xSeed, inputParams){
  new("ternaryFit", perturbationObj=perturbationObj, steadyStateObj=steadyStateObj, degreeObjMin=degreeObjMin, graphObjMin=graphObjMin, tableObjMin=tableObjMin, newScore=newScore, minScore=minScore, finalTemperature=finalTemperature, traces=traces, stageCount=stageCount, xSeed=xSeed, inputParams=inputParams)
}

## dim method
setMethod("dim", "ternaryFit", function(x) dim(x@perturbationObj))

## slot getters
setMethod("perturbationObj", "ternaryFit", function(x) x@perturbationObj)
setMethod("steadyStateObj", "ternaryFit", function(x) x@steadyStateObj)
setMethod("degreeObjMin", "ternaryFit", function(x) x@degreeObjMin)
setMethod("graphObjMin", "ternaryFit", function(x) x@graphObjMin)
setMethod("tableObjMin", "ternaryFit", function(x) x@tableObjMin)
setMethod("newScore", "ternaryFit", function(x) x@newScore)
setMethod("minScore", "ternaryFit", function(x) x@minScore)
setMethod("finalTemperature", "ternaryFit", function(x) x@finalTemperature)
setMethod("xSeed", "ternaryFit", function(x) x@xSeed)
setMethod("inputParams", "ternaryFit", function(x) x@inputParams)
setMethod("traces", "ternaryFit", function(x) x@traces)
setMethod("stageCount", "ternaryFit", function(x) x@stageCount)

## show method
setMethod("show", "ternaryFit", function(object) cat(class(object), "instance with", dim(object)[2], "perturbation experiments, measuring", dim(object)[2], "genes \n")) 

## validity method
setValidity("ternaryFit", function(object){
  nGene <- nrow(perturbationObj(object))
  nExperiment <- ncol(perturbationObj(object))

  if(nrow(steadyStateObj(object)) != nGene || ncol(steadyStateObj(object)) != nExperiment) return("'steadyStateObj' slot and 'pertubationObj' slot must have the same dimensions")
  indp <- which(perturbationObj(object)!=0, arr.ind=TRUE)
  if(any(steadyStateObj(object)[indp] != perturbationObj(object)[indp])) return("Non-zero elements of 'perturbationObj' slot must match corresponding elements of 'steadyStateObj'")
  if(length(degreeObjMin(object)) != nGene) return("'degreeObjMin' slot must have length equal to the number of rows of the 'perturbationObj' slot")
  if(nrow(graphObjMin(object)) != nGene || ncol(graphObjMin(object)) != nGene) return("'graphObjMin' slot must be a matrix with row and column lengths equal to the number of rows of the 'perturbationObj' slot")
  if(ncol(tableObjMin(object)) != nGene) return("'tableObjMin' slot must be a matrix with row length equal to the number of rows of the 'perturbationObj' slot")
  if(!is.numeric(newScore(object)) || length(newScore(object)) != 1 || is.na(newScore(object))) return("'newScore' slot must be a single numeric")
  if(!is.numeric(minScore(object)) || length(minScore(object)) != 1 || is.na(minScore(object))) return("'minScore' slot must be a single numeric")
  if(!is.numeric(finalTemperature(object)) || length(finalTemperature(object)) != 1 || is.na(finalTemperature(object))) return("'finalTemperature' slot must be a single numeric")
  if(!is.integer(xSeed(object)) || length(xSeed(object)) != 1 || is.na(xSeed(object))) return("'xSeed' slot must be a single integer")
  if(class(inputParams(object)) != "ternaryFitParameters") return("'inputParams' slot must be an object of class 'ternaryFitParameters")
  TRUE
})

## slot setters
setReplaceMethod("perturbationObj", "ternaryFit", function(x, value){x@perturbationObj <- value; validObject(x); x})
setReplaceMethod("steadyStateObj", "ternaryFit", function(x, value){x@steadyStateObj <- value; validObject(x); x})
setReplaceMethod("degreeObjMin", "ternaryFit", function(x, value){x@degreeObjMin <- value; validObject(x); x})
setReplaceMethod("graphObjMin", "ternaryFit", function(x, value){x@graphObjMin <- value; validObject(x); x})
setReplaceMethod("tableObjMin", "ternaryFit", function(x, value){x@tableObjMin <- value; validObject(x); x})
setReplaceMethod("newScore", "ternaryFit", function(x, value){x@newScore <- value; validObject(x); x})
setReplaceMethod("minScore", "ternaryFit", function(x, value){x@minScore <- value; validObject(x); x})
setReplaceMethod("finalTemperature", "ternaryFit", function(x, value){x@finalTemperature <- value; validObject(x); x})
setReplaceMethod("xSeed", "ternaryFit", function(x, value){x@xSeed <- value; validObject(x); x})
setReplaceMethod("inputParams", "ternaryFit", function(x, value){x@inputParams <- value; validObject(x); x})
setReplaceMethod("traces", "ternaryFit", function(x, value){x@traces <- value; validObject(x); x})
setReplaceMethod("stageCount", "ternaryFit", function(x, value){x@stageCount <- value; validObject(x); x})

