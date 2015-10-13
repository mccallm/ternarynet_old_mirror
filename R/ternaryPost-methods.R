ternaryPost <- function(perturbationObj, steadyStateObj, geneNames, experimentNames, scores, degreeObjs, graphObjs, tableObjs, inputParams){
  new("ternaryPost", perturbationObj=perturbationObj, steadyStateObj=steadyStateObj, geneNames=geneNames, experimentNames=experimentNames, scores=scores, degreeObjs=degreeObjs, graphObjs=graphObjs, tableObjs=tableObjs, inputParams=inputParams)
}

## dim method
setMethod("dim", "ternaryPost", function(x) dim(x@perturbationObj))

## slot getters
setMethod("perturbationObj", "ternaryPost", function(x) x@perturbationObj)
setMethod("steadyStateObj", "ternaryPost", function(x) x@steadyStateObj)
setMethod("geneNames", "ternaryPost", function(x) x@geneNames)
setMethod("experimentNames", "ternaryPost", function(x) x@experimentNames)
setMethod("scores", "ternaryPost", function(x) x@scores)
setMethod("degreeObjs", "ternaryPost", function(x) x@degreeObjs)
setMethod("graphObjs", "ternaryPost", function(x) x@graphObjs)
setMethod("tableObjs", "ternaryPost", function(x) x@tableObjs)
setMethod("inputParams", "ternaryPost", function(x) x@inputParams)

## show method
setMethod("show", "ternaryPost", function(object){
    cat(class(object), "instance with", dim(object)[2], "perturbation experiments, measuring", dim(object)[2], "genes \n")
    cat("Gene names:", geneNames(object),"\n")
    cat("Experiment names:", experimentNames(object),"\n")
    }) 

## validity method
setValidity("ternaryPost", function(object){
  TRUE
})

## slot setters
setReplaceMethod("perturbationObj", "ternaryPost", function(x, value){x@perturbationObj <- value; validObject(x); x})
setReplaceMethod("steadyStateObj", "ternaryPost", function(x, value){x@steadyStateObj <- value; validObject(x); x})
setReplaceMethod("geneNames", "ternaryPost", function(x, value){x@geneNames <- value; validObject(x); x})
setReplaceMethod("experimentNames", "ternaryPost", function(x, value){x@experimentNames <- value; validObject(x); x})
setReplaceMethod("scores", "ternaryPost", function(x, value){x@scores <- value; validObject(x); x})
setReplaceMethod("degreeObjs", "ternaryPost", function(x, value){x@degreeObjs <- value; validObject(x); x})
setReplaceMethod("graphObjs", "ternaryPost", function(x, value){x@graphObjs <- value; validObject(x); x})
setReplaceMethod("tableObjs", "ternaryPost", function(x, value){x@tableObjs <- value; validObject(x); x})
setReplaceMethod("inputParams", "ternaryPost", function(x, value){x@inputParams <- value; validObject(x); x})

