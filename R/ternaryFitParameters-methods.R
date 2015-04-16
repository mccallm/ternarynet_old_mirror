ternaryFitParameters <- function(perturbationType=as.integer(1),
                                 scoreType = as.integer(2),
                                 backupStage = as.integer(1),
                                 maxStage = as.integer(2000), 
                                 maxTransition = as.integer(2000),
                                 epsilon = 1e-05,
                                 beta0 = as.integer(3), 
                                 chi0 = 0.95,
                                 delta = 0.1,
                                 ne = as.integer(10000),
                                 m0 = as.integer(10000), 
                                 maxDegree = as.integer(1),
                                 pAddParent = 0.1,
                                 pExchangeParent = 0.1, 
                                 neighborDegree = as.integer(2),
                                 pNeighborhood = 0.01,
                                 rho = 0.25, 
                                 edgePenalty = 0){
  new("ternaryFitParameters", perturbationType=perturbationType, scoreType=scoreType, backupStage=backupStage, maxStage=maxStage, maxTransition=maxTransition, epsilon=epsilon, beta0=beta0, chi0=chi0, delta=delta, ne=ne, m0=m0, maxDegree=maxDegree, pAddParent=pAddParent, pExchangeParent=pExchangeParent, neighborDegree=neighborDegree, pNeighborhood=pNeighborhood, rho=rho, edgePenalty=edgePenalty)
}

## slot getters
setMethod("perturbationType", "ternaryFitParameters", function(x) x@perturbationType)
setMethod("scoreType", "ternaryFitParameters", function(x) x@scoreType)
setMethod("backupStage", "ternaryFitParameters", function(x) x@backupStage)
setMethod("maxStage", "ternaryFitParameters", function(x) x@maxStage)
setMethod("maxTransition", "ternaryFitParameters", function(x) x@maxTransition)
setMethod("epsilon", "ternaryFitParameters", function(x) x@epsilon)
setMethod("beta0", "ternaryFitParameters", function(x) x@beta0)
setMethod("chi0", "ternaryFitParameters", function(x) x@chi0)
setMethod("delta", "ternaryFitParameters", function(x) x@delta)
setMethod("ne", "ternaryFitParameters", function(x) x@ne)
setMethod("m0", "ternaryFitParameters", function(x) x@m0)
setMethod("maxDegree", "ternaryFitParameters", function(x) x@maxDegree)
setMethod("pAddParent", "ternaryFitParameters", function(x) x@pAddParent)
setMethod("pExchangeParent", "ternaryFitParameters", function(x) x@pExchangeParent)
setMethod("neighborDegree", "ternaryFitParameters", function(x) x@neighborDegree)
setMethod("pNeighborhood", "ternaryFitParameters", function(x) x@pNeighborhood)
setMethod("rho", "ternaryFitParameters", function(x) x@rho)
setMethod("edgePenalty", "ternaryFitParameters", function(x) x@edgePenalty)

## show method
setMethod("show", "ternaryFitParameters", function(object) cat(class(object), "instance \n")) 

## validity method
setValidity("ternaryFitParameters", function(object){
  ints <- c(perturbationType(object), scoreType(object), backupStage(object), maxStage(object), maxTransition(object), beta0(object), ne(object), m0(object), maxDegree(object), neighborDegree(object))
  if(length(ints) != 10 || any(!is.integer(ints)) || any(is.na(ints))) return("'perturbationType', 'scoreType', 'backupStage', 'maxStage', 'maxTransition', 'beta0', 'ne', 'm0', 'maxDegree', and 'neighborDegree' must all be single integers")
  nums <- c(epsilon(object), chi0(object), delta(object), pAddParent(object), pExchangeParent(object), rho(object), edgePenalty(object))
  if(length(nums) != 7 || any(!is.numeric(nums)) || any(is.na(nums))) return("'epsilon', 'chi0', 'delta', 'pAddParent', 'pExchangeParent', 'rho', and 'edgePenalty' must all be single numeric values")
  TRUE
})

## slot setters
setReplaceMethod("perturbationType", "ternaryFitParameters", function(x, value){x@perturbationType <- value; validObject(x); x})
setReplaceMethod("scoreType", "ternaryFitParameters", function(x, value){x@scoreType <- value; validObject(x); x})
setReplaceMethod("backupStage", "ternaryFitParameters", function(x, value){x@backupStage <- value; validObject(x); x})
setReplaceMethod("maxStage", "ternaryFitParameters", function(x, value){x@maxStage <- value; validObject(x); x})
setReplaceMethod("maxTransition", "ternaryFitParameters", function(x, value){x@maxTransition <- value; validObject(x); x})
setReplaceMethod("epsilon", "ternaryFitParameters", function(x, value){x@epsilon <- value; validObject(x); x})
setReplaceMethod("beta0", "ternaryFitParameters", function(x, value){x@beta0 <- value; validObject(x); x})
setReplaceMethod("chi0", "ternaryFitParameters", function(x, value){x@chi0 <- value; validObject(x); x})
setReplaceMethod("delta", "ternaryFitParameters", function(x, value){x@delta <- value; validObject(x); x})
setReplaceMethod("ne", "ternaryFitParameters", function(x, value){x@ne <- value; validObject(x); x})
setReplaceMethod("m0", "ternaryFitParameters", function(x, value){x@m0 <- value; validObject(x); x})
setReplaceMethod("maxDegree", "ternaryFitParameters", function(x, value){x@maxDegree <- value; validObject(x); x})
setReplaceMethod("pAddParent", "ternaryFitParameters", function(x, value){x@pAddParent <- value; validObject(x); x})
setReplaceMethod("pExchangeParent", "ternaryFitParameters", function(x, value){x@pExchangeParent <- value; validObject(x); x})
setReplaceMethod("neighborDegree", "ternaryFitParameters", function(x, value){x@neighborDegree <- value; validObject(x); x})
setReplaceMethod("pNeighborhood", "ternaryFitParameters", function(x, value){x@pNeighborhood <- value; validObject(x); x})
setReplaceMethod("rho", "ternaryFitParameters", function(x, value){x@rho <- value; validObject(x); x})
setReplaceMethod("edgePenalty", "ternaryFitParameters", function(x, value){x@edgePenalty <- value; validObject(x); x})
