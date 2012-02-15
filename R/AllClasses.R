setClass("ternaryFitParameters",
         representation(
            perturbationType="integer", # a single integer value
            scoreType="integer",        # a single integer value
            backupStage="integer",      # a single integer value
            maxStage="integer",         # a single integer value
            maxTransition="integer",    # a single integer value
            epsilon="numeric",          # a single numeric value
            beta0="integer",             # a single integer value
            chi0="numeric",             # a single numeric value
            delta="numeric",            # a single numeric value
            ne="integer",               # a single integer value
            m0="integer",               # a single integer value
            maxDegree="integer",        # a single integer value
            pAddParent="numeric",       # a single numeric value
            pExchangeParent="numeric",  # a single numeric value
            neighborDegree="integer",   # a single integer value
            pNeighborhood="numeric",    # a numeric vector of length ??
            rho="numeric",              # a single numeric value
            edgePenalty="numeric"      # a single numeric value
         )
)

setClass("ternaryFit",
         representation(
            perturbationObj="matrix",   # an nGene by nExperiment matrix
            steadyStateObj="matrix",    # an nGene by nExperiment matrix
            degreeObjMin="integer",     # an integer vector of length nGene
            graphObjMin="matrix",       # an nGene by nGene matrix
            tableObjMin="matrix",       # an nGene by tableWidth matrix
            newScore="numeric",         # a single numeric value
            minScore="numeric",         # a single numeric value
            finalTemperature="numeric", # a single numeric value
            traces="data.frame",        # a data frame of the fit traces
            stageCount="integer",       # a single integer value            
            xSeed="integer",            # a single integer value
            inputParams="ternaryFitParameters"
                                        # class containing input parameters
         )                        
)              

setClass("ternaryPost",
         representation(
            perturbationObj="matrix",   # an nGene by nExperiment matrix
            steadyStateObj="matrix",    # an nGene by nExperiment matrix
            scores="numeric",           # numeric vector
            degreeObjs="matrix",        # integer matrix
            graphObjs="array",          # integer array
            tableObjs="array",          # integer array
            inputParams="ternaryFitParameters"
                                        # class containing input parameters
         )                        
)              

                       
