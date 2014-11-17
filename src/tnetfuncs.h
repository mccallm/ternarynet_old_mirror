#ifndef TNET_HEADER
#define TNET_HEADER 1

// extern int globalSeed;
// int globalRand();
double uniformdist();
double stagemean(double *values, int n);
double stagevariance(double *values, int n);
int powi(int base, int exp);
int randomNeighborDegree(int nNeighbor, const double *pN);
bool IncrIndex(int *index, int idim, int ibase);
int ArrayToHash(int *iDegree, int nDegree, int nOutcomes);
void SortInteger(int n, int *iarray);
int RandomGeneExclude(int *excludeGenes, int nExcluded, int nGene);
void TestDimension(int *newDegree, int nOutcomes, int *newGraphComponent, int *newTableComponent);
void PerturbFunction(int *newDegree, int nOutcomes, int *newTableComponent);
void AddParent(int *newDegree, int proposalGene, int nGene, int nOutcomes, int *newGraphComponent, int *newTableComponent);
void PerturbGene(int proposalGene, int nGene, int maxDegree, int nOutcomes, double probAddParent, double probExchangeParent, int *newDegree, int *newTableComponent, int *newGraphComponent);
double errorFunction(double steadyStateValue, int predictedValue);
void ApplyOp(int nGene, int nOutcomes, int maxDegree, int *graphObj, int *tableObj, int *degreeObj, int *pathwayIn, int *pathwayOut);
double AttractorDistanceForced(int nGene, int nOutcomes, int maxDegree, int nExperiment, double edgePenalty, int *tableObj, int *graphObj, int *degreeObj, const double *steadyStateObj, const int *perturbationObj);
double initialTemp(double chi0, int m0,	int nOutcomes, int maxDegree, int nGene, int nExperiment, double edgePenalty, double pAddParent, double pExchangeParent, int neighborDegree, const double *pNeighborhood, int *degreeObj, int *graphObj, int *tableObj, const double *steadyStateObj, const int *perturbationObj);

#endif
