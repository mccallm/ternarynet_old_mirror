# include <iostream>
# include <fstream>
# include <vector>
# include <stdio.h>
# include <cstdlib>
# include <iomanip>
# include "math.h"
# include <algorithm>
# include <vector>

# include <R.h>
# include <Rdefines.h>
# include <Rinternals.h>
# include <Rmath.h>

# include "tnetfuncs.h"

using namespace std;


//int globalSeed=0;

//int globalRand(){
//	int u;
//cout << "g1 " << globalSeed << endl;
//	u = ((35 * globalSeed) + 192283) % 1048576;
//	globalSeed = u;	
//cout << "g2 " << globalSeed << endl;
//	return u;
//}


// function to generate random number of uniform distribution
double uniformdist(){
  //  double unrand;
  //unrand = ((double)globalRand())/( (double)(RAND_MAX)+(double)(1) );
  //unrand = ((double)globalRand())/( (double)(1048576)+(double)(1) );	
  //return unrand;
  return unif_rand();
}

double stagemean(double *values, int n){
  int i;
  double average=0.0;
  for (i=1;i <= n;i++){
    average = average + values[i];
  }
  average = average/n;
  return average;
}

double stagevariance(double *values, int n){
  int i;
  double s, ss,v;
  s = 0;
  ss = 0; 
  for (i=1;i <= n;i++){
    s = s + values[i];
    ss = ss + pow(values[i],2);
  }
  v = ss/(n-1)-pow(s,2)/(pow(n,2)-n);
  return v;
}

int powi(int base, int exp) {
	int n=1;
	for (int i=1; i <= exp; ++i) n = n*base;
	return n;
} 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int randomNeighborDegree(int nNeighbor, const double *pN) {

	int k;
        bool flag;
        double u;
        k=1;
        flag = true;
        while (flag) {
	        if (k >= nNeighbor) {
                	flag = false;
		} else { 
                	u = uniformdist();
                      	if (u <= pN[k-1]) {flag = false;} else {++k;};
		}
	} 

/*
cout << "nNeighbor = " << nNeighbor << " pN[0] = " << pN[0] << 
" pN[1] = " << pN[1] << " pN[2] = " << pN[2] << " pN[3] = " << pN[3] << 
" k = " << k << endl;	
*/
	return k;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool IncrIndex(int *index, int idim, int ibase) {

	int iplace;
	
	iplace=idim;
	while ( (index[iplace] == ibase) && (iplace > 1) ) {
		index[iplace]=1;		
		--iplace;
	}

	if ( (iplace ==1) && ( index[1] == ibase ) ) {
		for (int i=1; i <= idim; ++i) index[i]=1;
		return false;
	} else
	{
		++index[iplace];
		return true;
	}		

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


int ArrayToHash(int *iDegree, int nDegree, int nOutcomes)
{
	int iTemp;
        iTemp=0;
        for (int i=1; i <= nDegree; ++i) iTemp += powi(nOutcomes,i-1)*(iDegree[i]-1);
        return iTemp+1;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SortInteger(int n, int *iarray) {

	int iswap;

	if (n > 1) {
		for (int i=1; i <= (n-1); ++i) {
                	for (int j=1; j <= (n-i); ++j) {
                        	if (iarray[j] > iarray[j+1]) {
                                	iswap = iarray[j];
                                        iarray[j] = iarray[j+1];
                                        iarray[j+1] = iswap;
                                }
                        }
                }
        }

}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int RandomGeneExclude(int *excludeGenes, int nExcluded, int nGene) {

	int u, n;

	n=nGene - nExcluded;
	u = (int)floor(unif_rand() * n) + 1; //add (int)floor() 01-18-2012
	SortInteger(nExcluded, excludeGenes);
	for (int i=1; i <= nExcluded; ++i) {
		if  (u >= excludeGenes[i]) ++u;
	}
  	return u;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TestDimension(int *newDegree, int nOutcomes, int *newGraphComponent, int *newTableComponent) {

	int newSize, hashIndex1, hashIndex2, newReducedDegree;
	newSize = newDegree[0]+1;
	bool hasReducibleDimension;
	bool *isReducible = new bool[newSize];
	int *index = new int [newSize];
	int *tempIndex = new int[newSize]; 
	int *indexMap = new int [newSize];	
	int *reducedIndex = new int [newSize];		
	int *newReducedTableComponent = new int [powi(nOutcomes,newDegree[0]) + 1];
	

	//for (int i=1; i <= powi(nOutcomes,oldDegree); ++i) newTableComponent[i] = oldTableComponent[i];

	for (int i=1; i <= newDegree[0]; ++i) {
		index[i]=1;
		isReducible[i] = true;
	}

	do {
		for (int i=1; i <= newDegree[0]; ++i) {
			if ( (index[i] > 1) & isReducible[i] ) {
				for (int j=1; j <= newDegree[0]; ++j) tempIndex[j]=index[j];
				tempIndex[i]=1;
				hashIndex1=ArrayToHash(tempIndex, newDegree[0], nOutcomes);
				hashIndex2=ArrayToHash(index, newDegree[0], nOutcomes);			
				isReducible[i] = (newTableComponent[hashIndex2]==newTableComponent[hashIndex1]);	
			}
		}
		
		hasReducibleDimension = false;
		for (int i=1; i <= newDegree[0]; ++i) hasReducibleDimension = (hasReducibleDimension || isReducible[i]); 
	
	} while (IncrIndex(index, newDegree[0], nOutcomes) && hasReducibleDimension);	

	newReducedDegree=newDegree[0];
	if (hasReducibleDimension) {

		newReducedDegree=0;
		for (int i=1; i <= newDegree[0]; ++i) {
			if (!isReducible[i]) {
				++newReducedDegree;
				indexMap[newReducedDegree]=i;
				newGraphComponent[newReducedDegree]=newGraphComponent[i];		
			}
		}			
		if (newReducedDegree > 0) {
			for (int i=1; i <= newDegree[0]; ++i) {
				index[i] = 1;
				reducedIndex[i] = 1;
			}
			do {
				for (int i=1; i <= newReducedDegree; ++i) index[indexMap[i]]=reducedIndex[i];
				hashIndex1=ArrayToHash(index, newDegree[0], nOutcomes);
				hashIndex2=ArrayToHash(reducedIndex, newReducedDegree, nOutcomes);			
				newReducedTableComponent[hashIndex2]=newTableComponent[hashIndex1];	
			} while (IncrIndex(reducedIndex, newReducedDegree, nOutcomes));				
		} else {
			newReducedTableComponent[1]=2;
		}		
		
		newDegree[0]=newReducedDegree;
		for (int i=1; i <= powi(nOutcomes,newDegree[0]); ++i) newTableComponent[i]=newReducedTableComponent[i]; 		
		
	}

	delete [] isReducible;
	delete [] index;
	delete [] tempIndex; 
	delete [] indexMap;	
	delete [] reducedIndex;		
	delete [] newReducedTableComponent;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void PerturbFunction(int *newDegree, int nOutcomes, int *newTableComponent) {

	// ternary specific

	int u;
	int *index = new int [newDegree[0] + 1];
	int *excludedList = new int [2];

	for (int i=1; i <= newDegree[0]; ++i) index[i]=2;

        excludedList[1]=ArrayToHash(index, newDegree[0], nOutcomes);
        u=RandomGeneExclude(excludedList, 1, powi(nOutcomes,newDegree[0]));

//cout << "bp1" << endl;
        if (newTableComponent[u] != 2) {
        	newTableComponent[u]=2;
//cout << "bp2" << endl;
        } else
        {
                newTableComponent[u]=3;
                if (uniformdist() < 0.5) newTableComponent[u]=1;
//cout << "bp3" << endl;        
	}

	delete [] index;
	delete [] excludedList;

}        

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void AddParent(int *newDegree, int proposalGene, int nGene, int nOutcomes, int *newGraphComponent, int *newTableComponent) 
{
	int tempDegree, newParent, hashIndex1, hashIndex2;
	int *excludedGenes = new int [nGene + 1];	
	int *index = new int [newDegree[0] + 2];
	int *tempIndex = new int [newDegree[0] + 2];

	// Add New Parent 

	for (int i=1; i <= newDegree[0]; ++i) excludedGenes[i]=newGraphComponent[i];
	excludedGenes[newDegree[0]+1]=proposalGene; 
	tempDegree=newDegree[0]+1;
        newParent=RandomGeneExclude(excludedGenes, tempDegree, nGene);
	newGraphComponent[newDegree[0]+1]=newParent;
	++newDegree[0];
	
	// Extend Function Table

	for (int i=1; i <= newDegree[0]; ++i) index[i]=1;
	do {
		hashIndex1=ArrayToHash(index, newDegree[0]-1, nOutcomes);
		hashIndex2=ArrayToHash(index, newDegree[0], nOutcomes);			
		newTableComponent[hashIndex2]=newTableComponent[hashIndex1];			
	} while (IncrIndex(index, newDegree[0], nOutcomes));	

      	PerturbFunction(newDegree, nOutcomes, newTableComponent);

	delete [] excludedGenes;	
	delete [] index;
	delete [] tempIndex;
	
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void PerturbGene(int proposalGene, int nGene, int maxDegree, int nOutcomes, 
		double probAddParent, double probExchangeParent,
		int *newDegree, int *newTableComponent, int *newGraphComponent)  
{		

	int tempDegree, perturbationType, newParent, oldParent;
	bool addParentFlag, exchangeParentFlag;
	int *tempTableComponent = new int [powi(nOutcomes,maxDegree) + 1]; 
	int *tempGraphComponent = new int [maxDegree + 1];
	int *excludedGenes = new int [maxDegree + 2];	

	addParentFlag = (uniformdist() < probAddParent) && (newDegree[0] < maxDegree);
        exchangeParentFlag = (uniformdist() < probExchangeParent) && (newDegree[0] < (nGene-1));  // 05-22-2013
	if ( (newDegree[0] > 0) &&  (exchangeParentFlag && !addParentFlag) ) perturbationType=1; 	// 09-15-2011
	if ( (newDegree[0] == 0) || (addParentFlag) )  perturbationType=2; // 09-15-2011
	if ( (newDegree[0] > 0) && !exchangeParentFlag && !addParentFlag ) perturbationType=3;


//cout << "proposalGene = " << proposalGene << endl;
//cout << "addParentFlag = " << addParentFlag << endl;
//cout << "exchangeParentFlag = " << exchangeParentFlag << endl;
//cout << "newDegree[0] = " << newDegree[0] << endl;
//cout << "probAddParent = " << probAddParent << endl;
//cout << "probExchangeParent = " << probExchangeParent << endl;
//cout << "perturbationType = " << perturbationType << endl;


        if ( perturbationType==1 ) { 

                // Exchange Parent
                
		for (int i=1; i <= newDegree[0]; ++i) excludedGenes[i]=newGraphComponent[i];
		excludedGenes[newDegree[0]+1]=proposalGene; 
		tempDegree=newDegree[0]+1;
                newParent=RandomGeneExclude(excludedGenes, tempDegree, nGene);
                oldParent= (int)floor(unif_rand() * newDegree[0]) + 1; //add (int)floor() 01-18-2012
		newGraphComponent[oldParent]=newParent;
        }
        if ( perturbationType==2 ) {
                
		// Add Parent

		AddParent(newDegree, proposalGene, nGene, nOutcomes, newGraphComponent, newTableComponent);
        }       
	if ( perturbationType==3 ) {
        
		// Perturb Function Table

                PerturbFunction(newDegree, nOutcomes, newTableComponent);
                TestDimension(newDegree, nOutcomes, newGraphComponent, newTableComponent);
        }

	delete [] tempTableComponent; 
	delete [] tempGraphComponent;
	delete [] excludedGenes;	

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double errorFunction(double steadyStateValue, int predictedValue) {

	// ternary specific 
	
	double err, pv;
	//int ipause;	

	pv = predictedValue;
	if (predictedValue == 4) {
		err = 1;
	} else {
		err = fabs(pv - steadyStateValue);
		if (err > 1) err = 1;
	}
	
	return err;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int NewOutput(int OldState, int TableInput) {

	// increment changes in state 

	int NewState =0;
	int nIncr = 2; // 08-06-2012
	
	if (TableInput == 3 && OldState < nIncr) {
		NewState = OldState + 1;
	}
	if (TableInput == 1 && OldState > -nIncr ) {
		NewState = OldState -1;
	}		
	if (TableInput == 2 && OldState < 0) {
		NewState = OldState + 1;
	}		
	if (TableInput == 2 && OldState > 0) {
		NewState = OldState - 1;
	}		

	return NewState;

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			
void ApplyOp(int nGene, int nOutcomes, int maxDegree, int *graphObj, int *tableObj, int *degreeObj, int *pathwayIn, int *pathwayOut) {

// ternary specific 

// we don't need nIncr here. 

int tableWidth = powi(nOutcomes,maxDegree);
int hashIndex;
int *index = new int[nGene+1];

for (int i = 1; i <= nGene; ++i) {

	if (degreeObj[i] == 0) {
		
		pathwayOut[i]=NewOutput(pathwayIn[i], 2); 
		// pathwayOut[i]=2;
		
	} else {
		for (int j = 1; j<=degreeObj[i]; ++j) {
			// 08-15-2012 change index[] index from i to j
			index[j] = 2;
			if (pathwayIn[graphObj[nGene*(i-1)+j-1]] < 0) {index[j] = 1;}; 			
			if (pathwayIn[graphObj[nGene*(i-1)+j-1]] > 0) {index[j] = 3;}; // 
			// index[j]=pathwayIn[graphObj[nGene*(i-1)+j-1]];
		}
		hashIndex=ArrayToHash(index, degreeObj[i], nOutcomes);

		pathwayOut[i]=NewOutput(pathwayIn[i], tableObj[tableWidth*(i-1)+hashIndex-1]); 
                // pathwayOut[i] = tableObj[tableWidth*(i-1)+hashIndex-1];
        }
}

//delete hashIndex;
delete [] index;

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double AttractorDistanceForced(int nGene, int nOutcomes, int maxDegree, int nExperiment, double edgePenalty, 
	int *tableObj, int *graphObj, int *degreeObj, 
	const double *steadyStateObj, const int *perturbationObj) {

// ternary specific
        
int nIncr = 2; // 08-06-2012
int tempState = 0; // 08-06-2012

vector<int> pathwayList;
int tableWidth = powi(nOutcomes,maxDegree);
int *pathwayIn = new int[nGene+1];
int *pathwayOut = new int[nGene+1];
int *attractorSummary = new int[nGene+1];

int ListPos, ListSize, iTemp, jTemp;
double minScore;

int ipause;

minScore=0;

for (int iExperiment = 1; iExperiment <= nExperiment; ++iExperiment) {

        // clear and initialize pathway 

	pathwayList.clear();
        for (int i=1; i <= nGene; ++i) {
		
		tempState = 0;
		if (perturbationObj[(i-1)*nExperiment + iExperiment-1] == 1) {tempState = -nIncr;};
		if (perturbationObj[(i-1)*nExperiment + iExperiment-1] == 3) {tempState = nIncr;};
				
		pathwayList.push_back(tempState); // debug 09-21-2011
		pathwayIn[i]=tempState; // debug 09-21-2011
	}		

	ListPos=0; ListSize=1;

	while (ListPos == 0) {
        	++ListSize;
               ApplyOp(nGene, nOutcomes, maxDegree, graphObj, tableObj, degreeObj, pathwayIn, pathwayOut);			

		// apply perturbation

               for (int i=1; i <= nGene; ++i) {
		 if (perturbationObj[(i-1)*nExperiment + iExperiment-1] != 2) { // debug 09-21-2011
			tempState = 0;
			if (perturbationObj[(i-1)*nExperiment + iExperiment-1] == 1) {tempState = -nIncr;};
			if (perturbationObj[(i-1)*nExperiment + iExperiment-1] == 3) {tempState = nIncr;};
			pathwayOut[i]=tempState; // debug 09-21-2011
		 }
		 pathwayIn[i]=pathwayOut[i];
		 pathwayList.push_back(pathwayOut[i]);
	       }
		
		// search for cycle event

                iTemp=1;
                while ( (ListPos==0) && (iTemp < ListSize) ) {
                        ListPos=iTemp;
                        jTemp=1;
		        while ( (ListPos > 0) && (jTemp <= nGene) ) {
                                if (pathwayList[nGene*(iTemp-1)+jTemp-1] != pathwayOut[jTemp]) ListPos=0;
                                jTemp=jTemp+1;
                        }
                        iTemp=iTemp+1;
                }
        }

        for (int j=1; j <= nGene; ++j) attractorSummary[j]=2;
        for (int i=ListPos; i <= (ListSize-1); ++i) {
               	for (int j=1; j <= nGene; ++j) {
			if (attractorSummary[j] != 4) {
	                	if (pathwayList[nGene*(i-1)+j-1] < 0) {
        	                	if (attractorSummary[j] == 3 ) {attractorSummary[j]=4;}  else {attractorSummary[j]=1;}}
                        	if (pathwayList[nGene*(i-1)+j-1] > 0) {
        	                	if (attractorSummary[j] == 1 ) {attractorSummary[j]=4;}  else {attractorSummary[j]=3;}}
			}                                                
                }
	}
  
        for (int j=1; j <= nGene; ++j) minScore = minScore + errorFunction(steadyStateObj[(j-1)*nExperiment + iExperiment-1], attractorSummary[j]); 

}

// add complexity penalty 05-01-2011

for (int j=1; j <= nGene; ++j) minScore = minScore + edgePenalty*degreeObj[j]; 

pathwayList.clear();

delete [] pathwayIn;
delete [] pathwayOut;
delete [] attractorSummary;

return minScore;

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double initialTemp(double chi0, int m0,
	int nOutcomes, int maxDegree, int nGene, int nExperiment, double edgePenalty, 
	double pAddParent, double pExchangeParent, int neighborDegree, const double *pNeighborhood,
	int *degreeObj, int *graphObj, int *tableObj,
	const double *steadyStateObj, const int *perturbationObj) 
{
	int ipause; 

	double temperature, difference;
  	double d, dtemp, acceptprob;
  	int m1, m2, i;
	double newScore, oldScore;
	int u;


	int tableWidth = powi(nOutcomes,maxDegree);

	int proposalDegree, iProposal;
	int *proposalGene = new int[neighborDegree+1];
	int *oldDegree = new int[neighborDegree+1];
 	int *newDegree = new int[1];
	int *oldGraphComponent = new int[neighborDegree*nGene];
	int *newGraphComponent = new int[nGene+1];
	int *oldTableComponent = new int[neighborDegree*tableWidth]; 
	int *newTableComponent = new int[tableWidth+1]; 

  	temperature = 0;
  	m1 = 0;
  	m2 = 0;
  	d = 0;

	/*
	  cout << "steadyStateObj = " << endl; 
	  for (int ii = 1; ii <= nGene; ++ii) {
	  for (int jj = 1; jj <= nExperiment; jj++) {
	  cout << steadyStateObj[nExperiment*(ii-1) + jj-1] << " ";	
	  }
	  cout << endl;
	  }	
	  cin >> ipause;
	*/

	oldScore = AttractorDistanceForced(nGene, nOutcomes, nExperiment, maxDegree, edgePenalty,
		tableObj, graphObj, degreeObj, steadyStateObj, perturbationObj); 
			
  	for (int i=1; i<=m0; i++){
	
		proposalDegree = randomNeighborDegree(neighborDegree, pNeighborhood);

//cout << "proposalDegree = " << proposalDegree << endl;
//cout << "neighborDegree = " << neighborDegree << endl;
//cout << "pNeighborhood[0] = " << pNeighborhood[0] << endl; 

 		for (iProposal = 1; iProposal <= proposalDegree; ++iProposal) {
		  u = unif_rand();
//cout << "u = " << u << endl; 
			proposalGene[iProposal] = (int)floor(u * nGene) + 1; //add (int)floor() 01-18-2012
//cout << "proposalGene[iProposal] = " << proposalGene[iProposal] << endl; 
			oldDegree[iProposal] = degreeObj[proposalGene[iProposal]];
			//cout << "BP1" << endl;
			newDegree[0] = degreeObj[proposalGene[iProposal]];
			//cout << "BP2" << endl;
			for (int j=1; j <= oldDegree[iProposal]; ++j) {
				oldGraphComponent[(iProposal-1)*nGene+j-1]=graphObj[nGene*(proposalGene[iProposal]-1)+j-1];
				newGraphComponent[j]=graphObj[nGene*(proposalGene[iProposal]-1)+j-1];
			}		
			//cout << "BP3" << endl;		
			for (int j=1; j <= powi(nOutcomes, oldDegree[iProposal]); ++j) {
				oldTableComponent[(iProposal-1)*tableWidth+j-1]=tableObj[(proposalGene[iProposal]-1)*tableWidth+j-1];
				newTableComponent[j]=tableObj[(proposalGene[iProposal]-1)*tableWidth+j-1]; 		
			}
			//cout << "BP4" << endl;
			PerturbGene(proposalGene[iProposal], nGene, maxDegree, nOutcomes, 
				pAddParent, pExchangeParent, newDegree, newTableComponent, newGraphComponent);
			//cout << "BP5" << endl;
			degreeObj[proposalGene[iProposal]] = newDegree[0];
			//cout << "BP6" << endl;
			for (int j=1; j <= newDegree[0]; ++j) {
				graphObj[nGene*(proposalGene[iProposal]-1)+j-1]=newGraphComponent[j];}
			//cout << "BP7" << endl;
			for (int j=1; j <= powi(nOutcomes, newDegree[0]); ++j) {
				tableObj[(proposalGene[iProposal]-1)*tableWidth+j-1] = newTableComponent[j];}
			//cout << "BP8" << endl;
		}
		/*
		  cout << "Pre-accept (R) = " << endl; 
		  cout << "graphObj (R) = " << endl; 
		  for (int ii = 1; ii <= nGene; ++ii) {
		  for (int jj = 1; jj <= nGene; jj++) {
		  //cout << graphObj[nGene*(ii-1) + jj-1] << " ";	
		  }
		  //cout << endl;
		  }	
		  //cout << "degreeObj (non R) = " << endl;
		  for (int ii = 1; ii <= nGene; ++ii) {//cout << degreeObj[ii] << " ";}
		  //cout << endl;
		  //cout << "tableObj (non R) = " << endl; 
		  for (int ii = 1; ii <= nGene; ++ii) {
		  for (int jj = 1; jj <= tableWidth; jj++) {
		  //cout << tableObj[tableWidth*(ii-1) + jj-1] << " ";	
		  }
		  //cout << endl;
		  }	
		  //cin >> ipause;
		  */
		
		newScore = AttractorDistanceForced(nGene, nOutcomes, maxDegree, nExperiment, edgePenalty, 
		tableObj, graphObj, degreeObj, steadyStateObj, perturbationObj); 

		//cout << "newScore = " << newScore << endl;

		difference = newScore - oldScore;

		if (difference <= 0){
      			m1 = m1 + 1; 
		}
    		else{
      			m2 = m2 + 1;
      			d =  ((m2 - 1) * d + fabs(difference)) / m2; 
			//cout << "temp:" << temperature << endl;
      			if (temperature != 0){
        			acceptprob = exp((-difference)/temperature); 
//cout << "bp4" << endl;
        			if (acceptprob < uniformdist()){
//cout << "bp5" << endl;
					newScore = oldScore;
					for (iProposal = proposalDegree; iProposal >= 1; --iProposal) {
						degreeObj[proposalGene[iProposal]] = oldDegree[iProposal];
						for (int j=1; j <= oldDegree[iProposal]; ++j) {
							graphObj[(proposalGene[iProposal]-1)*nGene+j-1] 
								= oldGraphComponent[(iProposal-1)*nGene+j-1];}
						for (int j=1; j <= powi(nOutcomes, oldDegree[iProposal]); ++j) { 
							tableObj[(proposalGene[iProposal]-1)*tableWidth+j-1] 
								= oldTableComponent[(iProposal-1)*tableWidth+j-1];}
					}
				}
      			}
    		}    
		/*
//cout << "difference= " << difference << " m1= " << m1 << " m2= " << m2 << " temperature= " << temperature << endl;
//cout << "Post-accept (R) = " << endl; 
//cout << "graphObj (R) = " << endl; 
for (int ii = 1; ii <= nGene; ++ii) {
	for (int jj = 1; jj <= nGene; jj++) {
			//cout << graphObj[nGene*(ii-1) + jj-1] << " ";	
	}
	//cout << endl;
}	
//cout << "degreeObj (non R) = " << endl;
for (int ii = 1; ii <= nGene; ++ii) {//cout << degreeObj[ii] << " ";}
//cout << endl;
//cout << "tableObj (non R) = " << endl; 
for (int ii = 1; ii <= nGene; ++ii) {
	for (int jj = 1; jj <= tableWidth; jj++) {
			//cout << tableObj[tableWidth*(ii-1) + jj-1] << " ";	
	}
	//cout << endl;
}	
//cin >> ipause;
*/

		//cout << "m2=" << m2 << endl;
		//cout << "m1=" << m1 << endl;
		//cout << "chi0=" << chi0 << endl;

		// update temperature
    		dtemp = m2*chi0 - m1*(1 - chi0);

		//cout << "dtemp=" << dtemp << endl;

    		if (m2 > 0 && dtemp > 0) {
      			temperature = d / (log(m2 / dtemp));
    		}
		oldScore=newScore;
		

  	}

	//cout << "did we make it" << endl;
	
	delete [] proposalGene;
	delete [] oldDegree;
 	delete [] newDegree;
	delete [] oldGraphComponent;
	delete [] newGraphComponent;
	delete [] oldTableComponent; 
	delete [] newTableComponent; 

  return temperature;
}
