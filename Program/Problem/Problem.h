// *******************************************************************
//      file with specific functions to solve the TSP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Global Variables
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     k - number of local search heuristics
extern int n;                               // size of the vector solution

//---------------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC --------------------

// struct with node informations
struct TNode								
{
	int id;
	double x;
	double y;
}; 

//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  ----------

static std::vector <std::vector <double> > dist;	// matrix with Euclidean distance
static std::vector <TNode> node;					// vector of TSP nodes


//----------------------- IMPLEMENTATION OF FUNCTIONS  -------------------------------


/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
*************************************************************************************/
void ReadData(char nameTable[])
{ 
    char name[200] = "../Instances/";
    strcat(name,nameTable);

    FILE *arq;
    arq = fopen(name,"r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n",name);
        getchar();
        exit(1);
    }

    // => read data

    // read instance head
    char temp[100];
    fgets(temp, sizeof(temp), arq);
    fgets(temp, sizeof(temp), arq);
    fgets(temp, sizeof(temp), arq);
    fgets(temp, sizeof(temp), arq);
    fgets(temp, sizeof(temp), arq);
    fgets(temp, sizeof(temp), arq);

    // read node informations
    int nAux = 0;
    node.clear();
    TNode nodeTemp;

    while (!feof(arq))
    {
    	fscanf(arq, "%d %lf %lf", &nodeTemp.id, &nodeTemp.x, &nodeTemp.y);
    	node.push_back(nodeTemp);

    	nAux++;
    }
    fclose(arq);

    // calculate the euclidean distance
    dist.clear();
    dist.resize(nAux, std::vector<double>(nAux));

    for (int i=0; i<nAux; i++)
    {
    	for (int j=i; j<nAux; j++)
    	{
    		dist[i][j] = dist[j][i] = (floor (sqrt( (node[j].x - node[i].x) * (node[j].x - node[i].x) +
    										        (node[j].y - node[i].y) * (node[j].y - node[i].y) ) + 0.5 ) )/1.0;
    	}
    }
    
    n = nAux;
}

/************************************************************************************
 Method: CalculateFitness
 Description: calculate the fitness of a chromossom s
*************************************************************************************/
double CalculateFitness(TSol s)
{
    // calculate objective function
    s.ofv = 0;
    for (int i=0; i<n; i++){
        s.ofv += dist[s.vec[i%n].sol][s.vec[(i+1)%n].sol];
    }

    return s.ofv;
}

void Dec1(TSol &s) // Cheapest Insertion
{
    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < n)
		    s.vec[j].sol = j;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);

    TVecSol aux = s.vec[n];

    // order list of candidates
    TSol sC = s;

    // partial route with three points
    s.vec.resize(3);

    // construct a solution with cheapest insertion
    for (int i = 3; i<n; i++)
    {
        // find the cheapest position to insert the i-th point of sC
        int bestPosition = 0;
        float costBest = INFINITO;
        float costInsertion = 0;
        for (unsigned int j = 0; j<s.vec.size(); j++)
        {
            if (j == s.vec.size()-1)
            {
                // cost to insert between i-1 and 0
                costInsertion = dist[s.vec[j].sol][sC.vec[i].sol] + dist[sC.vec[i].sol][s.vec[0].sol] - dist[s.vec[j].sol][s.vec[0].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
            else
            {
                // cost to insert between i and i+1
                costInsertion = dist[s.vec[j].sol][sC.vec[i].sol] + dist[sC.vec[i].sol][s.vec[j+1].sol] - dist[s.vec[j].sol][s.vec[j+1].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
        }

        // insert the i-th point in the cheapest position
        s.vec.insert(s.vec.begin()+bestPosition+1,sC.vec[i]);
    }

    // last RK
    s.vec.push_back(aux);

    s.ofv = CalculateFitness(s);
}


/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem()
{
    //specific problem
    dist.clear();
    node.clear();
}

#endif