#ifndef _PROBLEM_H
#define _PROBLEM_H

#define INFINITO 999999999

#include <math.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "../Main/Data.h"

// Variables declared in main.cpp
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int n;                               // size of cromossoms

//------ DEFINITION OF TYPES OF PROBLEM SPECIFIC --------

// struct with node informations
struct TNode								
{
	int id;
	double x;
	double y;
}; 

//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------

static std::vector <std::vector <double> > dist;	// matrix with Euclidean distance
static std::vector <TNode> node;					// vector of TSP nodes



//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------


/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
*************************************************************************************/
void ReadData(char nameTable[]);

/************************************************************************************
 Method: Decoder()
 Description: Convert a random key solution in a real problem solution
*************************************************************************************/
void Decoder(TSol &s);

/************************************************************************************
 Method: LocalSearch
 Description: RVND
*************************************************************************************/
void LocalSearch(TSol &s);


/************************************************************************************
 Method: CalculateFitness
 Description: calculate the fitness of a chromossom s
*************************************************************************************/
double CalculateFitness(TSol s);

/************************************************************************************
 Method: Dec1
 Description: sort decoder 
*************************************************************************************/
void Dec1(TSol &s);

/************************************************************************************
 Method: Dec2
 Description: 2-opt decoder 
*************************************************************************************/
void Dec2(TSol &s);

/************************************************************************************
 Method: Dec3
 Description: Cheapest Insertion decoder 
*************************************************************************************/
void Dec3(TSol &s);

/************************************************************************************
 Method: Dec4
 Description: k-Farthest Insertion decoder 
*************************************************************************************/
void Dec4(TSol &s);

/************************************************************************************
 Method: Dec5
 Description: k-Nearest Insertion decoder 
*************************************************************************************/
void Dec5(TSol &s);

/************************************************************************************
 Method: LS1
 Description: 2-Opt
*************************************************************************************/
void LS1(TSol &s);

/************************************************************************************
 Method: LS2
 Description: NodeInsertion
*************************************************************************************/
void LS2(TSol &s);

/************************************************************************************
 Method: LS3
 Description: NodeExchange
*************************************************************************************/
void LS3(TSol &s);

/************************************************************************************
 Method: LS4
 Description: OrOpt2
*************************************************************************************/
void LS4(TSol &s);

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem();

#endif