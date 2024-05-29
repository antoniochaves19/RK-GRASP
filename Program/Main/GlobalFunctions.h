#ifndef _GlobalFunctions_H
#define _GlobalFunctions_H

/************************************************************************************
 Method: RANDOMICO
 Description: Generate a double random number between min and max
*************************************************************************************/
static double randomico(double min, double max);

/************************************************************************************
 Method: IRANDOMICO
 Description: Generate a int random number between min and max
*************************************************************************************/
static int irandomico(int min, int max);

/************************************************************************************
 Method: sortByFitness
 Description: Sort TSol by objective function
*************************************************************************************/
static bool sortByFitness(const TSol &lhs, const TSol &rhs);

/************************************************************************************
 Method: updateBestSolution
 Description: Update the best solution found during the search process
*************************************************************************************/
static void updateBestSolution(TSol s, TSol &bestSolution, struct timeval &Tbest);

/************************************************************************************
 Method: CreateInitialSolutions
 Description: Create a initial random solution
*************************************************************************************/
static void CreateInitialSolutions(TSol &s);

/************************************************************************************
 Method: ShakeSolution
 Description: shake the current solution
*************************************************************************************/
static void ShakeSolution(TSol &s, float betaMin, float betaMax);

/************************************************************************************
 Method: CretePoolSolutions
 Description: create a pool of solutions with different solutions
*************************************************************************************/
static void CretePoolSolutions();

/************************************************************************************
 Method: UpdatePoolSolutions
 Description: update the pool with different solutions
*************************************************************************************/
static void UpdatePoolSolutions(TSol s);

/************************************************************************************
 Method: RandomSelectElement
 Description: Consider the set of points in S that are integer steps (of size h) away 
 from x. Define the projection of the points in S onto the hyper-sphere centered at x  
 of radius h. The h-neighborhood of the point x is defined as the set of points in Bh.
 The algorithm randomly selects points in Bh, one at a time. 
*************************************************************************************/
static TSol hNeighborhood(TSol y, float h); //, float theta

/************************************************************************************
 Method: GridSearch
 Description: Generates a neighborhood and determines at which points in the neighbor-
 hood, if any, the objective function improves. If an improving point is found, it is 
 made the current point and the local search continues from the new solution.
*************************************************************************************/
static void GridSearch(TSol &x, float h); //, float theta

/************************************************************************************
 Method: NelderMeadSearch
 Description: The Nelder–Mead method is a numerical method used to find the minimum 
 of an objective function in a multidimensional space. It is a direct search method 
 based on function comparison.
*************************************************************************************/
static void NelderMeadSearch(TSol &x1);

/************************************************************************************
 Method: Random Variable Neighborhood Descent
 Description: RVND
*************************************************************************************/
static void RVND(TSol &s, float h);

/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
*************************************************************************************/
static void ReadData(char nameTable[]);

/************************************************************************************
 Method: Decoder()
 Description: Convert a random key solution in a real problem solution
*************************************************************************************/
static void Decoder(TSol &s);

/************************************************************************************
 Method: Dec1
 Description:define by users
*************************************************************************************/
static void Dec1(TSol &s);

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem();

/************************************************************************************
 Method: sortByRk
 Description: Sort TSol by random-keys
*************************************************************************************/
bool sortByRk(const TVecSol &lhs, const TVecSol &rhs);

#endif