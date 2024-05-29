// *******************************************************************
//      file with specific functions to solve a Problem
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Variables declared in main.cpp
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int n;                               // size of cromossoms

//----------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC -----------------------


//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------



//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------

/************************************************************************************
 Method: ReadData
 Description: read the input data
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
    

    fclose(arq);

    // define the size of the solution vector
    // n = 1;
}

/************************************************************************************
 Method: CalculateFitness
 Description: calculate the fitness of a chromossom s
*************************************************************************************/
double CalculateFitness(TSol s)
{    
    s.ofv = 0;

    return s.ofv;
}

/************************************************************************************
 Method: Decoders
 Description: users need to implement at least one decoder, DecK (K = [1,2,3,4,5])
*************************************************************************************/

void Dec1(TSol &s) 
{
    // create a solution of the problem

    s.ofv = CalculateFitness(s);
}

void Dec2(TSol &s){}
void Dec3(TSol &s){}
void Dec4(TSol &s){}
void Dec5(TSol &s){}

/************************************************************************************
 Method: Local Search Heuristics
 Description: users need to implement local searchs heuristics if he/she use the RVND,
 LSK (K = [1,2,3,4,5])
*************************************************************************************/

void LS1(TSol &s){}
void LS2(TSol &s){}
void LS3(TSol &s){}
void LS4(TSol &s){}
void LS5(TSol &s){}

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem(){}

#endif