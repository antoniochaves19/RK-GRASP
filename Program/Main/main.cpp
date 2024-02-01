#include <sys/time.h>
#include <math.h>
#include <cstring>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cstring>
#include <string.h>
#include <sys/time.h>
#include <utility>  // pair
#include <numeric>  // iota
#include <map>
#include <limits>
#include <random>
#include <chrono>
#include <iomanip> //graph
#include <sstream> //graph
#include <fstream> //graph

#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

#include "Data.h" 
#include "GlobalVariables.h" 
#include "Output.h"
#include "../Problem/Problem.h"
#include "../MH/Method.h"
#include "../MH/RKGRASP.h"
#include "../MH/BRKGA_QL.h"


// /************************************************************************************
// 								MAIN FUNCTION AREA
// *************************************************************************************/
int main(int argc, char *argv[ ])
{ 
    // read scenario name file 
    char nameScenario[256];
    strncpy(nameScenario,argv[1],255);
    int method = atoi(argv[2]);  

    // file with test instances and input data
	FILE *arqProblems; 
    arqProblems = fopen (nameScenario, "r"); 

    if (arqProblems == NULL){
        printf("\nERROR: File %s not found\n", nameScenario);
        exit(1);
    }

    //read first line of testScenario file
    if (fgets(nameTable, sizeof(nameTable), arqProblems) == NULL) {
        printf("\nERROR: File %s not found\n", nameTable);
        exit(1);
    }
    
    // best solution that is saved in out file
    TSol sBest;
    sBest.flag = 0;
    sBest.label = 0;
    sBest.ofv = 0;
    sBest.similar = 0;
    sBest.promising = 0;
    sBest.vec.clear();

	// run the RKO-MH for all test instances
	while (!feof(arqProblems))
	{
		// read the name of instances, debug mode, local search module, maximum time, maximum number of runs, maximum number of threads
        if (fscanf(arqProblems,"%s %d %d %d %d %d %d %f", nameTable, &debug, &numDecoders, &numLS, &MAXTIME, &MAXRUNS, &MAX_THREADS, &OPTIMAL) == 0) {
            printf("\nERROR: File %s not found\n", nameTable);
            exit(1);
        }
        strcpy(instance,nameTable);
        
        double foBest = INFINITY,
               foAverage = 0;

        float timeBest = 0,
              timeTotal = 0;

        std::vector <double> ofvs;
        ofvs.clear();

        // best solutions found in MAXRUNS
        sBest.ofv = INFINITY;

		// run RKO MaxRuns for each instance
        printf("\n\nInstance: %s [Threads %d]\nRun: ", instance, MAX_THREADS);
        for (int j=0; j<MAXRUNS; j++)
        {
            // fixed seed
            if (debug == 1)
                srand(j+1); 
            else
                srand(time(NULL)); 

            printf("%d ", j+1);
            
            gettimeofday(&Tstart, NULL);
            gettimeofday(&Tend, NULL);
            gettimeofday(&Tbest, NULL);

            // best solution found in this run
            bestSolution.ofv = INFINITY;

            // execute the RKO method
            switch (method) 
            {
            case 1:
                strcpy(nameMH,"BRKGA-QL");
                BRKGA_QL(nameTable);
                break;

            case 2:
                strcpy(nameMH,"RK-GRASP");
                RKGRASP(nameTable);
                break;
            
            default:
                break;
            } 

            gettimeofday(&Tend, NULL);

            // store the best solution found in MAXRUNS
            if (bestSolution.ofv < sBest.ofv)
                sBest = bestSolution;

            // calculate best and average results
            if (bestSolution.ofv < foBest)
                foBest = bestSolution.ofv;

            foAverage += bestSolution.ofv;

            // fitness of each solution found in the runs
            ofvs.push_back(bestSolution.ofv);

            timeBest += ((Tbest.tv_sec  - Tstart.tv_sec) * 1000000u + Tbest.tv_usec - Tstart.tv_usec) / 1.e6;
            timeTotal += ((Tend.tv_sec  - Tstart.tv_sec) * 1000000u + Tend.tv_usec - Tstart.tv_usec) / 1.e6; 
        }

        // create a .csv file with average results
        foAverage = foAverage / MAXRUNS;
        timeBest = timeBest / MAXRUNS;
        timeTotal = timeTotal / MAXRUNS;

        if (!debug)
        {
        	WriteSolution(nameMH, sBest, n, timeBest, timeTotal, instance);
        	WriteResults(nameMH, foBest, foAverage, ofvs, timeBest, timeTotal, instance);
        }
        else
        {
            WriteSolutionScreen(nameMH, sBest, n, timeBest, timeTotal, instance);
        }
    }

    fclose(arqProblems);
    return 0;
}
