#ifndef _GLOBALVARIABLES_H
#define _GLOBALVARIABLES_H

// /************************************************************************************
// 								GLOBAL VARIABLES
// *************************************************************************************/

// Input File
char instance[100];                         // name of instance
int debug = 1;                              // 0 - run mode      		    1 - debug mode
int numDecoders = 1;                        // number of decoders
int numLS = 1;                              // 0 - without local search     > k - number of local search heuristics
int MAXTIME = 1;                            // maximum runtime
int MAXRUNS =  1;                           // maximum number of runs of the method
unsigned MAX_THREADS = 40;            		// number of threads
float OPTIMAL = 0;                          // optimal solution (if it is known) or a lower bound

int n;                                      // size of the solution

// Run
struct timeval Tstart, Tend, Tbest;         // computational time (unix systems)
char nameTable[256];                        // name of instance
char nameMH[256];                           // name of the metaheuristic
TSol bestSolution;                          // best solution found

#endif