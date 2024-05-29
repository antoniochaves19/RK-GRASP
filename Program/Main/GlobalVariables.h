#ifndef _GLOBALVARIABLES_H
#define _GLOBALVARIABLES_H

// /************************************************************************************
// 								GLOBAL VARIABLES
// *************************************************************************************/

#define INFINITO 999999999

// Input File
char instance[100];                         // name of instance
int debug = 1;                              // 0 - run mode      		    1 - debug mode
int numDecoders = 1;                        // number of decoders
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

std::vector <TSol> pool;                    // pool of best solutions with diversity
int sizePool = 20;                          // number of solution in the pool

#endif