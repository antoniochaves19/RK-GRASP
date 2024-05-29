#ifndef _MultiStart_H
#define _MultiStart_H

// Global Variables
extern int debug;                           // 0 - run mode      		    1 - debug mode
extern int numDecoders;                     // number of decoders
extern int MAXTIME;                         // maximum runtime
extern float OPTIMAL;                       // optimal solution (if it is known)
extern struct timeval Tstart, Tend, Tbest;  // computational time (unix systems)  

extern int n;                               // size of cromossoms
extern TSol bestSolution;                   // best solution found in the A-BRKGA


void MultiStart()
{
    //Multi Start
    static TSol s;                              // solução corrente
    static TSol sMelhor;                        // melhor solução

    static int IterT;                           // iteracao corrente

    float currentTime = 0;              // computational time of the search process
    IterT = 0;

    // run the search process until stop criterion
    while(1)
    {
        // Create the initial solution with random keys 
        sMelhor.ofv = INFINITY;
       
		CreateInitialSolutions(s); 
		Decoder(s); 

		if (s.ofv < sMelhor.ofv)
		{
			sMelhor = s;
		}
        
        // update the best solution found
        updateBestSolution(sMelhor, bestSolution, Tbest);


        if (debug)
            printf("\n%d \tFO = %.2lf \t melhorFO = %.1lf", IterT, sMelhor.ofv, bestSolution.ofv);
			

        IterT++;

        // terminate the search process in MAXTIME
        gettimeofday(&Tend, NULL);
        currentTime = ((Tend.tv_sec  - Tstart.tv_sec) * 1000000u + Tend.tv_usec - Tstart.tv_usec) / 1.e6; 

        // stop criterium
        if (currentTime >= MAXTIME || floorf(bestSolution.ofv*10)/10 <= OPTIMAL){  
            break;
        }
    }

    // free memory with problem data
    FreeMemoryProblem();
}

#endif