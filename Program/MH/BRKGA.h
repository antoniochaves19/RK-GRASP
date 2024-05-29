#ifndef _BRKGA_H
#define _BRKGA_H

// Global Variables
extern int debug;                           // 0 - run mode      		    1 - debug mode
extern int numDecoders;                     // number of decoders
extern int MAXTIME;                         // maximum runtime
extern float OPTIMAL;                       // optimal solution (if it is known)
extern struct timeval Tstart, Tend, Tbest;  // computational time (unix systems)  
extern unsigned MAX_THREADS;                // number of threads

extern int n;                               // size of cromossoms
extern TSol bestSolution;                   // best solution found in the MH

/************************************************************************************
 Method: PARAMETRICUNIFORMCROSSOVER
 Description: create a new offspring with parametric uniform crossover
*************************************************************************************/
static TSol ParametricUniformCrossover(int elitesize, int popSize, double rhoe, std::vector <TSol> &Pop);

/************************************************************************************
			                  GENERAL FUNCTIONS
*************************************************************************************/
static void BRKGA(int method, int control)
{
    static int p;          	                        // size of population
    static double pe;              	                // fraction of population to be the elite-set
    static double pm;          	                    // fraction of population to be replaced by mutants
    static double rhoe;             	            // probability that offspring inherit an allele from elite parent

    static std::vector <TSol> Pop;                  // current population
    static std::vector <TSol> PopInter;             // intermediary population

    // offline control
    if (control == 0){
        p    = 1000;
        pe   = 0.20;                                                         
        pm   = 0.20;                                                    
        rhoe = 0.75;
    }
    
    // initialize population
    Pop.clear(); 
    PopInter.clear(); 
    Pop.resize(p);
    PopInter.resize(p);

    // Create the initial chromosomes with random keys
    for (int i=0; i<p; i++)
    {
        TSol ind;
        CreateInitialSolutions(ind); 
        Decoder(ind);
        Pop[i] = PopInter[i] = ind;
    }
    // Pop[0] = pool[0];
    
    // sort population in increase order of fitness
    sort(Pop.begin(), Pop.end(), sortByFitness);

    // save the best solution found
    updateBestSolution(Pop[0], bestSolution, Tbest);
    
    // useful local variables
    int numGenerations = 0;             // number of generations
    int bestGeneration = 0;             // generation in which found the best solution
    double bestFitness = Pop[0].ofv;    // best fitness found in past generation
    double bestOFV = 0;                 // best fitness found in the crossover
    double delta = 0;                   // difference between the best offspring and the best solution found
    float currentTime = 0;              // computational time of the search process
    
    // run the evolutionary process until stop criterion
    while(1)
    {
    	// number of generations
        numGenerations++;

        // The 'Pe' best chromosomes are maintained, so we just copy these into PopInter:
        #pragma omp parallel for num_threads(MAX_THREADS)
        for (int i=0; i<(int)(p*pe); i++){
            // copy the chromosome for next generation
            PopInter[i] = Pop[i]; 
        }  

        // We'll mate 'p - pe - pm' pairs; initially, i = pe, so we need to iterate until i < p - pm:
        bestOFV = INFINITY;
        #pragma omp parallel for num_threads(MAX_THREADS)
        for (int i = (int)(p*pe); i < p - (int)(p*pm); i++){
            // Parametric uniform crossover with mutation
            PopInter[i] = ParametricUniformCrossover((int)(p*pe), p, rhoe, Pop);
 
            // Calculate the fitness of new chromosomes
            Decoder(PopInter[i]); 

            if (PopInter[i].ofv < bestOFV){
                bestOFV = PopInter[i].ofv;
            }
        }
        
        // We'll introduce 'pm' mutants:
        #pragma omp parallel for num_threads(MAX_THREADS)
        for (int i = p - (int)(p*pm) - (int)(p*pe); i < p; i++){
            CreateInitialSolutions(PopInter[i]);
            Decoder(PopInter[i]);
        }  
                
        // Update the current population
        Pop = PopInter;   

        // Sort population in increase order of fitness
        sort(Pop.begin(), Pop.end(), sortByFitness);
        updateBestSolution(Pop[0], bestSolution, Tbest);

        // We improve the best fitness in the current population 
        delta = 0;
        if (Pop[0].ofv < bestFitness){
            delta = Pop[0].ofv - bestFitness;
            bestFitness = Pop[0].ofv;
            bestGeneration = numGenerations;
        }

        // print screen 
        if (debug){
            printf("\nGeneration: %4d \t %.1lf (%.2lf)  \t %.1lf | \t [%4d, %.2lf, %.2lf, %.2lf] \t",
                        numGenerations, bestSolution.ofv, bestSolution.vec[n].rk, bestFitness,
                        p, pe, pm, rhoe);
        }

        // ******************************* SHAKING *****************************
        if ((numGenerations - bestGeneration) > n) {
            
            if (debug) 
                printf("\n\nShaking elite and reset non-elite...\n\n");
            else
                srand(time(NULL)); 

            // reset the number of generations without improvement
            bestGeneration = numGenerations;

            for (int i=0; i<p; i++){

                // Shake the elite set
                if (i < (int)(p*pe)){
                    ShakeSolution(Pop[i], rhoe, rhoe);
                }

                // reset the non-elite chromosomes
                else{
                    CreateInitialSolutions(Pop[i]);
                }
                
                // decode the new solution
                Decoder(Pop[i]);
            }

            sort(Pop.begin(), Pop.end(), sortByFitness);
            updateBestSolution(Pop[0], bestSolution, Tbest);
            bestFitness = Pop[0].ofv;
        }
        // *********************************************************************

        // terminate the evolutionary process in MAXTIME
        gettimeofday(&Tend, NULL);
        currentTime = ((Tend.tv_sec  - Tstart.tv_sec) * 1000000u + Tend.tv_usec - Tstart.tv_usec) / 1.e6; 
        
        // stop criterium
        if (currentTime >= MAXTIME || floorf(bestSolution.ofv*10)/10 <= OPTIMAL){  
            break;
        }
    }

    // free memory with problem data
    FreeMemoryProblem();

    // free memory of BRKGA components
    Pop.clear();
    PopInter.clear();
}

static TSol ParametricUniformCrossover(int eliteSize, int popSize, double rhoe, std::vector <TSol> &Pop)
{	
	TSol s;

    int eliteParent = irandomico(0, eliteSize - 1);                 // one chromosome from elite set
    int nonEliteParent = irandomico(eliteSize, popSize-1);          // one chromosome from nonelite  population

    // best fit parent is parent elite
    if (Pop[eliteParent].ofv > Pop[nonEliteParent].ofv){
        int temp = eliteParent;
        eliteParent = nonEliteParent;
        nonEliteParent = temp;
    }

	// create a new offspring
	s.vec.resize(n+1);

    // Mate: including decoder gene in the n-th rk 
    for(int j = 0; j < n+1; j++)
    {
        //copy alelos of top chromossom of the new generation
        if (randomico(0,1) < rhoe){
            s.vec[j].rk = Pop[eliteParent].vec[j].rk;
        }
        else{
            s.vec[j].rk = Pop[nonEliteParent].vec[j].rk;
        }
    }

    // set the flag of local search as zero
    s.flag = 0;

    return s;
}


#endif