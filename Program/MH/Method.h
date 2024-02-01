#ifndef _Method_H
#define _Method_H

// gerador de números pseudoaleatório Mersenne Twister
static std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

/************************************************************************************
 Method: RANDOMICO
 Description: Generate a double random number between min and max
*************************************************************************************/
static double randomico(double min, double max)
{
    if (!debug)
        return std::uniform_real_distribution<double>(min, max)(rng);
    else
	    return ((double)(rand()%10000)/10000.0)*(max-min)+min;
}

/************************************************************************************
 Method: IRANDOMICO
 Description: Generate a int random number between min and max
*************************************************************************************/
static int irandomico(int min, int max)
{
	return (int)randomico(0,max-min+1.0) + min;
}

/************************************************************************************
 Method: sortByFitness
 Description: Sort TSol by objective function
*************************************************************************************/
static bool sortByFitness(const TSol &lhs, const TSol &rhs) { return lhs.ofv < rhs.ofv; }

/************************************************************************************
 Method: updateBestSolution
 Description: Update the best solution found during the search process
*************************************************************************************/
static void updateBestSolution(TSol s, TSol &bestSolution, struct timeval &Tbest)
{
    // save the best solution found 
    if (s.ofv < bestSolution.ofv)
    {
        bestSolution = s;
        gettimeofday(&Tbest, NULL);
    }
}

/************************************************************************************
 Method: CreateInitialSolutions
 Description: Create a initial random solution
*************************************************************************************/
static void CreateInitialSolutions(TSol &s)
{
	TVecSol aux;

    s.vec.clear();

	// create a random-key for each allelo (consider decoder type in the n-th random-key)
	for (int j = 0; j < n+1; j++)
	{
        aux.rk  = randomico(0,1);  // random value between [0,1[
        aux.sol = 0;
        s.vec.push_back(aux);
	}

    // flag to control the local search memory
    s.flag = 0;
}

/************************************************************************************
 Method: ShakeSolution
 Description: shake the current solution
*************************************************************************************/
static void ShakeSolution(TSol &s, float betaMin, float betaMax)
{
    int shaking_type = 0.0;
    int intensity = n*randomico(betaMin,betaMax);
    for(int k = 0; k < intensity; k++) {
        shaking_type = irandomico(1,4);
        int i = irandomico(0, n - 1);
        if(shaking_type == 1){
            // Invert value
            if (s.vec[i].rk > 0.0001)
                s.vec[i].rk = 1.0 - s.vec[i].rk;
            else
                s.vec[i].rk = 0.9999;
        }
        else 
        if (shaking_type == 2){
            // Swap two random positions
            int j = irandomico(0, n - 1);
            double temp = s.vec[i].rk;
            s.vec[i].rk = s.vec[j].rk;
            s.vec[j].rk = temp;
        }
        else
        if(shaking_type == 3){
            // Change to random value
            s.vec[i].rk = randomico(0,1);
        }
        i = irandomico(0, n - 2);
        if(shaking_type == 4){
            // Swap with neighbor
            double temp = s.vec[i].rk;
            s.vec[i].rk = s.vec[i+1].rk;
            s.vec[i+1].rk = temp;
        }
    }
}

#endif