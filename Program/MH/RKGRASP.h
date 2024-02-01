#ifndef _CGRASP_H
#define _CGRASP_H

// Variables declared in the main
extern int debug;                           // 0 - run mode      		    1 - debug mode
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int MAXTIME;                         // maximum runtime
extern float OPTIMAL;                       // optimal solution (if it is known)
extern struct timeval Tstart, Tend, Tbest;  // computational time (unix systems)  

extern int n;                               // size of random-key vector
extern TSol bestSolution;                   // best solution found for any method


/************************************************************************************
			              Declaration of the Functions
*************************************************************************************/

/************************************************************************************
 Method: ConstrutiveGreedyRandomized
 Description: It takes as input a solution vector s. Initially, the procedure allows 
 all coordinates of s to change (i.e. they are called unfixed). In turn, a discrete 
 line search is performed in each unfixed coordinate direction i of s with the other 
 n − 1 coordinates of s held at their current values.
*************************************************************************************/
static void ConstrutiveGreedyRandomized(TSol &s, float h, float alpha);

/************************************************************************************
 Method: LineSearch
 Description: Consider the line search in direction ei , where vector ei has zeros 
 in all components except the i-th, where it has value one. The objective function 
 is evaluated at points x + k·h·ei for k = 0,1,−1,2,−2,... such that li≤xi+k·h≤ui.
 Let k∗ the value of k that minimizes f(x+k·h·ei) subject to li ≤ xi +k·h ≤ ui
*************************************************************************************/
static void LineSearch(TSol s, float h, int i, double &bestZ, double &bestF);

/************************************************************************************
 Method: GridSearch
 Description: Generates a neighborhood and determines at which points in the neighbor-
 hood, if any, the objective function improves. If an improving point is found, it is 
 made the current point and the local search continues from the new solution.
*************************************************************************************/
static void GridSearch(TSol &x, float h);

/************************************************************************************
 Method: RandomSelectElement
 Description: Consider the set of points in S that are integer steps (of size h) away 
 from x. Define the projection of the points in S onto the hyper-sphere centered at x  
 of radius h. The h-neighborhood of the point x is defined as the set of points in Bh.
 The algorithm randomly selects points in Bh, one at a time. 
*************************************************************************************/
static TSol RandomlySelectElement(TSol x, float h);

/************************************************************************************
 Method: NelderMead
 Description: The Nelder–Mead method is a numerical method used to find the minimum 
 of an objective function in a multidimensional space.
*************************************************************************************/
static TSol NelderMead(TSol x1, TSol x2, TSol x3, float h);

/************************************************************************************
 Method: SwapRK
 Description: swap two rk of positions
*************************************************************************************/
static void SwapRK(TSol &s, float h);


/************************************************************************************
			                  Implementation
*************************************************************************************/

extern void RKGRASP(char nameTable[256])
{
    // C-GRASP parameters
    static TSol s;                              // solution
    static TSol sCurrent;                       // current solution
    static TSol sBest;                          // best solution of C-GRASP

    static double alpha;                        // greedy rate

    static float h;                             // grid dense
    static float hs = 0.12500;                  // start grid dense
    static float he = 0.00098;                  // end grid dense

    // free memory with problem data
    FreeMemoryProblem();

    //read data of the instance
    ReadData(nameTable);

    // computational time of the search process
    float currentTime = 0;       

    // run the search process until stop criterion (maxTime)
    int iter = 0;

    // create an initial solution
    CreateInitialSolutions(s); 
    Decoder(s);

    sBest = s;
    updateBestSolution(sBest, bestSolution, Tbest);

    while(1)
    {
        h = hs;

        while (h >= he)
        {
            iter++;

            // current solution
            sCurrent = s;

            // generate a random value for alpha
            alpha = randomico(0.0,0.99);  

            // choose a decoder randomly
            s.vec[n].rk = randomico(0,1);

            // construct a greedy randomized solution
            ConstrutiveGreedyRandomized(s, h, alpha);
            double ofvConstrutive = s.ofv;

            if (debug) printf("\nIter = %d \t h = %.4lf \t Constructive (%.3lf) = %8.2lf", iter, h, alpha, floorf(s.ofv*10)/10);
            
            // apply local search in current solution
            int ls = 4; //change the ls to use different local searchs

            if (ls == 1) GridSearch(s, h);
            if (ls == 2) s = NelderMead(sBest, sCurrent, s, h);
            if (ls == 3) SwapRK(s, h);
            if (ls == 4) LocalSearch(s); //available in the Problem.h

            if (debug) printf("\t Local search = %8.2lf (%8.2lf)", floorf(s.ofv*10)/10, ofvConstrutive - s.ofv);

            if (s.ofv < sBest.ofv){
                // update the best solution found by C-GRASP
                sBest = s;
                
                // update the best global solution found
                updateBestSolution(sBest, bestSolution, Tbest);

                if (debug) printf("**");
            }
            else{
                // make grid more dense
                h = h/2;
            }

            if (debug) printf("\t Best solution = %8.2lf", floorf(bestSolution.ofv*10)/10);

            // terminate the search process in MAXTIME
            gettimeofday(&Tend, NULL);
            currentTime = ((Tend.tv_sec  - Tstart.tv_sec) * 1000000u + Tend.tv_usec - Tstart.tv_usec) / 1.e6; 

            // stop criterium
            if (currentTime >= MAXTIME || floorf(bestSolution.ofv*10)/10 <= OPTIMAL){  
                break;
            }
        } 

        // create a new random seed
        if (!debug)
            srand(time(NULL)); 

        // restart the current solution 
        CreateInitialSolutions(s); 

        // decode the new solution
        Decoder(s);

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

static void ConstrutiveGreedyRandomized(TSol &s, float h, float alpha)
{
    std::vector<int> UnFixed(n);                // stores random-keys not yet fixed
    std::vector<int> chosenRK;                  // stores the random-keys that will be searched in the line search
    std::vector<int> RCL;                       // stores the best candidate solutions
    std::vector<double> z(n);                   // stores the best value of random-key i
    std::vector<double> g(n,INFINITY);          // stores the fo value with the random-key z_i

    double min, max;
    double betaMin = 0.5, 
           betaMax = 0.8;

    // initialize the points on the chromosome that can be changed
    for (int i = 0; i < n; i++) {
        UnFixed[i] = i;
    }

    // build a solution through perturbations in the current solution and choosing one of the best
    double intensity = randomico(betaMin, betaMax);
    
    // while (!UnFixed.empty())
    for (int j=0; j<n*intensity; j++)
    {
        // create a list of candidate solutions perturbing a (not yet 'fixed') rk of the current solution
        min = INFINITY;
        max = -INFINITY;

        // choose the subset of random keys that will be searched
        chosenRK.clear();
        int kMax = UnFixed.size() * 0.5;
        if (kMax < 2) 
            kMax = UnFixed.size();

        chosenRK = UnFixed;
        std::random_shuffle(chosenRK.begin(), chosenRK.end());
        chosenRK.resize(kMax);

        // line search
        z.clear();
        z.resize(n);
        g.clear();
        g.resize(n,INFINITY);

        #pragma omp parallel for num_threads(MAX_THREADS) 
        for (int k=0; k<kMax; k++) 
        {
            int i = chosenRK[k];
            TSol sAux = s;
            LineSearch(sAux, h, i, z[i], g[i]);
        }

        for(int i=0; i<n; i++)
        {
            if (min > g[i])
                min = g[i];

            if (max < g[i] && g[i] != INFINITY)
                max = g[i];
        }

        // create RCL
        RCL.clear();
        double threshold = min + alpha * (max - min);

        for (int i=0; i<n; i++)
        { 
            if (g[i] <= threshold)
            {
                RCL.push_back(i);
            }
        }

        // randomly select one of the best candidates to continue construction
        if (!RCL.empty()){
            int x = irandomico(0, (int)(RCL.size())-1);
            int kCurrent = RCL[x];  // index of the random key that will be fixed
            
            // update the current solution
            s.vec[kCurrent].rk = z[kCurrent];
            s.ofv = g[kCurrent];

            // remove rk k from the UnFixed set
            for (int i=0; i<(int)UnFixed.size(); i++)
            {
                if (UnFixed[i] == kCurrent) 
                {
                    UnFixed.erase(UnFixed.begin()+i);
                    break;
                }
            }
        }
    }

    // update the solution found in the constructive phase
    Decoder(s);
}

static void LineSearch(TSol s, float h, int i, double &bestZ, double &bestF)
{
    std::vector<double> rk;

    // find the best solution in the line
    bestZ = 0;
    bestF = INFINITY;

    // generate k as the possible random keys of gene i, k = rk_i * tau | tau = {0, 1, -1, 2, -2, ...}
    double tau = 0;
    rk.push_back(s.vec[i].rk + tau * h);
    for (int j=0; j<(int)(1.0/h)+1; j+=2){
        tau++;

        if (((s.vec[i].rk + tau * h) >= 0) && ((s.vec[i].rk + tau * h) < 1))
            rk.push_back(s.vec[i].rk + tau * h);

        if (((s.vec[i].rk + (-1*tau) * h) >= 0) && ((s.vec[i].rk + (-1*tau) * h) < 1))
            rk.push_back(s.vec[i].rk + (-1*tau) * h);
    }

    // (sample greedy) This method is similar to the greedy algorithm, but instead of selecting the best among all possible options, 
    // it only considers q < m possible insertions (chosen uniformly at random) in each iteration. The most profitable among those is selected. 
    int q = ceil(log2((int)(1.0/h))) + 1; //1; //

    if (q > (int)rk.size())
        q = rk.size();

    // choose a subset with q rks to calculate the decoder
    std::random_shuffle (rk.begin(), rk.end());

    // calculate the quality of the solution s with rk j
    for (int j=0; j<q; j++)
    {  
        s.vec[i].rk = rk[j];   
        Decoder(s);

        if (s.ofv < bestF){
            bestZ = s.vec[i].rk;
            bestF = s.ofv;
        }
    }

    rk.clear();
}

static void GridSearch(TSol &x, float h)
{
    int numGridPoints = n * floor(1.0/h);

    // initialize the best solution with the current solution x (rk and ofv)
    TSol xBest = x;

    #pragma omp parallel for num_threads(MAX_THREADS) 
    for (int i=0; i <= numGridPoints; i++ ) 
    {   
        // generate a solution neighboring xBest
        TSol xLocal = RandomlySelectElement(xBest, h);
        
        // decoder
        Decoder(xLocal);

        #pragma omp critical
        {
            if (xLocal.ofv < xBest.ofv){
                xBest = xLocal;
            }
        }
    }

    // return the best solution
    x = xBest;
}

static TSol RandomlySelectElement(TSol y, float h)
{
    double norm;

    // direction vector
    std::vector<double> Sh;      
    Sh.resize(n);

    // calculate S_h and norm
    norm = 0.0;
    for (int i=0; i<n; i++)
    {
        Sh[i] = y.vec[i].rk + irandomico(-100, 100) * h; 

        norm += pow(Sh[i] - y.vec[i].rk,2);
    }
    norm = sqrt(norm);

    // generate B_h (solution neighboring y)
    for (int i=0; i<n; i++)
    {
        y.vec[i].rk = y.vec[i].rk + h * ((Sh[i] - y.vec[i].rk) / norm);

        if (y.vec[i].rk < 0)
            y.vec[i].rk = randomico(0,1);
        else
        if (y.vec[i].rk >= 1.0)
            y.vec[i].rk = randomico(0,1);
    }

    Sh.clear();
    return y;
}

static TSol NelderMead(TSol x1, TSol x2, TSol x3, float h)
{
    const float ALPHA = 1.0;
    const float GAMMA = 2.0;
    const float PHI = 0.5;
    const float SIGMA = 0.5;

    // internal points
    TSol point_r;
    TSol point_e;
    TSol point_c;
    TSol centroid;

    int iter_count = 0;
    int eval_count = 0;
 
    // initial simplex has size nDims + 1
    int nDims = 2;
    std::vector<TSol> simplex(nDims+1);
    simplex[0] = x1;
    simplex[1] = x2;
    simplex[2] = x3;

    // perturb the current solution
    ShakeSolution(simplex[1],0.10,0.20);
    Decoder(simplex[1]);

    // sort points in the simplex so that simplex[0] is the point having
    // minimum fx and simplex[n] is the one having the maximum fx
    sort(simplex.begin(), simplex.end(), sortByFitness);
    
    // compute the simplex centroid
    centroid.vec.resize(n+1);
    for (int j = 0; j < nDims; j++) 
    {
        for (int i = 0; i < n; i++) 
        {
            centroid.vec[i].rk += simplex[j].vec[i].rk / (nDims);
        }
    }

    iter_count++; 

    int numIter = n * floor(1.0/h);
    // continue minimization until stop conditions are met
    while (eval_count < numIter) 
    {
        int shrink = 0;

        // reflection point (r)
        point_r.vec.resize(n+1);
        for (int i = 0; i < n; i++) 
        {
            point_r.vec[i].rk = (centroid.vec[i].rk + ALPHA * (centroid.vec[i].rk - simplex[nDims].vec[i].rk))/(ALPHA + 1);

            if ((point_r.vec[i].rk < 0) || (point_r.vec[i].rk >= 1))
                point_r.vec[i].rk = randomico(0,1);
        }
        Decoder(point_r);
        eval_count++;

        // point_r is better than the bestSolution
        if (point_r.ofv < simplex[0].ofv) 
        {
            // expansion point (e)
            point_e.vec.resize(n+1);
            for (int i = 0; i < n; i++) 
            {
                point_e.vec[i].rk = (centroid.vec[i].rk + GAMMA * (point_r.vec[i].rk - centroid.vec[i].rk))/(GAMMA + 1);
                
                if ((point_e.vec[i].rk < 0) || (point_e.vec[i].rk >= 1))
                    point_e.vec[i].rk = randomico(0,1);
            }
            Decoder(point_e);
            eval_count++;

            if (point_e.ofv < point_r.ofv) 
            {
                // expand
                // if (debug) printf("\nexpand          ");
                simplex[nDims] = point_e;
            } 
            else 
            {
                // reflect
                // if (debug) printf("\nreflect         ");
                simplex[nDims] = point_r;
            }
        } 
        // point_r is NOT better than the bestSolution
        else 
        {    
            // point_r is better than the second best solution
            if (point_r.ofv < simplex[nDims-1].ofv) 
            {
                // reflect
                // if (debug) printf("\nreflect         ");
                simplex[nDims] = point_r;
            } 
            else 
            {
                // point_r is better than the worst solution
                if (point_r.ofv < simplex[nDims].ofv) 
                {
                    // contraction point (c)
                    point_c.vec.resize(n+1);
                    for (int i = 0; i < n; i++) 
                    {
                        point_c.vec[i].rk = (centroid.vec[i].rk + PHI * (point_r.vec[i].rk - centroid.vec[i].rk))/(PHI + 1);

                        if ((point_c.vec[i].rk < 0) || (point_c.vec[i].rk >= 1))
                            point_c.vec[i].rk = randomico(0,1);
                    }
                    Decoder(point_c);
                    eval_count++;

                    if (point_c.ofv < point_r.ofv) 
                    {
                        // contract outside
                        // if (debug) printf("\ncontract out    ");
                        simplex[nDims] = point_c;
                    } 
                    else 
                    {
                        // shrink
                        // if (debug) printf("\nshrink1         ");
                        shrink = 1;
                    }
                } 
                else 
                {
                    // contraction point (c)
                    point_c.vec.resize(n+1);
                    for (int i = 0; i < n; i++) 
                    {
                        point_c.vec[i].rk = (centroid.vec[i].rk + PHI * (simplex[nDims].vec[i].rk - centroid.vec[i].rk))/(PHI + 1);

                        if ((point_c.vec[i].rk < 0) || (point_c.vec[i].rk >= 1))
                            point_c.vec[i].rk = randomico(0,1);
                    }
                    Decoder(point_c);
                    eval_count++;

                    if (point_c.ofv < simplex[nDims].ofv) 
                    {
                        // contract inside
                        // if (debug) printf("\ncontract in     ");
                        simplex[nDims] = point_c;
                    } 
                    else {
                        // shrink
                        // if (debug) printf("\nshrink2          ");
                        shrink = 1;
                    }
                }
            }
        }
        if (shrink) {
            for (int j = 1; j < nDims+1; j++) 
            {
                for (int i = 0; i < n; i++) 
                {
                    simplex[j].vec[i].rk = (simplex[0].vec[j].rk + SIGMA * (simplex[j].vec[i].rk - simplex[0].vec[i].rk))/(SIGMA + 1);

                    if ((simplex[j].vec[i].rk < 0) || (simplex[j].vec[i].rk >= 1))
                            simplex[j].vec[i].rk = randomico(0,1);
                }
                Decoder(simplex[j]);
                eval_count++;
            }

            // sort
            sort(simplex.begin(), simplex.end(), sortByFitness);
        }
        else
        {
            for (int i = nDims - 1; i >= 0 && simplex[i + 1].ofv < simplex[i].ofv; i--) 
            {
                TSol temp;
                temp = simplex[i+1];
                simplex[i+1] = simplex[i];
                simplex[i] = temp;
            }
        }

        // compute the simplex centroid
        centroid.vec.resize(n+1);
        for (int j = 0; j < nDims; j++) 
        {
            for (int i = 0; i < n; i++) 
            {
                centroid.vec[i].rk += simplex[j].vec[i].rk / (nDims);
            }
        }

        iter_count++;
    }

    return simplex[0];
}

static void SwapRK(TSol &s, float h)
{
    // swap rk's
    std::vector<double> rk1;
    std::vector<double> rk2;

    for (int i=0; i<n-1; i++)
    {
        for (int j = i+1; j<n; j++)
        {
            rk1.push_back(s.vec[i].rk);
            rk2.push_back(s.vec[j].rk);
        }
    }

    #pragma omp parallel for num_threads(MAX_THREADS)
    for (int k=0; k<(int)rk1.size(); k++)
    {
        TSol sViz = s;

        // exchange a rk i for a rk j
        double temp = sViz.vec[rk1[k]].rk;
        sViz.vec[rk1[k]].rk = sViz.vec[rk2[k]].rk;
        sViz.vec[rk2[k]].rk = temp;

        Decoder(sViz);
        
        // update the current solution if it has improved and continue the search from sViz
        #pragma omp critical
        {
            if (sViz.ofv < s.ofv){
                s = sViz;
            }
        }
    }  
}

#endif