#ifndef _Method_H
#define _Method_H

// pseudo-random number generator Mersenne Twister
static std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

/************************************************************************************
 Method: sortByRk
 Description: Sort TSol by random-keys
*************************************************************************************/
bool sortByRk(const TVecSol &lhs, const TVecSol &rhs) { return lhs.rk < rhs.rk; }

/************************************************************************************
 Method: sortByFitness
 Description: Sort TSol by objective function
*************************************************************************************/
static bool sortByFitness(const TSol &lhs, const TSol &rhs) { return lhs.ofv < rhs.ofv; }

/************************************************************************************
 Method: RANDOMICO
 Description: Generate a double random number between min and max
*************************************************************************************/
static double randomico(double min, double max)
{
    return std::uniform_real_distribution<double>(min, max)(rng);
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
    int intensity = (int)(n * randomico(betaMin, betaMax)) + 1;
    if (intensity < 1) intensity = 1;
    for(int k = 0; k < intensity; k++) {
        shaking_type = irandomico(1,4);
        // shaking_type = 1;
        int i = irandomico(0, n - 1);
        if(shaking_type == 1){
            // Change to random value
            s.vec[i].rk = randomico(0,1);
        }
        else
        if(shaking_type == 2){
            // Invert value
            if (s.vec[i].rk > 0.0001)
                s.vec[i].rk = 1.0 - s.vec[i].rk;
            else
                s.vec[i].rk = 0.9999;
        }
        else 
        if (shaking_type == 3){
            // Swap two random positions
            int j = irandomico(0, n - 1);
            double temp = s.vec[i].rk;
            s.vec[i].rk = s.vec[j].rk;
            s.vec[j].rk = temp;
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

/************************************************************************************
 Method: Decoder
 Description: Convert a random key solution in a real problem solution
*************************************************************************************/
void Decoder(TSol &s)
{
    // copy the random-key sequence of current solution 
    TSol temp = s;

    // define decoder function based in the random-key of position n+1
    // int dec = floor(s.vec[n].rk*numDecoders)+1;
    int dec = 1;
    switch (dec)
    {
        case 1: 
            Dec1(s);
            break;

        default:
            break;
    }

    // return initial random-key sequence and maintain the solution sequence
    for (int i=0; i<n; i++){
        s.vec[i].rk = temp.vec[i].rk;
    }
}

/************************************************************************************
 Method: hNeighborhood
 Description: Consider the set of points in v that are integer steps (of size h) away 
 from x. Define the projection of the points in v onto the hyper-sphere centered at 
 x of radius h. The h-neighborhood of the point x is defined as the set of points 
 in Bh. The algorithm randomly selects points in Bh, one at a time. 
*************************************************************************************/
static TSol hNeighborhood(TSol x, float h) //, float theta
{
    // direction vector
    std::vector<double> v;      
    v.resize(n);

    // define a random order for the neighors
    std::vector<int> RKorder(n);
    for (int i = 0; i < n; i++){
        RKorder[i] = i;
    }
    std::random_shuffle (RKorder.begin(), RKorder.end());

    double norm = 0.0;
    for (int i=0; i<n; i++)
    {
        if (randomico(0.0,1.0) <= 0.5)
            v[i] = irandomico(1, ceil((1.0 - x.vec[i].rk) / h) ); 
        else
            v[i] = -1 * irandomico(1, ceil((x.vec[i].rk) / h) ); 
        
        // v[i] = irandomico(-1, 1); 
        norm += pow(v[i] * h, 2);
    }
    norm = sqrt(norm);
    if (norm == 0) norm = 0.0001;

    // neighbor solution
    TSol xBest;
    xBest = x;
    xBest.ofv = INFINITY;
    int numIter = n*exp(-0.7);
    // #pragma omp parallel for num_threads(MAX_THREADS)
    for (int i=0; i<numIter; i++)
    {
        x.vec[RKorder[i]].rk = x.vec[RKorder[i]].rk + (1.0/norm) * h * v[i] * h;      

        // if value without [0,1), create a new random key
        if (x.vec[RKorder[i]].rk < 0 || x.vec[RKorder[i]].rk >= 1.0){
            x.vec[RKorder[i]].rk = randomico(0,1);
        }

        Decoder(x);
        if (x.ofv < xBest.ofv)
            xBest = x;
        else
            x = xBest;
    }
    
    v.clear();
    return xBest;
}

/************************************************************************************
 Method: GridSearch
 Description: Generates a neighborhood and determines at which points in the neighbor-
 hood, if any, the objective function improves. If an improving point is found, it is 
 made the current point and the local search continues from the new solution.
*************************************************************************************/
static void GridSearch(TSol &x, float h) //, float theta
{
    int numGridPoints = 1; //floor(n *  (1.0/h));
    int numExaminedPoints = 0;

    // set the best solution found as current solution x
    TSol xBest = x;
    // theta = 1.0; 

    while (numExaminedPoints < numGridPoints) 
    {   
        numExaminedPoints++;

        // create a neighbor solution in the h-Neighborhood
        TSol y = hNeighborhood(xBest, h); 

        // printf("\n%lf \t %lf", xBest.ofv, y.ofv);

        if (y.ofv < xBest.ofv){
            xBest = y;
            // numExaminedPoints = 0;
        }
    }

    // retornar a melhor solucao
    x = xBest;
}

/************************************************************************************
 Method: UX
 Description: uniform crossover
*************************************************************************************/
static TSol UX(TSol &s1, TSol &s2, double factor)
{	
	TSol s;

	// create a new solution
	s.vec.resize(n+1);

    // Mate: including decoder gene in the n-th rk 
    for(int j = 0; j < n; j++)
    {
        // mutation
        if (randomico(0,1) < 0.02){
            s.vec[j].rk = randomico(0,1);
        }

        //copy alelos of top chromossom of the new generation
        else{
            if (randomico(0,1) < 0.5){
                s.vec[j].rk = s1.vec[j].rk;
            }
            else{
                if (factor == -1)
                    s.vec[j].rk = 1.0 - s2.vec[j].rk;
                else
                    s.vec[j].rk = s2.vec[j].rk;
            }
        }
    }

    // set the flag of local search as zero
    s.flag = 0;
    return s;
}

/************************************************************************************
 Method: NelderMeadSearch
 Description: The Nelder–Mead method is a numerical method used to find the minimum 
 of an objective function in a multidimensional space. It is a direct search method 
 based on function comparison.
*************************************************************************************/
static void NelderMeadSearch(TSol &x1)
{
    // elite points
    int k1, k2;
    do {
        k1 = irandomico(0,pool.size()-1);
        k2 = irandomico(0,pool.size()-1);
    }
    while (k1 == k2);
    TSol x2 = pool[k1];
    TSol x3 = pool[k2];

    // internal points
    TSol x_r;
    TSol x_e;
    TSol x_c;
    TSol x0;

    TSol xBest = x1;

    int iter_count = 0;
    int eval_count = 0;

    // sort points in the simplex so that x1 is the point having
    // minimum fx and x3 is the one having the maximum fx
    if (x1.ofv > x2.ofv) {
        TSol temp = x1;
        x1 = x2;
        x2 = temp;
    }

    if (x1.ofv > x3.ofv) {
        TSol temp = x1;
        x1 = x3;
        x3 = temp;
    }

    if (x2.ofv > x3.ofv) {
        TSol temp = x2;
        x2 = x3;
        x3 = temp;
    }
    
    // compute the simplex centroid
    x0 = UX(x1, x2, 1);
    Decoder(x0);
    if (x0.ofv < xBest.ofv) xBest = x0;

    iter_count++; 

    // continue minimization until stop conditions are met
    int maxIter = n*exp(-2);;
    while (iter_count <= maxIter) 
    {
        int shrink = 0;

        // if (debug) printf("\nIteration %04d     ", iter_count); //getchar();

        // reflection point (r)
        x_r = UX(x0, x3, -1);
        Decoder(x_r);
        if (x_r.ofv < xBest.ofv) xBest = x_r;
        eval_count++;

        // point_r is better than the bestSolution
        if (x_r.ofv < x1.ofv) 
        {
            // expansion point (e)
            x_e = UX(x_r, x0, -1);
            Decoder(x_e);
            if (x_e.ofv < xBest.ofv) xBest = x_e;
            eval_count++;

            if (x_e.ofv < x_r.ofv) 
            {
                // expand
                // if (debug) printf("\nexpand          ");
                x3 = x_e;
            } 
            else 
            {
                // reflect
                // if (debug) printf("\nreflect         ");
                x3 = x_r;
            }
        } 
        // x_r is NOT better than the bestSolution
        else 
        {    
            // point_r is better than the second best solution
            if (x_r.ofv < x2.ofv) 
            {
                // reflect
                // if (debug) printf("\nreflect         ");
                x3 = x_r;
            } 
            else 
            {
                // point_r is better than the worst solution
                if (x_r.ofv < x3.ofv) 
                {
                    // contraction point (c)
                    x_c = UX(x_r, x0, 1);
                    Decoder(x_c);
                    if (x_c.ofv < xBest.ofv) xBest = x_c;
                    eval_count++;

                    if (x_c.ofv < x_r.ofv) 
                    {
                        // contract outside
                        // if (debug) printf("\ncontract out    ");
                        x3 = x_c;
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
                    x_c = UX(x0, x3, 1);
                    Decoder(x_c);
                    if (x_c.ofv < xBest.ofv) xBest = x_c;
                    eval_count++;

                    if (x_c.ofv < x3.ofv) 
                    {
                        // contract inside
                        // if (debug) printf("\ncontract in     ");
                        x3 = x_c;
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
            x2 = UX(x1, x2, 1);
            Decoder(x2);
            // if (x2.ofv < xBest.ofv) xBest = x2;
            eval_count++;

            x3 = UX(x1, x3, 1);
            Decoder(x3);
            // if (x3.ofv < xBest.ofv) xBest = x3;
            eval_count++;
        }

        // sort
        if (x1.ofv > x2.ofv) {
            TSol temp = x1;
            x1 = x2;
            x2 = temp;
        }

        if (x1.ofv > x3.ofv) {
            TSol temp = x1;
            x1 = x3;
            x3 = temp;
        }

        if (x2.ofv > x3.ofv) {
            TSol temp = x2;
            x2 = x3;
            x3 = temp;
        }

        // compute the simplex centroid
        x0 = UX(x1, x2, 1);
        Decoder(x0);
        if (x0.ofv < xBest.ofv) xBest = x0;

        iter_count++;
    }

    // debug = 1;
    x1 = xBest;
}

/************************************************************************************
 Method: SwapLS
 Description: swap local search
*************************************************************************************/
static void SwapLS(TSol &s)
{
    // define a random order for the neighors
    std::vector<int> RKorder(n);
    for (int i = 0; i < n; i++){
        RKorder[i] = i;
    }
    std::random_shuffle(RKorder.begin(), RKorder.end());
    
    TSol sBest = s;
    // #pragma omp parallel for num_threads(MAX_THREADS)
    for(int i = 0; i < n-1; i++) {
        for(int j = i+1; j < n; j++) {        
            // Swap positions i and j
            double temp = s.vec[RKorder[i]].rk;
            s.vec[RKorder[i]].rk = s.vec[RKorder[j]].rk;
            s.vec[RKorder[j]].rk = temp;

            Decoder(s);

            if (s.ofv <= sBest.ofv){
                sBest = s;
            }
            else{
                s = sBest;
            }
        }
    }
}

/************************************************************************************
 Method: InvertLS
 Description: Invert local search
*************************************************************************************/
static void InvertLS(TSol &s)
{
    // define a random order for the neighors
    std::vector<int> RKorder(n);
    for (int i = 0; i < n; i++){
        RKorder[i] = i;
    }
    std::random_shuffle (RKorder.begin(), RKorder.end());

    TSol sBest = s;
    // #pragma omp parallel for num_threads(MAX_THREADS)
    for(int i = 0; i < n; i++) {
        // invert the random-key value
        if (s.vec[RKorder[i]].rk > 0.00001)
            s.vec[RKorder[i]].rk = 1.0 - s.vec[RKorder[i]].rk;
        else
            s.vec[RKorder[i]].rk = 0.99999;

        Decoder(s);

        if (s.ofv <= sBest.ofv){
            sBest = s;
        }
        else{
            s = sBest;
        }
    }
}

/************************************************************************************
 Method: Farey LS
 Description: Farey local search
*************************************************************************************/
static void FareyLS(TSol &s)
{
    // define a random order for the neighors
    std::vector<int> RKorder(n);
    for (int i = 0; i < n; i++){
        RKorder[i] = i;
    }
    std::random_shuffle (RKorder.begin(), RKorder.end());

    std::vector<double> F = {0.00, 0.142857, 0.166667, 0.20, 0.25, 0.285714, 0.333333, 0.40, 0.428571, 0.50, 
                             0.571429, 0.60, 0.666667, 0.714286, 0.75, 0.80, 0.833333, 0.857143, 0.99999};
    TSol sBest = s;
    // #pragma omp parallel for num_threads(MAX_THREADS)
    for(int i = 0; i < n; i++) {
        for (int j=0; j<(int)F.size()-1; j++){
            // gerar um valor aleatorio entre dois intervalos da sequencia de Farey
            s.vec[RKorder[i]].rk = randomico(F[j], F[j+1]);

            Decoder(s);

            if (s.ofv <= sBest.ofv){
                sBest = s;
            }
            else{
                s = sBest;
            }
        }
    }
}

/************************************************************************************
 Method: RVND
 Description: Random Variable Neighborhood Descent
*************************************************************************************/
void RVND(TSol &s, float h)
{
    // ***** we use a Random Variable Neighborhood Descent (RVND) as local search ****
    int numLS = 5;

    // current neighborhood
	int k = 1;

    // predefined number of neighborhood moves
    std::vector <int> NSL;
    std::vector <int> NSLAux;
    
    for (int i=1; i<=numLS; i++)
    {
        NSL.push_back(i);
        NSLAux.push_back(i);
    }

    int numIter = 0;
	while (!NSL.empty())
	{
        // current objective function
        double foCurrent = s.ofv;
        numIter++;

        // randomly choose a neighborhood
        int pos = rand() % NSL.size();
        k = NSL[pos];

        switch (k)
        {
            case 1: 
                SwapLS(s); // LS1(s); 
                break;

            case 2: 
                InvertLS(s); // LS2(s); 
                break;

            case 3: 
                FareyLS(s); // LS3(s); 
                break;

            case 4: 
                GridSearch(s, h);// LS4(s); 
                break;

            case 5: 
                NelderMeadSearch(s); // LS5(s); 
                break;

            default:
                break;
        }

        // we better the current solution
        if (s.ofv < foCurrent){
            // refresh NSL
            NSL.clear();
            NSL = NSLAux;
        }
        else{
            // Remove N(n) from NSL
            NSL.erase(NSL.begin()+pos);
        }
	} //end while
}

/************************************************************************************
 Method: CretePoolSolutions
 Description: create a pool of solutions with different solutions
*************************************************************************************/
static void CretePoolSolutions()
{
    pool.resize(sizePool);
    #pragma omp parallel for num_threads(MAX_THREADS) 
    for (int i = 0; i < sizePool; i++){
        CreateInitialSolutions(pool[i]);
        Decoder(pool[i]);
    }
    
    // apply local search in the pool solutions
    #pragma omp parallel for num_threads(MAX_THREADS)
    for (int i = 0; i < sizePool; i++){
        // RVND(pool[i],0.00781);
        FareyLS(pool[i]);
    }

    // sort pool in increase order of fitness
    sort(pool.begin(), pool.end(), sortByFitness);

    // verify if exists similar solutions in the pool
    for (int i = sizePool-1; i >0 ; i--){
        if (pool[i].ofv == pool[i-1].ofv){
            ShakeSolution(pool[i], 0.2, 0.5);
            Decoder(pool[i]);
        }
    }
    sort(pool.begin(), pool.end(), sortByFitness);
    
    updateBestSolution(pool[0],bestSolution,Tbest);
}

/************************************************************************************
 Method: UpdatePoolSolutions
 Description: update the pool with different solutions
*************************************************************************************/
static void UpdatePoolSolutions(TSol s)
{
    int difSol = 1;
    for (int i = pool.size()-1; i >= 0; i--){
        if (s.ofv == pool[i].ofv){
            difSol = 0;
            break;
        }
    }
    if (s.ofv < pool[sizePool-1].ofv && difSol == 1){
        pool[sizePool-1] = s;

        // sort pool in increase order of fitness
        sort(pool.begin(), pool.end(), sortByFitness);
    }
}

#endif
