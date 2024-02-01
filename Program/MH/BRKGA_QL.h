#ifndef _BRKGA_QL_H
#define _BRKGA_QL_H

// Variables declared in the main
extern int debug;                           // 0 - run mode      		    1 - debug mode
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int MAXTIME;                         // maximum runtime
extern float OPTIMAL;                       // optimal solution (if it is known)
extern struct timeval Tstart, Tend, Tbest;  // computational time (unix systems)  
extern unsigned MAX_THREADS;                // number of threads

extern int n;                                      // size of cromossoms
extern TSol bestSolution;                          // best solution found in the BRKGA-QL


// Q-Learning parameters
static double epsilon;                             // greed choice possibility
static double lf;                                  // learning factor
static double df;                                  // discount factor
static double R;                                   // reward
static double qTotal;                              // q*

// list of actions RL
static int sizeP []    = {233, 377, 610, 987, 1597, 2584};
static double Pe[]     = {0.10, 0.15, 0.20, 0.25, 0.30}; 
static double Pm[]     = {0.01, 0.02, 0.03, 0.04, 0.05}; 
static double Rhoe[]   = {0.55, 0.60, 0.65, 0.70, 0.75, 0.80}; 

// number of parameters in Q-table
static const int par = 4;

// actions
static int a0 = 0;                                 // p (current action)
static int a1 = 0;                                 // pe (current action)
static int a2 = 0;                                 // pm (current action)
static int a3 = 0;                                 // rhoe (current action)

static float Qmax = 0;

static std::vector <std::vector <TQ> > Q;          // Q-Table

// QL parameters (auxiliar)
static float epsilon_max = 1.0;                    // maximum epsilon 
static float epsilon_min = 0.1;                    // minimum epsilon
static int Ti = 1;                                 // number of epochs performed
static int restartEpsilon = 1;                     // number of restart epsilon



/************************************************************************************
			                  DECLARATION
*************************************************************************************/

/************************************************************************************
 Method: BRKGA_QL
 Description: the evolutionary process of the BRKGA-QL
*************************************************************************************/
static void BRKGA_QL(char nameTable[256]);

/************************************************************************************
 Method: InitiateQTable()
 Description: Initiate the Q-Table with random values
*************************************************************************************/
static void InitiateQTable();

/************************************************************************************
 Method: ChooseAction()
 Description: Choose actions and update the parameters
*************************************************************************************/
static void ChooseAction(int &p, double &pe, double &pm, double &rhoe);

/************************************************************************************
 Method: UpdatePopulationSize()
 Description: Update the population size with new value of p
*************************************************************************************/
static void UpdatePopulationSize(int p, double pe, double pm, double rhoe, std::vector <TSol> &Pop, std::vector <TSol> &PopInter);

/************************************************************************************
 Method: UpdateQLParameters(currentTime)
 Description: Update the parameters of the Q-Learning method
*************************************************************************************/
static void SetQLParameters(float currentTime);

/************************************************************************************
 Method: UpdateQTable()
 Description: Update the values of Q-Table
*************************************************************************************/
static void UpdateQTable();

/************************************************************************************
 Method: ChaoticInd
 Description: create a solution between a mutant and a elite individual
*************************************************************************************/
static void ChaoticInd(TSol &s, int rhoe);

/************************************************************************************
 Method: PARAMETRICUNIFORMCROSSOVER
 Description: create a new offspring with parametric uniform crossover
*************************************************************************************/
static TSol ParametricUniformCrossover(int elitesize, int popSize, double pm, double rhoe, std::vector <TSol> &Pop);

/************************************************************************************
 Method: PEARSON CORRELATION
 Description: calculate the Pearson correlation coefficient between two chromossoms
*************************************************************************************/
static double PearsonCorrelation(std::vector <TVecSol> s1, std::vector <TVecSol> s2);

/************************************************************************************
 Metodo: IC(TSol Pop)
 Description: apply clustering method to find promising solutions in the population
*************************************************************************************/
static void IC(int p, double pe, std::vector <TSol> &Pop);

/************************************************************************************
 Method: LP
 Description: Apply Label Propagation to find communities in the population
*************************************************************************************/
static void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas, std::vector <TSol> &Pop);

/************************************************************************************
 Method: PROMISINGLP
 Description: Find the promising solutions to represent the communities
*************************************************************************************/
static void PromisingLP(int p, double pe, std::vector <TSol> &Pop);


/************************************************************************************
			                  IMPLEMENTATION
*************************************************************************************/

static void BRKGA_QL(char nameTable[256]) 
{
    // BRKGA parameters
    static int p;          	                        // size of population
    static double pe;              	                // fraction of population to be the elite-set
    static double pm;          	                    // fraction of population to be replaced by mutants
    static double rhoe;             	            // probability that offspring inherit an allele from elite parent

    // BRKGA variables
    static std::vector <TSol> Pop;                  // current population
    static std::vector <TSol> PopInter;             // intermediary population

    // free memory with problem data
    FreeMemoryProblem();

    //read data of the instance
    ReadData(nameTable);
 
    // initialize Q-Table
    InitiateQTable();

    // number of restart epsilon
    restartEpsilon = 1;  

    // maximum epsilon  
    epsilon_max = 1.0;  
    
    // initialize population
    Pop.clear();  
    PopInter.clear(); 

    // define the population size with the higher value of P
    p = sizeP[sizeof(sizeP)/sizeof(sizeP[0]) - 1]; 

    Pop.resize(p);
    PopInter.resize(p);

    // Create the initial chromosomes with random keys
    #pragma omp parallel for num_threads(MAX_THREADS)
    for (int i=0; i<p; i++)
    {
        TSol ind;
        CreateInitialSolutions(ind); 
        Decoder(ind);
        Pop[i] = PopInter[i] = ind;
    }
    
    // sort population in increase order of fitness
    sort(Pop.begin(), Pop.end(), sortByFitness);

    // save the best solution found
    updateBestSolution(Pop[0], bestSolution, Tbest);
    
    // useful local variables
    int numGenerations = 0;             // number of generations
    int bestGeneration = 0;             // generation in which found the best solution
    double bestFitness = Pop[0].ofv;    // best fitness found in past generation
    float currentTime = 0;              // computational time of the search process
    int sumLS = 0;                      // number of local search applied in each generation
    int noImprov = 0;                   // number of generations without improvement in the best solution
    
    // run the evolutionary process until stop criterion
    while(1)
    {
    	// number of generations
        numGenerations++;

        // number of generations without improvement in the best solution
        noImprov++;

        // *************************** BRKGA-QL **************************************
        // set Q-Learning parameters                                              //**
        SetQLParameters(currentTime);                                             //**
        //                                                                        //**
        // choose a action (value) for each state (parameter)                     //**
        ChooseAction(p, pe, pm, rhoe);                                            //**
        //                                                                        //**
        // update population size                                                 //**
        UpdatePopulationSize(p, pe, pm, rhoe, Pop, PopInter);                                    //**
        // ***************************************************************************


        // **************************** BRKGA ****************************************
        // define parameters for classic BRKGA                                    //**
        //p       = 987;                                                          //**
        //pe      = 0.20;                                                         //**
        //pm      = 0.05;                                                         //**
        //rhoe    = 0.70;                                                         //** 
        // ***************************************************************************

        // if (debug){
        //     FILE *arquivo;
        //     arquivo = fopen("../Results/Parametros.csv","a");
        //     fprintf(arquivo, "\n%d \t %d \t %.3lf \t %.3lf \t %.3lf",numGenerations, p, pe, pm, rhoe);
        //     fclose(arquivo);
        // }


        // The 'Pe' best chromosomes are maintained, so we just copy these into PopInter:
        for (int i=0; i<(int)(p*pe); i++){
            // copy the chromosome for next generation
            PopInter[i] = Pop[i]; 
        }  

        // We'll mate 'P - Pe' pairs; initially, i = p*pe, so we need to iterate until i < p:
        #pragma omp parallel for num_threads(MAX_THREADS)
        for (int i = (int)(p*pe); i < p; i++){            
            // Parametric uniform crossover with mutation
            PopInter[i] = ParametricUniformCrossover((int)(p*pe), p, pm, rhoe, Pop);
 
            // Calculate the fitness of new chromosomes
            Decoder(PopInter[i]); 
        }
                
        // Update the current population
        Pop = PopInter;   

        // Sort population in increase order of fitness
        sort(Pop.begin(), Pop.end(), sortByFitness);
        updateBestSolution(Pop[0], bestSolution, Tbest);

        // Reward
        R = 0;

        // We improve the best fitness in the current population 
        if (Pop[0].ofv < bestFitness){

            // The reward function is based on improvement of the current best fitness and binary reward
             R = 1 + 100000*(((bestFitness/Pop[0].ofv)-1)/(p));

            // The immediate reward values are limited by some constant.
            if (R > 2) R = 2;

            bestFitness = Pop[0].ofv;
            bestGeneration = numGenerations;
            noImprov = 0;
        }

        // Print reward of each generation
        // if (debug) {
        //     FILE *arquivo;
        //     arquivo = fopen("../Results/Reward.txt","a");
        //     fprintf(arquivo, "\n%d \t %.3lf",numGenerations, R);
        //     fclose(arquivo);
        // }
        
        // Update the Q-Table values
        if (R > 0){
            UpdateQTable();
        }


        // ********************* LOCAL SEARCH IN COMMUNITIES *******************
        sumLS = 0;

        // Verify if there are local search heuristics 
        if (numLS > 0){   

            //apply local search when BRKGA found a new better solution or n*pe generations without improvements
            if (R >= 1 || noImprov > (int)n*pe){

                // restart the count of generations without improvements (local search)
                noImprov = 0;

                // Identify commuties in the Elite with Label Propagation method
	            IC(p, pe, Pop);

	            std::vector <int> promisingSol; 
                promisingSol.clear();

	            for (int i=0; i < (int)(p*pe); i++) {
                    
                    // insert the individual index in the promising list
	                if (Pop[i].promising == 1){
	                	promisingSol.push_back(i);
	                }
                    
                    // generate caotic individual (crossover between one elite and one mutant)
                    else if (i > 0){
                        ChaoticInd(Pop[i], rhoe);
                        Decoder(Pop[i]);

                        // set flag as 0 to permit new local search
                        Pop[i].flag = 0;
                    }
	            }

	            #pragma omp parallel for num_threads(MAX_THREADS)
                for (unsigned int i=0; i < promisingSol.size(); i++){

                    // local search not influence the evolutionary process
                    LocalSearch(Pop[promisingSol[i]]);

                    // set flag as 1 to prevent new local search in the same solution
                    Pop[promisingSol[i]].flag = 1;

                    #pragma omp critical
                    {   
                        updateBestSolution(Pop[promisingSol[i]], bestSolution, Tbest);                    
                    }
                }

                // if (promisingSol.size() > 1)
                //     PR(Pop[promisingSol[0]],Pop[promisingSol[1]]);

                sumLS = promisingSol.size();
                promisingSol.clear();

                sort(Pop.begin(), Pop.end(), sortByFitness);
                updateBestSolution(Pop[0], bestSolution, Tbest);
	        }
	    }
        // *********************************************************************


        // ******************************* SHAKING *****************************
        if ((numGenerations - bestGeneration) > 5*n) {
            
            if (debug) 
                printf("\n\nShaking elite and reset non-elite...\n\n");
            else
                srand(time(NULL)); 

            // reset the number of generations without improvement
            bestGeneration = numGenerations;

            // Shake the elite set
            float shaking_type = 0.0;
            int intensity = n*rhoe;
            #pragma omp parallel for num_threads(MAX_THREADS)
            for(int e = 0; e < (int)(p*pe); e++) {
                for(int k = 0; k < intensity; k++) {
                    shaking_type = irandomico(1,4);
                    int i = irandomico(0, n - 1);
                    if(shaking_type == 1){
                        // Invert value
                        if (Pop[e].vec[i].rk > 0.0001)
                            Pop[e].vec[i].rk = 1.0 - Pop[e].vec[i].rk;
                        else
                            Pop[e].vec[i].rk = 0.9999;
                    }
                    else 
                    if (shaking_type == 2){
                        // Swap two random positions
                        int j = irandomico(0, n - 1);
                        double temp = Pop[e].vec[i].rk;
                        Pop[e].vec[i].rk = Pop[e].vec[j].rk;
                        Pop[e].vec[j].rk = temp;
                    }
                    else
                    if(shaking_type == 3){
                        // Change to random value
                        Pop[e].vec[i].rk = randomico(0,1);
                    }
                    i = irandomico(0, n - 2);
                    if(shaking_type == 4){
                        // Swap with neighbor
                        double temp = Pop[e].vec[i].rk;
                        Pop[e].vec[i].rk = Pop[e].vec[i+1].rk;
                        Pop[e].vec[i+1].rk = temp;
                    }
                }
                Decoder(Pop[e]);
            }

            // reset the non-elite chromosomes
            #pragma omp parallel for num_threads(MAX_THREADS)
            for (int i=(int)(p*pe); i<p; i++){
                CreateInitialSolutions(Pop[i]);
                Decoder(Pop[i]);
            }

            sort(Pop.begin(), Pop.end(), sortByFitness);
            updateBestSolution(Pop[0], bestSolution, Tbest);
            bestFitness = Pop[0].ofv;
        }
        // *********************************************************************

        // print screen 
        if (debug){
            printf("\nGeneration: %3d [%4d - %3d(%.2lf) (%3d) - (%.2lf) - (%.2lf)] \t %.1lf (%.2lf)  \t %.1lf [%.4lf] \t %.4lf \t %.4lf",
                        numGenerations, p, (int)(p*pe), pe, sumLS, pm, rhoe, bestSolution.ofv, bestSolution.vec[n].rk, bestFitness, R, epsilon, lf);
        }

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

    // free memory of BRKGA-QL components
    Pop.clear();
    PopInter.clear();
    Q.clear();
}

static void InitiateQTable()
{
    // initialize the Q-Table values at 0
    Q.clear();
    Q.resize(par);

    qTotal = 0.0;

    // rho_e
    for (unsigned int j=0; j<sizeof(Rhoe)/sizeof(Rhoe[0]); j++)
    {
        TQ aux;
        aux.S = 0;
        aux.pVar = Rhoe[j];
        aux.q = 0;
        aux.k = 0;
        aux.kImp = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    // p
    for (unsigned int j=0; j<sizeof(sizeP)/sizeof(sizeP[0]); j++)
    {
        TQ aux;
        aux.S = 1;
        aux.pVar = sizeP[j];
        aux.q = 0;
        aux.k = 0;
        aux.kImp = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    // pm
    for (unsigned int j=0; j<sizeof(Pm)/sizeof(Pm[0]); j++)
    {
        TQ aux;
        aux.S = 2;
        aux.pVar = Pm[j];
        aux.q = 0;
        aux.k = 0;
        aux.kImp = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    // pe
    for (unsigned int j=0; j<sizeof(Pe)/sizeof(Pe[0]); j++)
    {
        TQ aux;
        aux.S = 3;
        aux.pVar = Pe[j];
        aux.q = 0;
        aux.k = 0;
        aux.kImp = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }                                    
}

static void ChooseAction(int &p, double &pe, double &pm, double &rhoe)
{
    // choose actions for each state from Q-Table using epsilon-Greedy policy
    for (int i=0; i<par; i++)
    {
        int a = 0, 
            aAux = 0;

        // set variable a with the current action
        switch (i)
        {
            case 0:
                a = a0;
                break;
        
            case 1:
                a = a1;
                break;
            
            case 2:
                a = a2;
                break; 

            case 3:
                a = a3;
                break;
        }       
                
        // found actions with the highest value of Q(i,-).q 
        double bQ = -INFINITY;  
        for (unsigned int j=0; j<Q[i].size(); j++)
        {
            if (Q[i][j].q > bQ)
            {
                bQ = Q[i][j].q;
                aAux = j;
            }
            else
            if (Q[i][j].q == bQ && randomico(0,1) >= 0.5) // trie randomly
            {
                aAux = j;
            }

            // update the best future reward
            if (Q[i][j].q > Qmax)
                Qmax = Q[i][j].q;
        }

        // epsilon-greedy policy
        if (randomico(0,1) <= 1-epsilon) 
        {
            // choose the action with highest Q value
            a = aAux;
        }
        else
        {
            // choose a randonly selected action (value)
            a = irandomico(0,Q[i].size()-1);
        }

        // set new action
        switch (i)
        {
            case 0:
                a0 = a;
                break;
        
            case 1:
                a1 = a;
                break;
            
            case 2:
                a2 = a;
                break; 

            case 3:
                a3 = a;
                break;
        }

        // update number of choices state i and action a
        Q[i][a].k++;
    }
        
    // update parameters with actions 
    rhoe    = Q[0][a0].pVar;   
    p       = Q[1][a1].pVar;
    pm      = Q[2][a2].pVar;
    pe      = Q[3][a3].pVar;  
}

static void SetQLParameters(float currentTime)
{
    // **** define epsilon ****
    static const double PI = 3.14159265;            // pi

    // restart epsilon once Ti epochs are performed (Ti is 10% of the runtime)
    Ti = MAXTIME * 0.1;
    if (currentTime >= restartEpsilon * Ti){
        restartEpsilon++;

        // cosine decay with warm restart
        epsilon_max = epsilon_max - 0.1;
        if (epsilon_max < epsilon_min)
            epsilon_max = epsilon_min;
        epsilon = epsilon_max;
    }
    else {
        epsilon = epsilon_min + 0.5 * (epsilon_max - epsilon_min) * (1 + cos((((int)currentTime%Ti)/(float)(Ti))*PI));
    }
    
    // *** define learning rate ***

    // initialy, a higher priority is given to the newly gained information (exploration mode)
    // then, we decrement lf and have a higher priority for the existing information in Q-Table (exploitation mode)
    lf = 1 - (0.9 * currentTime / MAXTIME); 

    // *** define discount rate ***

    // we look for a higher, long-term reward
    df = 0.8;
}

static void UpdateQTable()
{ 
    // set the current action for each state
    for (int s=0; s<par; s++)
    {
        int a; 

        switch (s)
        {
            case 0:
                a = a0;
                break;
            
            case 1:
                a = a1;
                break;

            case 2:
                a = a2;
                break;

            case 3:
                a = a3;
                break;
        }
        
        qTotal -= Q[s][a].q;

        // Q-Learning
        // Q(s,a) is incremented when the action leads to a state, in which there exists an action such that the best possible Q-value and
        // the reward R is greater than current value of Q(s,a).
        // i.e., the old value of Q(s,a) was too pessimistic 
        // df*Qmax is the target Q-value

        Q[s][a].q = Q[s][a].q + lf*(R + df*Qmax - Q[s][a].q); 
       
        Q[s][a].kImp++;
        qTotal += Q[s][a].q;
    }
}

static void UpdatePopulationSize(int p, double pe, double pm, double rhoe, std::vector <TSol> &Pop, std::vector <TSol> &PopInter)
{
    // *** define the new population size

    // size of the current population
    int oldPsize = Pop.size();

    // proportional pruning 
    if (oldPsize > p){

        // copy the current population
        PopInter = Pop;

        // define new size of Pop
        Pop.resize(p);

        // select the elite chromosomes
        for (int i=0; i<(int)(p*pe); i++){
            // copy p*pe best chromosomes
            Pop[i] = PopInter[i];
        }

        // select the non-elite chromosomes
        int pos = (int)(pe*oldPsize);
        for (int i=(int)(p*pe); i<p; i++){
            // copy the chromosome
            Pop[i] = PopInter[pos];
            pos++;
        }

        // clean intermediate population
        PopInter.clear();
        PopInter.resize(p);
    }
    
    // generate new chromosomes 
    else if (oldPsize < p){

        // define new size of Pop
        Pop.resize(p);

        #pragma omp parallel for num_threads(MAX_THREADS)
        for (int k = oldPsize; k < p; k++)
        {
        	Pop[k] = ParametricUniformCrossover((int)(oldPsize*pe), oldPsize-1, pm, rhoe, Pop);
            Decoder(Pop[k]);
        }

        // sort new population
        sort(Pop.begin(), Pop.end(), sortByFitness);
        updateBestSolution(Pop[0], bestSolution, Tbest);
        
        // clean intermediate population
        PopInter.clear();
        PopInter.resize(p);
    }
}

static void ChaoticInd(TSol &s, int rhoe)
{
    // generate a caotic individual
    for (int k=0; k<n+1; k++)
    {      
        if (randomico(0,1) > rhoe)
           s.vec[k].rk = randomico(0,1);
    }

    // set the flag of local search as zero
    s.flag = 0;
}

static TSol ParametricUniformCrossover(int eliteSize, int popSize, double pm, double rhoe, std::vector <TSol> &Pop)
{	
	TSol s;

    int eliteParent = irandomico(0, eliteSize - 1);                 // one chromosome from elite set
    int nonEliteParent = irandomico(0, popSize-1);                  // one chromosome from entire population

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
        // mutation
        if (randomico(0,1) < pm)
        {
            s.vec[j].rk = randomico(0,1);
        }
        else
        {
            //copy alelos of top chromossom of the new generation
            if (randomico(0,1) < rhoe){
                s.vec[j].rk = Pop[eliteParent].vec[j].rk;
            }
            else{
                s.vec[j].rk = Pop[nonEliteParent].vec[j].rk;
            }
        }
    }

    // set the flag of local search as zero
    s.flag = 0;

    return s;
}

static double PearsonCorrelation(std::vector <TVecSol> X, std::vector <TVecSol> Y)
{
    double correlation = 0;
    double sumXY = 0;
    double sumX2 = 0;
    double sumY2 = 0;
    double sumX = 0;
    double sumY = 0;

    for(int j=0; j<n; j++)
    {
        sumX += X[j].rk;
        sumX2 += X[j].rk * X[j].rk;
        sumXY += X[j].rk * Y[j].rk;
        sumY += Y[j].rk;
        sumY2 += Y[j].rk * Y[j].rk;
    }

    //Pearson
    correlation= ((n*sumXY) - (sumX*sumY) ) / (sqrt( (n*sumX2 - sumX*sumX) * (n*sumY2 - sumY*sumY) ));
    return correlation;
}

static void IC(int p, double pe, std::vector <TSol> &Pop) 
{
    int Tpe = (int)p*pe;
    std::vector<std::vector<std::pair<int, double> > > listaArestas(Tpe, std::vector<std::pair<int, double> >());

    static double sigma = 0.6;                      // pearson correlation factor

	// create weighted (pearson correlation) graph
	int entrouAresta = 0;
	double pearson = 0.0;
	for (int i = 0; i < Tpe - 1; i++) {
		for (int j = i + 1; j < Tpe; j++)
		{
			pearson = PearsonCorrelation(Pop[i].vec, Pop[j].vec);
			if (pearson > sigma) {
				entrouAresta++;
				listaArestas[i].push_back(std::make_pair(j, pearson));
				listaArestas[j].push_back(std::make_pair(i, pearson));
			}
            else{
                entrouAresta += 5;
            }
		}
	}

	// apply clustering method
	LP(listaArestas, Pop);

	PromisingLP(p, pe, Pop);
    listaArestas.clear();
}

static void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas, std::vector <TSol> &Pop)
{
    int nk = listaArestas.size();

	// Create vector with visit order
	std::vector<int> ordemVisita(nk);
	iota(ordemVisita.begin(), ordemVisita.end(), 0);

	// initialize each node with its own label
	for (int i = 0; i < nk; i++)
		Pop[i].label = i;

	int iteracao = 1;
	int labelVizinho, melhorLabel;
	double melhorPeso;
	std::map<int, double> totalLabels;
	std::map<int, double>::iterator it;

	int movimentos = 1;
	while (movimentos) 
    {
		movimentos = 0;
		random_shuffle(ordemVisita.begin(), ordemVisita.end());
        for (std::vector<int>::size_type idVertice=0; idVertice <ordemVisita.size(); idVertice++)
        {
			// Calculate the weigth of the labels
			totalLabels.clear();
            for (std::vector<std::pair<int, double> >::iterator itVizinho = listaArestas[idVertice].begin();
                itVizinho != listaArestas[idVertice].end(); ++itVizinho) {
                int idVizinho = itVizinho->first;
                labelVizinho = Pop[idVizinho].label;
                it = totalLabels.find(labelVizinho);
                if (it != totalLabels.end()) {
                    it->second += itVizinho->second;
                } else {
                    totalLabels[labelVizinho] = itVizinho->second;
                }
            }

			// Best label is itself initially
			melhorLabel = Pop[idVertice].label;
			melhorPeso = std::numeric_limits<double>::min();
            for (std::map<int, double>::iterator itTotais = totalLabels.begin(); itTotais != totalLabels.end(); ++itTotais) {
                if (itTotais->second > melhorPeso) {
                    melhorLabel = itTotais->first;
                    melhorPeso = itTotais->second;
                }
            }

			if (melhorLabel != Pop[idVertice].label) {
				Pop[idVertice].label = melhorLabel;
				movimentos = 1;
			}
		}
		iteracao++;
	}

    ordemVisita.clear();
}

static void PromisingLP(int p, double pe, std::vector <TSol> &Pop)
{
    int Tpe = (int)p*pe;
    std::vector <int> grupos;
	int tamanhoGrupos = 0;

	// initialize promisings solutions
	for (int i = 0; i < Tpe; i++)
		Pop[i].promising = 0;

	// save labels defined by LP in groups
	int achei;

    for (int i = 0; i < Tpe; i++)
	{
		achei = 0;
		for (unsigned int j = 0; j < grupos.size(); j++)
		{
			if (Pop[i].label == grupos[j])
				achei = 1;
		}
		if (achei == 0)
		{
			tamanhoGrupos++;
			grupos.push_back(Pop[i].label);
		}
	}

	// find the best solution in the group (with flag = 0)
	for (unsigned int j = 0; j < grupos.size(); j++)
	{
		double menorFO = INFINITY;
		int localMenor = -1;
		int local = -1;
		for (int i = 0; i < Tpe; i++)
		{
			if (Pop[i].label == grupos[j])
			{
				// find the best solution of the group
				if (local == -1)
					local = i;

				// we not apply local search in this solution yet
                if (Pop[i].ofv < menorFO && Pop[i].flag == 0) 
				{
					menorFO = Pop[i].ofv;
					localMenor = i;
				}
			}
		}

		if (localMenor == -1)
			localMenor = local;

		if (Pop[localMenor].label != -1)
			Pop[localMenor].promising = 1;
	}
}


#endif