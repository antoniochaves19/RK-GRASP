#ifndef _GRASP_H
#define _GRASP_H

// Variables declared in the main
extern int debug;                           // 0 - run mode      		    1 - debug mode
extern int numDecoders;                     // number of decoders
extern int MAXTIME;                         // maximum runtime
extern float OPTIMAL;                       // optimal solution (if it is known)
extern struct timeval Tstart, Tend, Tbest;  // computational time (unix systems)  

extern int n;                               // size of random-key vector
extern TSol bestSolution;                   // best solution found for any method
extern std::vector <TSol> pool;             // pool of best solutions with diversity

/************************************************************************************
			                  Implementation
*************************************************************************************/

/************************************************************************************
 Method: LineSearch
 Description: Consider the line search in direction ei , where vector ei has zeros 
 in all components except the i-th, where it has value one. The objective function 
 is evaluated at points x + k·h·ei for k = 0,1,−1,2,−2,... such that li≤xi+k·h≤ui.
 Let k∗ the value of k that minimizes f(x+k·h·ei) subject to li ≤ xi +k·h ≤ ui
*************************************************************************************/
static void LineSearch(TSol s, float h, int i, double &bestZ, double &bestF)
{
    std::vector<double> rk;

    // encontrar a melhor solucao na linha
    bestZ = 0;
    bestF = INFINITY;

    // encontrar a melhor solucao na linha (partindo do ponto corrente)
    // bestZ = s.vec[i].rk;
    // bestF = s.ofv;

    // gerar k como sendo as possiveis random keys do gene i, k = rk_i * tau | tau = {0, 1, -1, 2, -2, ...}
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
    int q = ceil(log2((int)(1.0/h))) + 1; 

    if (q > (int)rk.size())
        q = rk.size();

    // escolher um subconjunto com q rks para calcular o decoder
    std::random_shuffle (rk.begin(), rk.end());

    // calcular a qualidade da solucao s com a rk j
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

/************************************************************************************
 Method: ConstrutiveGreedyRandomized
 Description: It takes as input a solution vector s. Initially, the procedure allows 
 all coordinates of s to change (i.e. they are called unfixed). In turn, a discrete 
 line search is performed in each unfixed coordinate direction i of s with the other 
 n − 1 coordinates of s held at their current values.
*************************************************************************************/
static void ConstrutiveGreedyRandomized(TSol &s, float h, float alpha)
{
    // ** nao estou alterando o decoder
    // A fase construtiva inicia da solução corrente s e fará uma “perturbação” desta solução com uma intensidade entre Beta_min e Beta_max. 
    // Em cada iteração da fase construtiva escolhe-se aleatoriamente um conjunto de random keys que ainda não foram fixadas e para cada 
    // random key deste conjunto gera-se todos os possíveis valores de novas random keys baseado no valor de h. Escolhe-se aleatoriamente 
    // um destes valores e calcula o decoder (esta será a saída do line search). Com a lista de novas random keys geradas, construimos a RCL 
    // e selecionamos um valor aleatoriamente entre os melhores, fixando esta random key.

    std::vector<int> UnFixed(n);                // armazena as random-keys ainda nao fixadas
    std::vector<int> chosenRK;                  // armazena as random-keys que serao pesquisadas no line search
    std::vector<int> RCL;                       // armazena as melhores solucoes candidatas
    std::vector<double> z(n);                   // armazena o melhor valor da random-key i
    std::vector<double> g(n,INFINITY);          // armazena o valor da fo com a random-key z_i

    double min, max;
    double betaMin = 0.5, 
           betaMax = 0.8;

    // inicializar os pontos do cromossomo que podem ser alterados
    for (int i = 0; i < n; i++) {
        UnFixed[i] = i;
    }

    // construir uma solucao por meio de perturbacoes na solucao corrente e escolha de uma das melhores
    double intensity = randomico(betaMin, betaMax);
    
    // while (!UnFixed.empty())
    for (int j=0; j<n*intensity; j++)
    {
        // criar uma lista de solucoes candidatas perturbando uma rk (ainda não 'fixada') da solucao corrente
        min = INFINITY;
        max = -INFINITY;

        // escolher o subconjunto de random keys que serao pesquisadas
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

        #pragma omp parallel for num_threads(MAX_THREADS) //private(s)
        for (int k=0; k<kMax; k++) 
        {
            int i = chosenRK[k];

            TSol sAux = s;
            
            // linear search
            LineSearch(sAux, h, i, z[i], g[i]);

            // #pragma omp critical
            // {
            //     if (min > g[i])
            //         min = g[i];

            //     if (max < g[i])
            //         max = g[i];
            // }
        }

        for(int i=0; i<n; i++)
        {
            if (min > g[i])
                min = g[i];

            if (max < g[i] && g[i] != INFINITY)
                max = g[i];
        }

        // criar a RCL
        RCL.clear();
        double threshold = min + alpha * (max - min);

        for (int i=0; i<n; i++)
        { 
            if (g[i] <= threshold)
            {
                RCL.push_back(i);
            }
        }

        // selecionar aleatoriamente um dos melhores candidatos para continuar a construcao
        if (!RCL.empty()){
            int x = irandomico(0, (int)(RCL.size())-1);
            int kCurrent = RCL[x];  // indice da random key que sera fixado
            
            // atualizar a solucao corrente
            s.vec[kCurrent].rk = z[kCurrent];
            s.ofv = g[kCurrent];
            // printf("\tFoViz %d = %lf (%d, %lf)\n", j, s.ofv, k, z[k]);

            // retirar a rk k do conjunto UnFixed
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

/************************************************************************************
 Method: GRASP
 Description: Metaheurist Greedy Randomized Adpative Search Procedura.
*************************************************************************************/
extern void GRASP(int method, int control)
{
    // GRASP parameters
    static TSol s;                              // current solution
    static TSol sLine;                          // constructive solution
    static TSol sLineBest;                      // local search solution
    static TSol sBest;                          // best solution of RK-GRASP

    static double sigma;                        // greedy rate
    static float theta;                         // intensity perturbation in GridSearch

    static float h;                             // grid dense
    static float hs = 0.12500;                  // start grid dense
    static float he = 0.00098;                  // end grid dense

    int ls = 0;                                 // local search method 

    // computational time of the search process
    float currentTime = 0;       

    // offline control
    if (control == 0){
        // define parameters of GRASP
        sigma = randomico(0.1, 0.9);
        // ls = busca;
    }

    int iter = 0;

    // create an initial solution
    CreateInitialSolutions(s);
    Decoder(s);

    sBest = s;
    updateBestSolution(sBest, bestSolution, Tbest);

    // run the search process until stop criterion (maxTime)
    while(1)
    {
        h = hs;
        while (h >= he)
        {
            iter++;

            // offline control
            if (control == 0){
                sigma = randomico(0.1, 0.9);
            }

            // construct a greedy randomized solution
            sLine = s;
            ConstrutiveGreedyRandomized(sLine, h, sigma);

            if (debug) printf("\nIter = %d \t h = %.4lf \t Constructive (%.1lf) = %.2lf", iter, h, sigma, floorf(sLine.ofv*10)/10);
            
            // apply local search in current solution
            sLineBest = sLine;
            RVND(sLineBest, h);

            if (debug) printf("\t Local search (%d, %.3lf) = %.5lf (%8.2lf)", ls, theta, floorf(sLineBest.ofv*10)/10, sLine.ofv - sLineBest.ofv);

            // update the best solution found by GRASP
            if (sLineBest.ofv < sBest.ofv){
                sBest = sLineBest;
                
                // update the best global solution found
                updateBestSolution(sBest, bestSolution, Tbest);
                if (debug) printf("**");
            }
            // make grid more dense
            else{
                h = h/2;
            }

            // update the pool of solutions
            UpdatePoolSolutions(sLineBest);

            // accept criterion
            if (sLineBest.ofv < s.ofv){
                s = sLineBest;
            }
            else{
                // Metropolis criterion
                if ( randomico(0.0,1.0) < (exp(-(sLineBest.ofv - s.ofv)/(1000 - 1000*(currentTime / MAXTIME)))) ){ 
                    s = sLineBest;
                }
            }

            // if (debug) printf("\t Best solution = %lf", floorf(bestSolution.ofv*10)/10);
            if (debug) printf("\t Best solution = %lf", bestSolution.ofv);

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
        // CreateInitialSolutions(s); 

        // shake the current solution
        ShakeSolution(s,0.2,0.5);

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
}

#endif