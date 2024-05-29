// *******************************************************************
//      file with specific functions to solve the HMP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Global Variables
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int n;                               // size of cromossoms

//---------------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC --------------------


//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------

static std::vector <std::vector <int> > hij;	// matrix with edge weight
static std::vector <double> ti;				    // vector of nodes weigth
static int numNodes;							// number of nodes
static int numClusters;							// number of clusters
static double cap;								// cluster capacity
static int totalHij;                            // total edge weigth


//----------------------- IMPLEMENTATION OF FUNCTIONS  -------------------------------


/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
*************************************************************************************/
void ReadData(char nameTable[])
{ 
    char name[200] = "../Instances/";
    strcat(name,nameTable);

    FILE *arq;
    arq = fopen(name,"r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n",name);
        getchar();
        exit(1);
    }

    // => read data
    fscanf(arq, "%d %d %lf", &numNodes, &numClusters, &cap);

    // read graph informations
    ti.clear();
    ti.resize(numNodes);

    hij.clear();
    hij.resize(numNodes, std::vector<int>(numNodes));

    for (int i=0; i<numNodes; i++)
    {
    	fscanf(arq, "%lf", &ti[i]);
    }

    totalHij = 0;
    for (int i=0; i<numNodes; i++)
    {
        for (int j=0; j<numNodes; j++)
        {
    	    fscanf(arq, "%d", &hij[i][j]);
            totalHij += hij[i][j];
        }
    }

    fclose(arq);
    
    n = numNodes+1; // a ultima posicao (n-1) representa o numero de base que adicionamos inicialmente em rncs diferentes


    // imprimir
    // printf("\n\n%d %d %lf \n", numNodes, numClusters, cap);

    // for (int i=0; i<numNodes; i++)
    // {
    // 	printf("%lf\n", ti[i]);
    // }

    // for (int i=0; i<numNodes; i++)
    // {   
    //     printf("\n");
    //     for (int j=0; j<numNodes; j++)
    //     {
    // 	    printf("%d\t", hij[i][j]);
    //     }
    // }
    // printf("\n\n");
    // getchar();
}

/************************************************************************************
 Method: CalculateFitness
 Description: calculate the fitness of a chromossom s
*************************************************************************************/
int CalculateFitness(std::vector< std::vector<int> > RNC)
{
    // calculate objective function
    int soma = 0;

    // somar os pesos de cada RNC
    // printf("\nRNC = %d; \t", RNC.size());
    for (int k=0; k<(int)RNC.size(); k++)
    {
        // printf("%d \t", RNC[k].size());
        // total handover between base stations assigned to the same RNC
        for (int i=0; i<(int)RNC[k].size(); i++)
        {
            for (int j=0; j<(int)RNC[k].size(); j++)
            {
                soma += hij[RNC[k][i]][RNC[k][j]];
            }
        }
    }
    return totalHij - soma;

    // for (int k=0; k<(int)RNC.size(); k++)
    // {
    //     for (int i=0; i<(int)RNC[k].size(); i++)
    //     {
    //         for (int k1=0; k1<(int)RNC.size(); k1++)
    //         {
    //             if (k1 != k)
    //             {
    //                 for (int j=0; j<(int)RNC[k1].size(); j++)
    //                 {
    //                     soma += hij[RNC[k][i]][RNC[k1][j]];
    //                 }
    //             }
    //         }
    //     }
    //     // soma += ((double)(rand()%10000)/10000.0)*(100-10)+10;
    // }

    // return soma;
}

/************************************************************************************
 Method: Decoders
 Description: users need to implement at least one decoder, DecK (K = [1,2,3,4,5])
*************************************************************************************/

void Dec3(TSol &s) // best RNC
{
    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < n)
		    s.vec[j].sol = j;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);   

    // criar as estruturas dos RNCs
    std::vector< std::vector<int> > RNC;
    std::vector<double> capRNC;

    RNC.resize(numClusters);
    capRNC.resize(numClusters);
 

    //empacotar os nodes nos rncs
    int penalty = 0;
    for(int i=0; i<numNodes; i++)
    {
        int bestRNC = -1;
        int bestCost = -1;

        // alocar os primeiros numClusters nodes em RNCs separados
        if(i<numClusters)
        {
            capRNC[i] += ti[s.vec[i].sol];
            RNC[i].push_back(s.vec[i].sol);
        }
        else
        {
            //encontrar o melhor rnc k que possui capacidade para atender o node i
            int k=0;
            while(k < numClusters)
            {
                //verificar a capacidade disponivel do rnc k
                if (capRNC[k] + ti[s.vec[i].sol] <= cap)
                {
                    // calcular o custo de inserir o ponto i no rnc k
                    int cost = 0;
                    for (int j=0; j<(int)RNC[k].size(); j++)
                    {
                        cost += hij[RNC[k][j]][s.vec[i].sol] + hij[s.vec[i].sol][RNC[k][j]];
                    }

                    if (cost > bestCost)
                    {
                        bestRNC = k;
                        bestCost = cost;
                    }
                }
                k++;
            }

            //inserir o node i no melhor rnc aberto disponivel
            if (bestRNC >= 0)
            {
                // atualizar a capacidade
                capRNC[bestRNC] += ti[s.vec[i].sol];

                // inserir a base
                RNC[bestRNC].push_back(s.vec[i].sol);
            }

            // penalizar a solucao se um node nao foi agrupado
            else
            {
                penalty += 1000 + 100*(ti[s.vec[i].sol]+capRNC[bestRNC] - cap); 
            }
        }
    }

    s.ofv = CalculateFitness(RNC) + penalty;
}

void Dec2(TSol &s) // first RNC
{
    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < n)
		    s.vec[j].sol = j;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);   

    // criar as estruturas dos RNCs
    std::vector< std::vector<int> > RNC;
    std::vector<double> capRNC;

    RNC.resize(numClusters);
    capRNC.resize(numClusters,0);
 

    //empacotar os nodes nos rncs
    int penalty = 0;
    for(int i=0; i<numNodes; i++)
    {
        int primeiroRNC = -1;

        //encontrar o primeiro rnc j que possui capacidade para atender o node i
        int j=0;
        while(j < numClusters)
        {
            //verificar a capacidade disponivel do rnc j
            if (capRNC[j] + ti[s.vec[i].sol] <= cap)
            {
                primeiroRNC = j;
                break;
            }
            j++;
        }

        //inserir o node i no primeiro rnc aberto disponivel
        if (primeiroRNC >= 0 && j < numClusters)
        {
            capRNC[primeiroRNC] += ti[s.vec[i].sol];
            RNC[primeiroRNC].push_back(s.vec[i].sol);
        }

        // penalizar a solucao se um node nao foi agrupado
        else
        {
            penalty += 99999999; 
        }
    }

    s.ofv = CalculateFitness(RNC) + penalty;

    // s.ofv = ((double)(rand()%10000)/10000.0)*(1000000-1000)+1000;

    // printf("\nFO: %.0lf\n",s.ofv);
}

void Dec1(TSol &s) // best RNC without m initial alocation
{
    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < n-1)
		    s.vec[j].sol = j;
        else
        if (j == n-1)
            // s.vec[j].sol = (s.vec[n-1].rk * (numClusters)) + 1; // define o numero de bases alocadas em diferentes RNCs
            s.vec[j].sol = (s.vec[n-1].rk * (numClusters/2)) + numClusters/2; 
            // s.vec[j].sol = numClusters;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-2, sortByRk);   

    // criar as estruturas dos RNCs
    std::vector< std::vector<int> > RNC;
    std::vector<double> capRNC;

    RNC.resize(numClusters);
    capRNC.resize(numClusters);

    // definir o numero de RNC usados inicialmente
    int numRNC = s.vec[n-1].sol;

    //empacotar as bases nos rncs
    int penalty = 0;
    for(int i=0; i<numNodes; i++)
    {
        int bestRNC = -1;
        int bestCost = -1;

        // alocar os primeiros numRNC bases em RNCs separados
        if(i<numRNC)
        {
            capRNC[i] += ti[s.vec[i].sol];
            RNC[i].push_back(s.vec[i].sol);
        }
        else
        {
            //encontrar o melhor rnc k que possui capacidade para atender o node i
            int k=0;
            while(k < numClusters)
            {
                //verificar a capacidade disponivel do rnc k
                if (capRNC[k] + ti[s.vec[i].sol] <= cap)
                {
                    // calcular o custo de inserir o ponto i no rnc k
                    int cost = 0;
                    for (int j=0; j<(int)RNC[k].size(); j++)
                    {
                        cost += hij[RNC[k][j]][s.vec[i].sol] + hij[s.vec[i].sol][RNC[k][j]];
                    }

                    if (cost > bestCost)
                    {
                        bestRNC = k;
                        bestCost = cost;
                    }
                }
                k++;
            }

            //inserir o node i no melhor rnc aberto disponivel
            if (bestRNC >= 0)
            {
                // atualiza a capacidade 
                capRNC[bestRNC] += ti[s.vec[i].sol];

                // insere a base 
                RNC[bestRNC].push_back(s.vec[i].sol);
            }

            // penalizar a solucao se um node nao foi agrupado
            else
            {
                penalty += 100000 + ti[s.vec[i].sol] * 10; 
            }
        }
    }

    s.ofv = CalculateFitness(RNC) + penalty;
}
void Dec4(TSol &s){}
void Dec5(TSol &s){}

/************************************************************************************
 Method: Local Search Heuristics
 Description: users need to implement local searchs heuristics if he/she use the RVND,
 LSK (K = [1,2,3,4,5])
*************************************************************************************/

void LS1(TSol &s) 
{
    TSol sViz = s;
    // trocar as random keys de dois nodes
    for (int i = 0; i < numNodes-1; i++)
    {
        for (int j = i+1; j < numNodes; j++)
        {
            // trocar as random keys dos nodes i e j
            double temp = sViz.vec[i].rk;
            sViz.vec[i].rk = sViz.vec[j].rk;
            sViz.vec[j].rk = temp;

            Decoder(sViz);

            // atualiza a solucao corrente se melhorou e continua a busca a partir de sViz
            if (sViz.ofv < s.ofv){
                s = sViz;
            }
            // retorna para a solucao corrente
            else{
                sViz = s;
                // temp = sViz.vec[hi].rk;
                // sViz.vec[hi].rk = sViz.vec[pj].rk;
                // sViz.vec[pj].rk = temp;
            }
        }
    } 
}

void LS2(TSol &s){}
void LS3(TSol &s){}
void LS4(TSol &s){}
void LS5(TSol &s){}

/************************************************************************************
 Method: ImprimirSolHMP
 Description: print the problem solution
*************************************************************************************/
void ImprimirSolHMP(std::vector< std::vector<int> > RNC, std::vector<double> capRNC)
{
    for (int k=0; k<(int)RNC.size(); k++)
    {
        printf("\nRNC %d (%.3lf): ", k, capRNC[k]);
        // total handover between base stations assigned to the same RNC
        for (int i=0; i<(int)RNC[k].size(); i++)
        {
            printf("%d ", RNC[k][i]);
        }
    }
}

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem()
{
    //specific problem
    hij.clear();
    ti.clear();
}

#endif