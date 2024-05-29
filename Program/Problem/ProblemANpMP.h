// *******************************************************************
//      file with specific functions to solve the alpha-NpMP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Variables declared in main.cpp
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int n;                               // size of cromossoms

//------ DEFINITION OF TYPES OF PROBLEM SPECIFIC --------


//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------

static std::vector <std::vector <int> > dist;	// matrix with Euclidean distance

static int numNodes;								// number of nodes
static int numMedians;								// number of medians
static int alphaN;									// number of alpha neighbors


//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------

/************************************************************************************
 Method: ReadData
 Description: read the input data
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
    int aux;
    fscanf(arq, "%d %d %d", &numNodes, &aux, &numMedians);

    numMedians = 10;
    alphaN = 5;

    dist.clear();
    dist.resize(numNodes, std::vector<int>(numNodes, 9999999));
    for (int i = 0; i < numNodes; i++){
        dist[i][i] = 0;
    }

    // read distance informations
    while (!feof(arq))
    {
        int i, j, distance;
    	fscanf(arq, "%d %d %d", &i, &j, &distance);
    	dist[i-1][j-1] = distance;
        dist[j-1][i-1] = distance;

        // printf("\n%d \t %d \t %d", i, j, dist[i-1][j-1]);
    }
    fclose(arq);

    // Floyd-Warshall algorithm to compute the shortest paths between every pair of vertices
    /* Add all vertices one by one to the set of intermediate vertices.
      ---> Before start of a iteration, we have shortest distances between all
      pairs of vertices such that the shortest distances consider only the
      vertices in set {0, 1, 2, .. k-1} as intermediate vertices.
      ----> After the end of a iteration, vertex no. k is added to the set of
      intermediate vertices and the set becomes {0, 1, 2, .. k} */

    for (int k = 0; k < numNodes; k++){
        // Pick all vertices as source one by one
        for (int i = 0; i < numNodes; i++){
            // Pick all vertices as destination for the above picked source
            for (int j = 0; j < numNodes; j++){
                // If vertex k is on the shortest path from i to j, then update the value of dist[i][j]
                if (dist[i][k] + dist[k][j] < dist[i][j])
                    dist[i][j] = dist[i][k] + dist[k][j];
            }
        }
    }

    // print
    // for (int i = 0; i < numNodes; i++){
    //     printf("\n");
    //     for (int j = 0; j < numNodes; j++){
    //         printf("%d ", dist[i][j]);
    //     }
    // }
    // getchar();

    

    // define the size of the solution vector
    n = numMedians;
}

/************************************************************************************
 Method: CalculateFitness
 Description: calculate the fitness of a chromossom s
*************************************************************************************/
double CalculateFitness(TSol s)
{    
    s.ofv = 0;
    // std::vector <std::vector <int> > solution;
    // solution.resize(numNodes, std::vector<int>(alphaN));

    // for each node find the alpha nearest medians 
    for (int j=0; j<numNodes; j++)
    {
        std::vector<int> medians(numMedians,0);
        int minDist;
        int minIndex;

        // find the k nearest median
        for (int k=0; k<alphaN; k++)
        {
            minDist = 99999999;
            minIndex = 0;

            for (int i=0; i<numMedians; i++)
            {
                if ( (dist[s.vec[i].sol][j] < minDist) && (medians[i] == 0) )
                {
                    minDist = dist[s.vec[i].sol][j];
                    minIndex = i;
                }
            }

            // insert the alpha nearest median in the solution
            // solution[j][k] = s.vec[minIndex].sol;
            medians[minIndex] = 1;

            s.ofv += minDist;

            // printf("\n");
            // for (int i=0; i<numMedians; i++){
            //     printf("%d \t", medians[i]);
            // }
        }
        // getchar();
    }

    // printf("\n");
    // for (int i=0; i<numMedians; i++)
    // {
    //     printf("%d \t", s.vec[i].sol);
    // }
    // printf("\n");

    // for (int j=0; j<numNodes; j++)
    // {
    //     printf("\n");
    //     for (int k=0; k<alphaN; k++)
    //     {
    //         printf("%d \t", solution[j][k]);
    //     }
    // }
    // getchar();


    return s.ofv;
}

/************************************************************************************
 Method: Decoders
 Description: users need to implement at least one decoder, DecK (K = [1,2,3,4,5])
*************************************************************************************/
void Dec1(TSol &s) // mapping
{
    // create an initial solution of the problem
    std::vector<int> nodesList;

    // list of candidate nodes
    for (int j = 0; j < numNodes; j++){
        nodesList.push_back(j);
	}

    // define the p medians
    for (int i=0; i<n; i++)
    {
        // define the index of the nodeList
        int index = floor(s.vec[i].rk * (int)nodesList.size());

        // define the node that is a median
        s.vec[i].sol = nodesList[index];

        // remove this node of the nodeList
        nodesList.erase(nodesList.begin()+index);
    }

    s.ofv = CalculateFitness(s);
}

void Dec2(TSol &s){}
void Dec3(TSol &s){}
void Dec4(TSol &s){}
void Dec5(TSol &s){}

/************************************************************************************
 Method: Local Search Heuristics
 Description: users need to implement local searchs heuristics if he/she use the RVND,
 LSK (K = [1,2,3,4,5])
*************************************************************************************/
/*void LS1(TSol &s){
    // TSol sViz = s;

    // // inicializar o vetor de nos
    // std::vector<double> nodes;
    // nodes.resize(numNodes);
    // for (int i = 0; i < numNodes; i++)
    // {
    //     nodes[i] = i;
    // }
    
    // std::random_shuffle (nodes.begin(), nodes.end());

    // int melhorou = 1;

    // while (melhorou)
    // {
    //     melhorou = 0;
    
    //     // para todos os nos
    //     for (int i=0; i<numNodes; i++)
    //     {
    //         // escolher uma mediana para sair
    //         int j = irandomico(0, s.vec.size()-1);

    //         // trocar a chave j por um valor aleatorio no intervalo 
    //         sViz.vec[j].rk = randomico(nodes[i]*1.0/numNodes, (nodes[i]+1)*1.0/numNodes);

    //         Decoder(sViz);

    //         // atualiza a solucao corrente se melhorou e continua a busca a partir de sViz
    //         if (sViz.ofv < s.ofv){
    //             s = sViz;
    //             melhorou = 1;
    //         }
    //         // retorna para a solucao corrente
    //         else{
    //             sViz = s;
    //         }
    //     }  
    // }
}
void LS2(TSol &s){
    // TSol sViz = s;

    // // para cada mediana
    // for (int i=0; i<n; i++)
    // {
    //     // testar as numNodes faixas de valores
    //     for (int j = 0; j<numNodes; j++)
    //     {
    //         // trocar a chave i por um valor aleatorio no intervalo [j*1/numNodes, (j+1)*1/numNodes]
    //         sViz.vec[i].rk = randomico(j*1.0/numNodes, (j+1) * (1.0/numNodes));

    //         Decoder(sViz);

    //         // atualiza a solucao corrente se melhorou e continua a busca a partir de sViz
    //         if (sViz.ofv < s.ofv){
    //             s = sViz;
    //         }
    //         // retorna para a solucao corrente
    //         else{
    //             sViz = s;
    //         }
    //     }
    // }  
}
void LS3(TSol &s){}
void LS4(TSol &s){}
void LS5(TSol &s){}
*/

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/

void FreeMemoryProblem()
{
    //specific problem
    dist.clear();
}

#endif