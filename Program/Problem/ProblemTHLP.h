// *******************************************************************
//      file with specific functions to solve the THLP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Global Variables
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int n;                               // size of cromossoms


//---------------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC --------------------




//---------- DEFINITION OF CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  ------------

static std::vector <std::vector <double> > dist;	// matrix with Euclidean distance
static std::vector <std::vector <double> > flow;	// matrix with flow

static int nPontos;
static int nHubs;
static double alfaH;


//----------------------- IMPLEMENTATION OF FUNCTIONS  -------------------------------

// Sort TSol by sol
bool sortBySol(const TVecSol &lhs, const TVecSol &rhs) { return lhs.sol < rhs.sol; }

// Sort by column
bool sortcol(const std::vector<double>& v1, const std::vector<double>& v2){ return v1[0] < v2[0]; }

/************************************************************************************
 Method: find
 Description: Find set of vertex i
*************************************************************************************/
int find(int i, std::vector<int>& parent)
{
    while (parent[i] != i)
        i = parent[i];
    return i;
}

/************************************************************************************
 Method: union
 Description: does union of i and j. It returns false if i and j are already in same set.
*************************************************************************************/
void union1(int i, int j, std::vector<int>& parent)
{
    int a = find(i,parent);
    int b = find(j,parent);
    parent[a] = b;
}


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
    fscanf(arq, "%d", &nPontos);
    fscanf(arq, "%d", &nHubs);
    fscanf(arq, "%lf", &alfaH);

    n = nPontos + (nPontos - nHubs) + ((nHubs * (nHubs-1))/2);

    // imprimir
    // printf("\nPontos: %d \nHubs: %d \nAlfa: %.1lf \nn: %d", nPontos, nHubs, alfaH, n);
    

    //  euclidean distance
    dist.clear();
    dist.resize(nPontos, std::vector<double>(nPontos));

    flow.clear();
    flow.resize(nPontos, std::vector<double>(nPontos));


    // read node informations
    while (!feof(arq))
    {
        int i, j;
        double d = 0.0, f = 0.0;
    	fscanf(arq, "%d %d %lf %lf", &i, &j, &f, &d);

        flow[i][j] = f;
        dist[i][j] = d/1000.0;
    }
    fclose(arq);


    // imprimir
    // for (int i=0; i<nPontos; i++)
    // {
    //     printf("\n");
    //     for (int j = 0; j < nPontos; j++)
    //     {
    //         printf("%.2lf\t",dist[i][j]);
    //     }
    // }

    // for (int i=0; i<nPontos; i++)
    // {
    //     printf("\n");
    //     for (int j = 0; j < nPontos; j++)
    //     {
    //         printf("%.2lf\t",flow[i][j]);
    //     }
    // }

    // getchar();
}

/************************************************************************************
 Method: Decoders
 Description: users need to implement at least one decoder, DecK (K = [1,2,3,4,5])
*************************************************************************************/

void Dec1(TSol &s) 
{
    std::vector <int> a;							// assigment vector
    std::vector <int> h;							// hub vector
    std::vector <int> hP;						    // hub position vector 
    std::vector <std::vector <int> > Tp;			// matrix with tree
    std::vector <std::vector <double> > cT;		    // matrix with cost of tree

    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < nPontos){
		    s.vec[j].sol = j;                    // parte 1
        }
        else{
            if (j < (nPontos + nPontos - nHubs)) // parte 2 
                s.vec[j].sol = 0;
            
            else 
                if (j < n)                       // parte 3
                    s.vec[j].sol = j - (nPontos + nPontos - nHubs);
                else
                    s.vec[j].sol = -1;
        }
	}

    // imprimir
    // printf("\ns: ");
    // for (int j = 0; j < n+1; j++)
	// {
    //     printf("%d (%.2lf) ",s.vec[j].sol, s.vec[j].rk);
    // }
    

    // inicializar as estruturas
    a.clear();
    a.resize(nPontos,0);

    h.clear();
    h.resize(nHubs,0);

    hP.clear();
    hP.resize(nPontos,-1);


    Tp.clear();
    Tp.resize(nHubs, std::vector<int>(nHubs));

    // sort random-key vector until nPontos
    sort(s.vec.begin(), s.vec.begin()+nPontos, sortByRk);   

    // imprimir
    // printf("\ns: ");
    // for (int j = 0; j < n+1; j++)
	// {
    //     printf("%d (%.2lf) ",s.vec[j].sol, s.vec[j].rk);
    // }

    // definir os hubs
    for (int i = 0; i < nHubs; i++)
    {
        h[i] = s.vec[i].sol;
        hP[s.vec[i].sol] = i;
    }

    // printf("\nh[i]: ");
    // for (int i = 0; i < nHubs; i++)
    // {
    //     printf("%d ", h[i]);
    // }

    // printf("\nhP[i]: ");
    // for (int i = 0; i < nPontos; i++)
    // {
    //     printf("%d ", hP[i]);
    // }

    // alocar os pontos 
    int k = nPontos;
    for (int i = 0; i < nPontos; i++)
    {
        // o ponto eh um hub
        if (i < nHubs)
            a[s.vec[i].sol] = s.vec[i].sol;
        else
        {
            a[s.vec[i].sol] = h[floor(s.vec[k].rk * nHubs)];
            s.vec[k].sol = a[s.vec[i].sol];
            k++;
        }
    }

    // printf("\n\ns: ");
    // for (int j = 0; j < n+1; j++)
	// {
    //     printf("%d (%.2lf) ",s.vec[j].sol, s.vec[j].rk);
    // }

    // printf("\na[i]: ");
    // for (int i = 0; i < nPontos; i++)
    // {
    //     printf("[%d] %d, ", i, a[i]);
    // }
    
    
    // definir a arvore
    std::vector <std::vector <double> > tree;
    tree.clear();
    tree.resize((nHubs * (nHubs-1))/2, std::vector<double>(4));

    k = 0;
    int kI = nPontos + nPontos - nHubs;
    for (int i=0; i<nHubs-1; i++)
    {
        for (int j=i+1; j<nHubs; j++)
        {
            tree[k][0] = s.vec[kI+k].rk;
            tree[k][1] = h[i];
            tree[k][2] = h[j];
            tree[k][3] = k;
            k++;
        }
    }

    // printf("\n\ntree: ");
    // for (int i = 0; i < (nHubs * (nHubs-1))/2; i++)
    // {
    //     printf("\n%lf \t %.0lf \t%.0lf", tree[i][0], tree[i][1], tree[i][2]);
    // }

    // ordenar as chaves aleatorias dos arcos dos hubs
    std::sort(tree.begin(), tree.end(), sortcol);

    // gravar a informacao em s.vec.sol
    for (std::vector<int>::size_type i=0; i<tree.size(); i++)
    {
        s.vec[kI+i].sol = tree[i][3];
    }


    // Kruskal
    // Initialize sets of disjoint sets.
    std::vector <int> parent;
    parent.clear();
    parent.resize(nHubs);
    for (int i = 0; i < nHubs; i++){
        parent[i] = i;
    }

    // Include edges one by one
    int edge_count = 0;
    k = 0;
    while (edge_count < nHubs - 1) {
        int a = hP[tree[k][1]], 
            b = hP[tree[k][2]];

        // se o arco k nao forma ciclo, posso inserir 1 em Tp[a][b]
        if (find(a,parent) != find(b,parent))
        {
            union1(a,b,parent);
            Tp[a][b] = Tp[b][a] = 1;
            edge_count++;
        }

        // vai para o proximo arco
        k++;
    }


    // ******** poderia ter usado a mesma matriz para Tp e cT



    // printf("\n\ntree: ");
    // for (int i = 0; i < (nHubs * (nHubs-1))/2; i++)
    // {
    //     printf("\n%lf \t %.0lf \t%.0lf", tree[i][0], tree[i][1], tree[i][2]);
    // }

    // printf("\n\nTp: ");
    // for (int i = 0; i < nHubs; i++)
    // {
    //     printf("\n");
    //     for (int j = 0; j < nHubs; j++)
    //     {
    //         printf("%d \t", Tp[i][j]);
    //     }
    // }
    // getchar();

    // Cria uma matriz cT[V][V] que armazena as distâncias mais curtas entre todos os pares de hubs
    cT.clear();
    cT.resize(nHubs, std::vector<double>(nHubs));
    for (int i = 0; i < nHubs; i++) {
        for (int j = 0; j < nHubs; j++) {
            if (i != j){
                if (Tp[i][j] > 0)
                    cT[i][j] = dist[h[i]][h[j]];
                else
                    cT[i][j] = 100000000;
            }
            else
                cT[i][j] = 0;
        }
    }

    // Executa o algoritmo de Floyd
    for (int k = 0; k < nHubs; k++) { // Escolhe um vértice intermediário
        for (int i = 0; i < nHubs; i++) { // Escolhe o vértice de origem
            for (int j = 0; j < nHubs; j++) { // Escolhe o vértice de destino
                // Se a distância de i a j passando pelo vértice intermediário k for menor que a distância atual, atualiza a distância
                if (cT[i][k] + cT[k][j] < cT[i][j]) { 
                    cT[i][j] = cT[i][k] + cT[k][j];
                }
            }
        }
    }

    // printf("\n\ncT: ");
    // for (int i = 0; i < nHubs; i++)
    // {
    //     printf("\n");
    //     for (int j = 0; j < nHubs; j++)
    //     {
    //         printf("%.2lf \t", cT[i][j]);
    //     }
    // }
    // getchar();


    // Calcular a função objetivo
    double custo = 0;
    for (int i=0; i<nPontos; i++){
        for (int j=0; j<nPontos; j++){
                //encontrar o custo do caminho entre i e j
                double firstMile = dist[i][a[i]];	                  // custo de i até seu hub
                double lastMile = dist[a[j]][j]; 	                  // custo do hub de j até j

                double middleMile = alfaH * cT[hP[a[i]]][hP[a[j]]];   // custo entre hubs de i e j

                custo += (firstMile + lastMile + middleMile) * flow[i][j];
        }
    }

    // s.ofv = CalculateFitness(s);
    s.ofv = custo;
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

void LS1(TSol &s) // trocar hubs
{
    TSol sViz = s;
    for (int i=0; i<nHubs; i++)
    {
        for (int j =nHubs; j<nPontos; j++)
        {
            // trocar um hub hi por um ponto pj nao hub 
            int hi = sViz.vec[i].sol;
            int pj = sViz.vec[j].sol;

            double temp = sViz.vec[hi].rk;
            sViz.vec[hi].rk = sViz.vec[pj].rk;
            sViz.vec[pj].rk = temp;

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

void LS2(TSol &s) // trocar atribuicao dos pontos
{
    // definir o fator
    float fator = 1.0/nHubs;

    // 1o hub, i=0 => [i*fator, (i+1)*fator[ = [0.00, 0.33[
    // 2o hub, i=1 => [i*fator, (i+1)*fator[ = [0.33, 0.66[
    // 3o hub, i=2 => [i*fator, (i+1)*fator[ = [0.66, 0.99[

    TSol sViz = s;
    // double rk_anterior;
 
    // trocar a alocacao de cada ponto nao hub
    for (int i=nPontos; i<nPontos+nPontos-nHubs; i++)
    {
        // int hubAtual = hP[sViz.vec[i].sol];

        // escolher uma nova chave para o ponto
        for (int j=0; j<nHubs; j++)
        {
            // if (j != hubAtual)
            // {
                // gerar uma chave no intervalo j
                double min = j*fator;
                double max = (j+1)*fator;
                // rk_anterior = sViz.vec[i].rk;
                sViz.vec[i].rk = ((double)(rand()%10000)/10000.0)*(max-min)+min;
            // }

            Decoder(sViz);

            if (sViz.ofv < s.ofv){
                s = sViz;
                // hubAtual = hP[s.vec[i].sol];
            }
            else{
                sViz = s;
                // sViz.vec[i].rk = rk_anterior;
            }
        }
    }
}

void LS3(TSol &s) // trocar arvore
{
    TSol sViz = s;
    int inicio = nPontos+nPontos-nHubs;
    int fim = s.vec.size()-1;

    for (int i=inicio; i<inicio+nHubs-1; i++)
    {
        for (int j =inicio+nHubs-1; j<fim; j++)
        {
            // trocar uma arco i da arvore por um arco j nao selecionado (decoder verifica se nao forma ciclo)
            int ai = sViz.vec[i].sol;
            int aj = sViz.vec[j].sol;

            double temp = sViz.vec[ai].rk;
            sViz.vec[ai].rk = sViz.vec[aj].rk;
            sViz.vec[aj].rk = temp;

            Decoder(sViz);

            // atualiza a solucao corrente se melhorou e continua a busca a partir de sViz
            if (sViz.ofv < s.ofv){
                s = sViz;
            }
            // retorna para a solucao corrente
            else{
                sViz = s;
                // temp = sViz.vec[ai].rk;
                // sViz.vec[ai].rk = sViz.vec[aj].rk;
                // sViz.vec[aj].rk = temp;
            }
        }
    }  
}

void LS4(TSol &s){}
void LS5(TSol &s){}

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem()
{
    //specific problem
    dist.clear();
    flow.clear();
}

#endif
